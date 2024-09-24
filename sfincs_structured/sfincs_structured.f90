program sfincs_structured
   !
   ! use omp_lib
   !
   implicit none
   !
   integer :: num_args, ix
   character(len=12), dimension(:), allocatable :: args
   !
   integer :: n, m, nmax, mmax, it, nt
   integer*8  :: count0, count1, count_rate, count_max
   !
   real*4 :: g, zbunif, manning, dt, dx, dy, gnavg2, hu, hv, frc
   real*8 :: ttotal
   !
   logical :: copyq, momentum, continuity
   !
   integer, dimension(:,:), allocatable :: kcs
   integer, dimension(:,:), allocatable :: kcu
   integer, dimension(:,:), allocatable :: kcv
   !
   real*4, dimension(:,:), allocatable :: zs
   real*4, dimension(:,:), allocatable :: zb
   real*4, dimension(:,:), allocatable :: zbu
   real*4, dimension(:,:), allocatable :: zbv
   real*4, dimension(:,:), allocatable :: qu
   real*4, dimension(:,:), allocatable :: qv
   real*4, dimension(:,:), allocatable :: qu0
   real*4, dimension(:,:), allocatable :: qv0
   !
   real*4, parameter :: expo = 1.0 / 3.0
   !integer, parameter :: expo = 0   
   !
   !call acc_init( acc_device_nvidia )
   !
   ! Read input arguments
   !
   num_args = command_argument_count()
   allocate(args(num_args))
   do ix = 1, num_args
      call get_command_argument(ix,args(ix))
   end do
   !
   read(args(1),*)nmax
   read(args(2),*)mmax
   read(args(3),*)nt
   read(args(4),*)copyq
   read(args(5),*)momentum
   read(args(6),*)continuity
   !
   zbunif  = - 10.0
   manning = 0.02
   dx      = 10.0
   dy      = 10.0
   !
   g       = 9.81
   !
   dt = dx / sqrt(- g * zbunif)
   !
   gnavg2 = g * manning**2
   !
   allocate(kcs(nmax, mmax))
   allocate(kcu(nmax, mmax))
   allocate(kcv(nmax, mmax))
   allocate(zs(nmax, mmax))
   allocate(zbu(nmax, mmax))
   allocate(zbv(nmax, mmax))
   allocate(qu(nmax, mmax))
   allocate(qv(nmax, mmax))
   allocate(qu0(nmax, mmax))
   allocate(qv0(nmax, mmax))
   !
   kcs = 1
   kcu = 0
   kcv = 0
   zs = 0.0
   qu = 0.0
   qv = 0.0   
   qu0 = 0.0
   qv0 = 0.0   
   ! zb = zbunif
   zbu = zbunif
   zbv = zbunif
   !
   ! Set mask values
   !
   kcs(1:nmax, 1) = 0
   kcs(1:nmax, mmax) = 0
   kcs(1, 1:mmax) = 0
   kcs(nmax, 1:mmax) = 0
   !
   ! Compute bed level at u and v points and set mask
   !
   do m = 1, mmax - 1
      do n = 1, nmax - 1
         if (kcs(n, m) == 1 .and. kcs(n, m + 1)) then
            kcu(n, m) = 1		 
		   endif		   
	      if (kcs(n, m) == 1 .and. kcs(n + 1, m)) then
		      kcv(n, m) = 1		 
		   endif		   
      enddo
   enddo  	  
   !
   !$acc data, copyin( qu0, qv0, qu, qv, kcs, kcu, kcv, zs, zbu, zbv )
   !
   call system_clock(count0, count_rate, count_max)
   !
   do it = 1, nt   
      !
      !$acc parallel, present( qu0, qv0, qu, qv, kcs, kcu, kcv, zs, zbu, zbv ), num_gangs( 512 ), vector_length( 128 )
      !
      ! Momentum
      !
      if (copyq) then
      !$omp parallel &
      !$omp private ( n, m )
      !$omp do
      !$acc loop independent gang
      do m = 1, mmax
         !$acc loop independent vector
         do n = 1, nmax
            qu0(n, m) = qu(n, m)
            qv0(n, m) = qv(n, m)
         enddo		 
      enddo
      !$omp end do
      !$omp end parallel
      endif
      !
      if (momentum) then
      !$omp parallel &
      !$omp private ( n, m, hu, hv, frc )
      !$omp do
      !$acc loop independent gang
      do m = 1, mmax 
         !$acc loop independent vector
         do n = 1, nmax
            !
	         ! U
            ! 
            if (kcu(n, m) == 1) then
               !
               hu = 0.5 * ( zs(n, m + 1) + zs(n, m) ) - zbu(n, m)
               frc = - g * hu * (zs(n, m + 1) - zs(n, m)) / dx
               qu(n, m) = (qu0(n, m) + frc * dt) / (1.0 + gnavg2 * dt * abs(qu0(n, m)) / (hu**2 * hu**expo))
               !
            endif
	         !
            ! V
            !  
            if (kcv(n, m) == 1) then
               !
               hv = 0.5 * ( zs(n + 1, m) + zs(n, m) ) - zbv(n, m)
               frc = - g * hv * (zs(n + 1, m) - zs(n, m)) / dy
               qv(n, m) = (qv0(n, m) + frc * dt) / (1.0 + gnavg2 * dt * abs(qv0(n, m)) / (hv**2 * hv**expo))
               !
            endif
            !
         enddo		 
      enddo
      !$omp end do
      !$omp end parallel
      endif
      !
      ! Continuity
      !
      if (continuity) then
      !$omp parallel &
      !$omp private ( n, m )
      !$omp do
      !$acc loop independent gang
      do m = 1, mmax
         !$acc loop independent vector
         do n = 1, nmax
	         !
            if (kcs(n, m) == 1) then
               !
               zs(n, m) = zs(n, m) + dt * ( (qu(n, m - 1) - qu(n, m)) / dx + (qv(n - 1, m) - qv(n, m)) / dy ) 
               !
	 	      endif
            !
         enddo		 
      enddo    
      !$omp end do
      !$omp end parallel
      endif
      !
      !$acc end parallel
      !
   enddo
   !
   !$acc end data
   !
   call system_clock(count1, count_rate, count_max)
   !
   ttotal = 1.0 * (count1 - count0) / count_rate
   !
   write(*,'(f10.4)') ttotal
   !
   open(500, file='duration.txt')
   write(500,'(f10.4,a)')ttotal
   close(500)
   !
end program
