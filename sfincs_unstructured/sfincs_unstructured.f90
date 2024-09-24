program sfincs_unstructured
   !
   ! use omp_lib
   !
   implicit none
   !
   integer :: num_args, ix
   character(len=12), dimension(:), allocatable :: args
   !
   integer    :: n, m, nm, ip, nmu, num, nmax, mmax, it, nt, np, npuv
   integer*8  :: count0, count1, count_rate, count_max
   !
   real*4 :: g, zbunif, manning, dt, dx, dy, gnavg2, huv, frc, alfa, zsuv
   real*8 :: ttotal
   !
   logical :: copyq, momentum, continuity
   !
   integer, dimension(:,:), allocatable :: kcs0
   integer, dimension(:,:), allocatable :: nms
   integer, dimension(:), allocatable   :: kcs
   integer, dimension(:), allocatable   :: kcuv
   !
   integer, dimension(:), allocatable   :: iunm
   integer, dimension(:), allocatable   :: iunmu
   integer, dimension(:,:), allocatable   :: iunb
   !
   integer, dimension(:), allocatable   :: izmd
   integer, dimension(:), allocatable   :: izmu
   integer, dimension(:), allocatable   :: iznd
   integer, dimension(:), allocatable   :: iznu
   !
   real*4, dimension(:), allocatable :: zs
   real*4, dimension(:), allocatable :: zbuv
   real*4, dimension(:), allocatable :: quv
   real*4, dimension(:), allocatable :: quv0
   !
   integer, dimension(:,:), allocatable :: uv_index_uv   
   integer, dimension(:,:), allocatable :: uv_index_z  
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
   !copyq = .true.
   !momentum = .true.
   !continuity = .true.
   !  
   zbunif  = - 10.0
   manning = 0.02
   dx      = 10.0
   dy      = 10.0
   !
   g       = 9.81
   alfa    = 0.5
   !
   dt = alfa * dx / sqrt(- g * zbunif)
   !
   gnavg2 = g * manning**2
   !
   allocate(kcs0(nmax, mmax))
   allocate(nms(nmax, mmax))
   !
   kcs0 = 1
   !
   ! Set mask values
   !
   kcs0(1:nmax, 1) = 0
   kcs0(1:nmax, mmax) = 0
   kcs0(1, 1:mmax) = 0
   kcs0(nmax, 1:mmax) = 0
   !
   ! Turn into unstructured
   ! Count number of z and uv points
   !
   np = 0
   npuv = 0
   !
   do m = 1, mmax
      do n = 1, nmax
	      if (kcs0(n, m) == 1) then
            np = np + 1
            nms(n, m) = np
   	      if (kcs0(n, m + 1)) then
               npuv = npuv + 1
		      endif		   
    	      if (kcs0(n + 1, m)) then
               npuv = npuv + 1
            endif		   
         endif		   
      enddo
   enddo  	  
   !
   allocate(zs(np))
   allocate(kcs(np))
   allocate(kcuv(npuv))
   allocate(quv(npuv + 1))
   allocate(quv0(npuv + 1))
   allocate(zbuv(npuv))
   !
   allocate(iunm(npuv))
   allocate(iunmu(npuv))
   allocate(iunb(2,npuv))
   allocate(uv_index_uv(8,npuv))
   !
   allocate(izmd(np))
   allocate(izmu(np))
   allocate(iznd(np))
   allocate(iznu(np))
   !
   zs = 0.0
   quv = 0.0
   quv0 = 0.0
   zbuv = zbunif
   kcs = 0
   kcuv = 0
   !
   iunm  = 0
   iunmu = 0
   iunb  = 0
   uv_index_z  = 0
   uv_index_uv = npuv + 1
   izmd  = npuv + 1
   izmu  = npuv + 1
   iznd  = npuv + 1
   iznu  = npuv + 1
   !
   ! Compute bed level at u and v points and set mask
   !
   nm = 0
   ip = 0
   !
   do m = 1, mmax
      do n = 1, nmax
         !
	      if (kcs0(n, m) == 1) then
            !
            nm = nm + 1
            kcs(nm) = 1           
            !
	         if (kcs0(n, m + 1)) then
               !
               ip = ip + 1
               !
               kcuv(ip) = 1
               !
               iunm(ip) = nm
               iunb(1, ip) = nm
               nmu = nms(n, m + 1)
               iunmu(ip) = nmu
               iunb(2, ip) = nmu
               !
               izmu(nm) = ip
               izmd(nmu) = ip
               !
               !
            endif		   
            !
            if (kcs0(n + 1, m)) then
               !
               ip = ip + 1
               !
               kcuv(ip) = 1
               !
               iunm(ip) = nm
               iunb(1, ip) = nm
               num = nms(n + 1, m)
               iunmu(ip) = num
               iunb(2,ip) = num
               !
               iznu(nm) = ip
               iznd(num) = ip
               !
            endif
            !
         endif		   
         !
      enddo
   enddo
   !
   do ip = 1, npuv
      !
      nm  = iunm(ip)
      nmu = iunmu(ip)
      !
      uv_index_uv(1, ip) = izmd(nm)  ! unmd
      uv_index_uv(2, ip) = izmu(nmu) ! unmu
      uv_index_uv(3, ip) = izmu(nmu) ! undm
      uv_index_uv(4, ip) = izmu(nmu) ! unum
      !
      !
   enddo   
   !
   deallocate(kcs0)
   deallocate(nms)
   !
   !$acc data, copyin( quv0, quv, kcs, kcuv, zs, zbuv, iunb, izmd, izmu, iznd, iznu ) 
   !
   call system_clock(count0, count_rate, count_max)
   !   
   do it = 1, nt   
      !
      !$acc parallel, present( quv0, quv, kcs, kcuv, zs, zbuv, iunb, izmd, izmu, iznd, iznu ), num_gangs( 512 ), vector_length( 128 )
      !
      ! Momentum
      !
      if (copyq) then
      !$omp parallel &
      !$omp private ( ip )
      !$omp do
      !$acc loop independent, gang, vector
      do ip = 1, npuv
         quv0(ip) = quv(ip)
      enddo
      !$omp end do
      !$omp end parallel
      endif
      !
      if (momentum) then
      !$omp parallel &
      !$omp private ( ip, nm, nmu, huv, frc )
      !$omp do
      !$acc loop independent, gang, vector
      do ip = 1, npuv
         !
         if (kcuv(ip) == 1) then
            !
            nm  = iunb(1, ip)
            nmu = iunb(2, ip)
            !
            huv = 0.5 * ( zs(nmu) + zs(nm) ) - zbuv(ip)
            frc = - g * huv * (zs(nmu) - zs(nm)) / dx
            quv(ip) = (quv0(ip) + frc * dt) / (1.0 + gnavg2 * dt * abs(quv0(ip)) / (huv**2 * huv**expo))
            !
         endif
         !
      enddo
      !$omp end do
      !$omp end parallel
      endif
      !
      ! Continuity
      !
      if (continuity) then
      !$omp parallel &
      !$omp private ( nm )
      !$omp do
      !$acc loop independent, gang, vector
      do nm = 1, np
         !
         if (kcs(nm) == 1) then
            !
            zs(nm) = zs(nm) + dt * ( (quv(izmd(nm)) - quv(izmu(nm))) / dx + (quv(iznd(nm)) - quv(iznu(nm))) / dy ) 
            !
         endif
         !
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
