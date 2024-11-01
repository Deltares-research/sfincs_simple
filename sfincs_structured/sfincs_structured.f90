module constants
   real*4, device, parameter  :: g = 9.81
   real*4, device, parameter  :: zbunif = - 10.0
   real*4, device, parameter  :: manning = 0.02
   real*4, device, parameter  :: dx = 10.0
   real*4, device, parameter  :: dy = 10.0
   real*4, device, parameter  :: gnavg2 = g * manning * manning
   real*4, device, parameter  :: expo = 1.0 / 3.0
   real*4, device, parameter  :: dt =  dx / sqrt(- g * zbunif)

end module
   
   module mathOps
   
contains

  attributes(global) subroutine saxpy(x, y, a)
    implicit none
    real :: x(:), y(:)
    real, value :: a
    integer :: i, n
    n = size(x)
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    if (i <= n) y(i) = y(i) + a*x(i)
  end subroutine saxpy 
  
    attributes(global) subroutine U(nmax, mmax, dt, zs, zbu, kcu, qu0, qu)
    use constants
    implicit none
    integer, value, intent(in) :: nmax, mmax
    real, dimension(nmax,mmax), intent(in)  :: zs, zbu, qu0, kcu
    real, dimension(nmax,mmax), intent(out) :: qu
    real, dimension(nmax,mmax), shared :: zs_shared
    real, value, intent(in) :: dt
    integer :: n, m

    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    m = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    
   if (n <= nmax && m <= mmax) then
    zs_shared(n,m) = zs(n,m) !> global to shared copy
    
    call syncthreads() 
    
    hu = 0.5 * ( zs_shared(n, m + 1) + zs_shared(n, m) ) - zbu(n, m)
    frc = - g * hu * (zs_shared(n, m + 1) - zs_shared(n, m)) / dx
    qu(n, m) = kcu(n, m) * (qu0(n, m) + frc * dt) / (1.0 + gnavg2 * dt * abs(qu0(n, m)) / (hu*hu * hu**expo))
   end if
  end subroutine U
  
   attributes(global) subroutine V(nmax, mmax, dt, zs, zbv, kcv, qv0, qv)
    use constants
    implicit none
    integer, value, intent(in) :: nmax, mmax
    real, dimension(nmax,mmax), intent(in)  :: zs, zbv, kcv, qv0
    real, dimension(nmax,mmax), intent(out) :: qv
    real, dimension(nmax,mmax), shared :: zs_shared
    real, value, intent(in) :: dt
    integer :: n, m
    
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    m = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    
    if (n <= nmax && m <= mmax) then
      zs_shared(n,m) = zs(n,m) !> global to shared copy
      
      call syncthreads() 
      
      hv = 0.5 * ( zs_shared(n + 1, m) + zs_shared(n, m) ) - zbv(n, m)
      frc = - g * hv * (zs_shared(n + 1, m) - zs_shared(n, m)) / dy
      qv(n, m) = kcv(n, m) * (qv0(n, m) + frc * dt) / (1.0 + gnavg2  * dt * abs(qv0(n, m)) / (hv**2 * hv**expo))
    end if
   
  end subroutine V
  
     attributes(global) subroutine Z(nmax, mmax, dt, zs, qu, qv, kcs)
    use constants
    implicit none
    integer, value, intent(in) :: nmax, mmax
    real, dimension(nmax,mmax), intent(in)  :: qu, qv
    real, dimension(nmax,mmax), intent(inout) :: zs
    real, dimension(nmax,mmax), shared :: qu_shared, qv_shared
    real, value, intent(in) :: dt
    integer :: n, m
    
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    m = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    
    if (n <= nmax && m <= mmax) then
      qu_shared(n,m) = qu(n,m) !> global to shared copy
      qv_shared(n,m) = qv(n,m) !> global to shared copy
      
      call syncthreads() 
      
      zs(n, m) = zs(n, m) + kcs(n, m)* dt * ( (qu_shared(n, m - 1) - qu_shared(n, m)) / dx + (qv_shared(n - 1, m) - qv_shared(n, m)) / dy ) 
    end if
end subroutine Z
   
end module mathOps

   
   program sfincs_structured

   use cudafor
   use mathOps
   
   implicit none

   integer :: num_args, ix
   character(len=12), dimension(:), allocatable :: args
   integer :: n, m, nmax, mmax, it, nt
   integer*8  :: count0, count1, count_rate, count_max
   type(dim3) :: grid, tBlock
   real*8 :: ttotal
   logical :: copyq, momentum, continuity

   integer, device, dimension(:,:), allocatable :: kcs_d
   integer, device, dimension(:,:), allocatable :: kcu_d
   integer, device, dimension(:,:), allocatable :: kcv_d
   integer, dimension(:,:), allocatable :: kcs
   integer, dimension(:,:), allocatable :: kcu
   integer, dimension(:,:), allocatable :: kcv

   real*4, device, dimension(:,:), allocatable :: zs
   real*4, device, dimension(:,:), allocatable :: zb
   real*4, device, dimension(:,:), allocatable :: zbu
   real*4, device, dimension(:,:), allocatable :: zbv
   real*4, device, dimension(:,:), allocatable :: qu
   real*4, device, dimension(:,:), allocatable :: qv
   real*4, device, dimension(:,:), allocatable :: qu0
   real*4, device, dimension(:,:), allocatable :: qv0

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
   zbu = zbunif
   zbv = zbunif

   ! Set mask values
   kcs(1:nmax, 1) = 0
   kcs(1:nmax, mmax) = 0
   kcs(1, 1:mmax) = 0
   kcs(nmax, 1:mmax) = 0
   !
   ! Compute bed level at u and v points and set mask
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
   
   kcs_d = kcs
   kcu_d = kcu
   kcv_d = kcv
   
   call system_clock(count0, count_rate, count_max)
   
   tBlock = dim3(32,32,1)
   grid = dim3(ceiling(real(nmax)/tBlock%x),ceiling(real(mmax)/tBlock%x),1)
  
   do it = 1, nt   
      if (copyq) then
         qu0 = qu
         qv0 = qv
      endif
      if (momentum) then
       call U<<<grid, tBlock>>>(nmax, mmax, dt, zs, zbu, kcu_d, qu0, qu)
       call V<<<grid, tBlock>>>(nmax, mmax, dt, zs, zbv, kcv_d, qv0, qv)
      endif
      if (continuity) then
        call Z<<<grid, tBlock>>>(nmax, mmax, dt, zs, qu, qv, kcs_d)
      endif
   enddo
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
