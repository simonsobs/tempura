!////////////////////////////////////////////////////!
! * Noise spectrum of (real) quadratic estimator
!////////////////////////////////////////////////////!

module noise_spec
  use alkernel, only: kernels_lens, kernels_tau, kernels_rot, kernels_lenstau, get_lfac, clbb_est
  implicit none

  private kernels_lens, kernels_tau, kernels_rot, kernels_lenstau, get_lfac, clbb_est

contains


subroutine qtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,Nl,gtype,temp_arg)
!*  Noise spectrum of reconstructed fields from the temperature quadratic estimator (asymmetric case)
!*
!*  Args:
!*    :lmax (int)          : Maximum multipole of output noise spectrum
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipole for reconstruction
!*    :wx0/1 [l] (double)  : 1/(CXX+NXX) (0:rlmax)
!*    :wxy0/1 [l] (double) : CXY/(CXX+NXX) (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Multiplying square of L(L+1)/2, i.e., convergence ('k') or lensing potential ('', default)
!*
!*  Returns:
!*    :Nl [2,l] (double) : Noise spectrum (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: gtype
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: wx0, wx1, wxy0, wxy1, a0a1, b0b1, a0b1, a1b0
  double precision, intent(out), dimension(2,0:lmax) :: Nl
  !internal
  integer :: rL(2), l
  double precision, dimension(lmax) :: lk2
  double precision, dimension(:,:), allocatable :: W1, W2
  double precision, dimension(2,2,lmax) :: SG
  !opt4py :: gtype = ''

  rL = (/rlmin,rlmax/)

  !gradient-leg and C-inverse leg
  allocate(W1(2,rL(1):rL(2)),W2(2,rL(1):rL(2))); W1=0d0; W2=0d0

  do l = rlmin, rlmax
    W1(1,l) = wx0(l)  * wx1(l)  * a0a1(l)
    W2(1,l) = wxy0(l) * wxy1(l) * b0b1(l)
    W1(2,l) = wx0(l)  * wxy1(l) * a0b1(l)
    W2(2,l) = wxy0(l) * wx1(l)  * a1b0(l)
  end do

  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,gtype,lk2)
    call Kernels_lens(rL,W1(1,:),W2(1,:),SG(1,:,:),'S0')
    call Kernels_lens(rL,W1(2,:),W2(2,:),SG(2,:,:),'G0')
  case('amp','src')
    call get_lfac(lmax,gtype,lk2)
    call Kernels_tau(rL,W1(1,:),W2(1,:),SG(1,1,:),'S0')
    call Kernels_tau(rL,W1(2,:),W2(2,:),SG(2,1,:),'G0')
  end select

  Nl = 0d0
  do l = 1, lmax
    if (lk2(l)/=0d0)  Nl(1,l) = sum(SG(:,1,l))/lk2(l)
    if (lk2(l)/=0d0)  Nl(2,l) = sum(SG(:,2,l))/lk2(l)
  end do
  Nl(2,1) = 0d0

end subroutine qtt_asym


subroutine xtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,Nl,gtype,temp_arg)
!*  Noise spectrum of reconstructed fields from the temperature quadratic estimator (asymmetric case)
!*
!*  Args:
!*    :est (str)           : Estimator combination (lensamp,lenssrc,amplens,ampsrc,srclens,srcamp)
!*    :lmax (int)          : Maximum multipole of output noise spectrum
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipole for reconstruction
!*    :wx0/1 [l] (double)  : 1/(CXX+NXX) (0:rlmax)
!*    :wxy0/1 [l] (double) : CXY/(CXX+NXX) (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Multiplying square of L(L+1)/2, i.e., convergence ('k') or lensing potential ('', default)
!*
!*  Returns:
!*    :Nl [l] (double)   : Noise cross spectrum with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: gtype
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: wx0, wx1, wxy0, wxy1, a0a1, b0b1, a0b1, a1b0
  double precision, intent(out), dimension(0:lmax) :: Nl
  !internal
  integer :: rL(2), l
  double precision, dimension(lmax) :: lk2
  double precision, dimension(:,:), allocatable :: W1, W2
  double precision, dimension(2,lmax) :: SG
  !opt4py :: gtype = ''

  rL = (/rlmin,rlmax/)

  !gradient-leg and C-inverse leg
  allocate(W1(2,rL(1):rL(2)),W2(2,rL(1):rL(2))); W1=0d0; W2=0d0

  do l = rlmin, rlmax
    W1(1,l) = wx0(l)  * wx1(l)  * a0a1(l)
    W2(1,l) = wxy0(l) * wxy1(l) * b0b1(l)
    W1(2,l) = wx0(l)  * wxy1(l) * a0b1(l)
    W2(2,l) = wxy0(l) * wx1(l)  * a1b0(l)
  end do

  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lensamp','lenssrc')
    call get_lfac(lmax,gtype,lk2)
    call Kernels_lenstau(rL,W1(1,:),W2(1,:),SG(1,:),'S0')
    call Kernels_lenstau(rL,W1(2,:),W2(2,:),SG(2,:),'G0')
  case('amplens','srclens')
    call get_lfac(lmax,gtype,lk2)
    call Kernels_lenstau(rL,W1(1,:),W2(1,:),SG(1,:),'S0')
    call Kernels_lenstau(rL,W2(2,:),W1(2,:),SG(2,:),'G0')
  case('ampsrc')
    call get_lfac(lmax,gtype,lk2)
    call Kernels_tau(rL,W1(1,:),W2(1,:),SG(1,:),'S0')
    call Kernels_tau(rL,W2(2,:),W1(2,:),SG(2,:),'G0')
  case('srcamp')
    call get_lfac(lmax,gtype,lk2)
    call Kernels_tau(rL,W1(1,:),W2(1,:),SG(1,:),'S0')
    call Kernels_tau(rL,W2(2,:),W1(2,:),SG(2,:),'G0')
  end select

  Nl = 0d0
  do l = 1, lmax
    if (lk2(l)/=0d0)  Nl(l) = sum(SG(:,l))/dsqrt(lk2(l))
  end do

end subroutine xtt_asym


end module noise_spec



