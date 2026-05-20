!////////////////////////////////////////////////////!
! * Normalization of (real) quadratic estimator
!////////////////////////////////////////////////////!

module norm_quad
  use alkernel, only: kernels_lens, kernels_tau, kernels_rot, kernels_lenstau, get_lfac, clbb_est
  implicit none

  private kernels_lens, kernels_tau, kernels_rot, kernels_lenstau, get_lfac, clbb_est

contains


subroutine quad_tt(est,lmax,rlmin,rlmax,TT,fTT,OCT,Al,lfac)
!*  Normalization of reconstructed fields from the temperature quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :TT[l] (double)    : Theory TT spectrum, with bounds (0:rlmax)
!*    :fTT[l] (double)   : True TT spectrum (fTT=TT for forecast), with bounds (0:rlmax)
!*    :OCT[l] (double)   : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)        : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: TT, fTT, OCT
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: rL(2), l
  double precision, dimension(lmax) :: lk2
  double precision, dimension(4,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qtt): observed cltt is zero'
  end do

  W(1,:) = 1d0 / OCT(rlmin:rlmax)
  W(2,:) = TT(rlmin:rlmax)* fTT(rlmin:rlmax)/ OCT(rlmin:rlmax)
  W(3,:) = TT(rlmin:rlmax) / OCT(rlmin:rlmax)
  W(4,:) = fTT(rlmin:rlmax) / OCT(rlmin:rlmax)

  lk2 = 1d0
  SG  = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    call Kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'S0')
    call Kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'G0')
  case('amp')
    call Kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'S0')
    call Kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'G0')
  case('src')
    call Kernels_tau(rL,W(1,:),W(1,:),SG(1,1,:),'S0')
    call Kernels_tau(rL,W(1,:),W(1,:),SG(2,1,:),'G0')
    SG = SG/4d0
  end select

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  Al(2,1) = 0d0

end subroutine quad_tt


subroutine quad_tt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,Al,lfac)
!*  Normalization of reconstructed fields from the temperature quadratic estimator (asymmetric case)
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :glmin/glmax (int) : Minimum/Maximum multipole of gradient leg
!*    :llmin/llmax (int) : Minimum/Maximum multipole of C-inverse leg
!*    :rlmax (int)       : Minimum/Maximum multipole of TT
!*    :TT [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCTG [l] (double) : Observed TT spectrum for gradient leg, with bounds (0:glmax)
!*    :OCTL [l] (double) : Observed TT spectrum for C-inverse leg, with bounds (0:llmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmax, glmin, glmax, llmin, llmax
  double precision, intent(in), dimension(0:rlmax) :: TT
  double precision, intent(in), dimension(0:glmax) :: OCTG
  double precision, intent(in), dimension(0:llmax) :: OCTL
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: rL(2), l
  double precision, dimension(lmax) :: lk2
  double precision, dimension(:,:), allocatable :: WL, WG
  double precision, dimension(2,2,lmax) :: SG

  if (max(llmax,glmax)/=rlmax)  stop 'error (qtt): max(glmax,llmax) should be lmax of TT'

  do l = glmin, glmax
    if (OCTG(l)==0d0) stop 'error (qtt): observed cltt is zero for gradient leg'
  end do

  do l = llmin, llmax
    if (OCTL(l)==0d0) stop 'error (qtt): observed cltt is zero for C-inverse leg'
  end do

  rL = (/min(glmin,llmin),rlmax/)

  !gradient-leg and C-inverse leg
  allocate(WL(2,rL(1):rL(2)),WG(2,rL(1):rL(2))); WL=0d0; WG=0d0

  select case(est)
  case('lens','amp')
    do l = glmin, glmax
      WG(1,l) = TT(l)**2 / OCTG(l)
      WG(2,l) = TT(l) / OCTG(l)
    end do
    do l = llmin, llmax
      WL(1,l) = 1d0 / OCTL(l)
      WL(2,l) = TT(l) / OCTL(l)
    end do
  case('src')
    do l = glmin, glmax
      WG(1,l) = 1d0 / OCTG(l)
    end do
    do l = llmin, llmax
      WL(1,l) = 1d0 / OCTL(l)
    end do
  end select

  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    call Kernels_lens(rL,WL(1,:),WG(1,:),SG(1,:,:),'S0')
    call Kernels_lens(rL,WL(2,:),WG(2,:),SG(2,:,:),'G0')
  case('amp')
    call Kernels_tau(rL,WL(1,:),WG(1,:),SG(1,1,:),'S0')
    call Kernels_tau(rL,WL(2,:),WG(2,:),SG(2,1,:),'G0')
  case('src')
    call Kernels_tau(rL,WL(1,:),WG(1,:),SG(1,1,:),'S0')
    call Kernels_tau(rL,WL(1,:),WG(1,:),SG(2,1,:),'G0')
    SG = SG/4d0
  end select

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  Al(2,1) = 0d0

end subroutine quad_tt_asym


subroutine quad_te(est,lmax,rlmin,rlmax,TE,fTE,OCT,OCE,Al,lfac)
!*  Normalization of reconstructed fields from the TE quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,rot,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :TE [l] (double)   : Theory TE spectrum, with bounds (0:rlmax)
!*    :fTE [l] (double)  : True TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: TE,fTE, OCT, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(8,rlmin:rlmax) :: W
  double precision, dimension(4,2,lmax) :: SG

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qte): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qte): observed clee is zero'
  end do

  W(1,:) = 1d0/OCT(rlmin:rlmax)
  W(2,:) = TE(rlmin:rlmax)*fTE(rlmin:rlmax)/OCE(rlmin:rlmax)
  
  W(3,:) = TE(rlmin:rlmax)/OCT(rlmin:rlmax)
  W(4,:) = fTE(rlmin:rlmax)/OCE(rlmin:rlmax)
  
  W(5,:) = 1d0/OCE(rlmin:rlmax)
  W(6,:) = TE(rlmin:rlmax)*fTE(rlmin:rlmax)/OCT(rlmin:rlmax)

  ! this term is necessary if TE /= fTE
  W(7,:) = fTE(rlmin:rlmax)/OCT(rlmin:rlmax)
  W(8,:) = TE(rlmin:rlmax)/OCE(rlmin:rlmax)
  
  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'S0')
    call kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'Gc')
    call kernels_lens(rL,W(5,:),W(6,:),SG(3,:,:),'Sp')
    call kernels_lens(rL,W(7,:),W(8,:),SG(4,:,:),'Gc')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'S0')
    call kernels_tau(rL,W(3,:),W(4,:),SG(2,1,:),'Gc')
    call kernels_tau(rL,W(5,:),W(6,:),SG(3,1,:),'Sp')
    call kernels_tau(rL,W(7,:),W(8,:),SG(4,1,:),'Gc')
  case('rot')
    call kernels_rot(rL,W(5,:),W(6,:),SG(3,1,:),'Sp')
  case('src')
    call kernels_tau(rL,W(1,:),W(5,:),SG(1,1,:),'S0')
    call kernels_tau(rL,W(1,:),W(5,:),SG(2,1,:),'Gc')
    call kernels_tau(rL,W(1,:),W(5,:),SG(3,1,:),'Sp')
    SG = SG/4d0
    SG(2,:,:) = 2d0*SG(2,:,:)
  end select
  
  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(2,1) = 0d0
  case('rot')
    Al(1,1) = 0d0
  end select

end subroutine quad_te


subroutine quad_tb(est,lmax,rlmin,rlmax,TE,fTE,OCT,OCB,Al,lfac)
!*  Normalization of reconstructed fields from the TB quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,rot,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :TE [l] (double)   : Theory TE spectrum, with bounds (0:rlmax)
!*    :fTE [l] (double)  : True TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double)  : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: TE,fTE, OCT, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(2,rlmin:rlmax) :: W
  double precision, dimension(2,lmax) :: SG

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qtb): observed cltt is zero'
    if (OCB(l)==0d0) stop 'error (qtb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCB(rlmin:rlmax)
  W(2,:) = TE(rlmin:rlmax)*fTE(rlmin:rlmax) / OCT(rlmin:rlmax)

  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lens(rL,W(1,:),W(2,:),SG,'Sm')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,:),'Sm')
  case('rot')
    call kernels_rot(rL,W(1,:),W(2,:),SG(1,:),'Sm')
  case('src')
    call kernels_tau(rL,W(1,:),W(1,:),SG(1,:),'Sm')
    SG = SG/4d0
  end select

  Al = 0d0
  do l = 1, lmax
    if (SG(1,l)/=0d0)  Al(1,l) = lk2(l)/SG(1,l)
    if (SG(2,l)/=0d0)  Al(2,l) = lk2(l)/SG(2,l)
  end do
  select case(est)
  case('lens','amp','src')
    Al(1,1) = 0d0
  end select

end subroutine quad_tb


subroutine quad_ee(est,lmax,rlmin,rlmax,EE,fEE,OCE,Al,lfac)
!*  Normalization of reconstructed amplitude modulation from the EE quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,rot,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)   : Theory EE spectrum, with bounds (0:rlmax)
!*    :fEE [l] (double)  : True EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE,fEE, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(4,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (qee): observed clee is zero'
  end do

  W(1,:) = 1d0/OCE(rlmin:rlmax)
  W(2,:) = EE(rlmin:rlmax)*fEE(rlmin:rlmax) / OCE(rlmin:rlmax)
  W(3,:) = EE(rlmin:rlmax) / OCE(rlmin:rlmax)
  W(4,:) = fEE(rlmin:rlmax) / OCE(rlmin:rlmax)

  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sp')
    call kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'Gp')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('rot')
    call kernels_rot(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_rot(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('src')
    call kernels_tau(rL,W(1,:),W(1,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(1,:),W(1,:),SG(2,1,:),'Gp')
    SG = SG/4d0
  end select


  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(2,1) = 0d0
  case('rot')
    Al(1,1) = 0d0
  end select

end subroutine quad_ee


subroutine quad_eb(est,lmax,rlmin,rlmax,EE,fEE,OCE,OCB,BB,fBB,Al,lfac)
!*  Normalization of reconstructed fields from the EB quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,rot,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)   : Theory EE spectrum, with bounds (0:rlmax)
!*    :fEE [l] (double)  : True EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double)  : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optionals): 
!*    :BB [l] (double)   : Theory BB spectrum, with bounds (0:rlmax)
!*    :fBB [l] (double)  : True BB spectrum, with bounds (0:rlmax)
!*    :lfac (str)        : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE,fEE, BB,fBB, OCE, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(8,rlmin:rlmax) :: W
  double precision, dimension(4,2,lmax) :: SG

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (qeb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCE(rlmin:rlmax)
  W(2,:) = BB(rlmin:rlmax)*fBB(rlmin:rlmax)/ OCB(rlmin:rlmax)
  W(3,:) = EE(rlmin:rlmax)/OCE(rlmin:rlmax)
  W(4,:) = fBB(rlmin:rlmax)/OCB(rlmin:rlmax)
  W(5,:) = 1d0/OCB(rlmin:rlmax)
  W(6,:) = EE(rlmin:rlmax)*fEE(rlmin:rlmax) / OCE(rlmin:rlmax)

  ! if fEE/=EE or fBB=BB, this term is necessary
  W(7,:) = fEE(rlmin:rlmax)/OCE(rlmin:rlmax)
  W(8,:) = BB(rlmin:rlmax)/OCB(rlmin:rlmax)
  
  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    if (sum(BB)/=0d0) then
      call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sm')
      call kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'Gm')
      call kernels_lens(rL,W(7,:),W(8,:),SG(4,:,:),'Gm')
    end if
    call kernels_lens(rL,W(5,:),W(6,:),SG(3,:,:),'Sm')
  case('amp')
    if (sum(BB)/=0d0) then
      call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sm')
      call kernels_tau(rL,W(3,:),W(4,:),SG(2,1,:),'Gm')
      call kernels_tau(rL,W(7,:),W(8,:),SG(4,1,:),'Gm')
    end if
    call kernels_tau(rL,W(5,:),W(6,:),SG(3,1,:),'Sm')
  case('rot')
    if (sum(BB)/=0d0) then
      call kernels_rot(rL,W(1,:),W(2,:),SG(1,1,:),'Sm')
      call kernels_rot(rL,W(3,:),W(4,:),SG(2,1,:),'Gm')
      call kernels_rot(rL,W(7,:),W(8,:),SG(4,1,:),'Gm')
    end if
    call kernels_rot(rL,W(5,:),W(6,:),SG(3,1,:),'Sm')
  case('src')
    call kernels_tau(rL,W(1,:),W(5,:),SG(1,1,:),'Sm')
    call kernels_tau(rL,W(1,:),W(5,:),SG(2,1,:),'Gm')
    call kernels_tau(rL,W(1,:),W(5,:),SG(3,1,:),'Sm')
    SG = SG/4d0
    SG(2,:,:) = 2d0*SG(2,:,:)
  end select

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(1,1) = 0d0
  end select

end subroutine quad_eb


subroutine quad_bb(est,lmax,rlmin,rlmax,BB,fBB,OCB,Al,lfac)
!*  Normalization of reconstructed fields from the BB quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,rot,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :BB [l] (double)   : Theory BB spectrum, with bounds (0:rlmax)
!*    :fBB [l] (double)  : True BB spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double)  : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: BB,fBB, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(4,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCB(l)==0d0) stop 'error (qbb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCB(rlmin:rlmax)
  W(2,:) = BB(rlmin:rlmax)*fBB(rlmin:rlmax) / OCB(rlmin:rlmax)
  W(3,:) = BB(rlmin:rlmax) / OCB(rlmin:rlmax)
  W(4,:) = fBB(rlmin:rlmax) / OCB(rlmin:rlmax)

  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sp')
    call kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'Gp')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('rot')
    call kernels_rot(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_rot(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('src')
    call kernels_tau(rL,W(1,:),W(1,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(1,:),W(1,:),SG(2,1,:),'Gp')
    SG = SG/4d0
  end select

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(2,1) = 0d0
  case('rot')
    Al(1,1) = 0d0
  end select

end subroutine quad_bb


subroutine rfunc_ttte(est,lmax,rlmin,rlmax,tCTT,fCTE,A,B,R_p,R_m,lfac)
!*  Correlation between TT and TE weights, R_L^{TT,TE}[A,B]
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :tCTT [l] (double) : Theory TT spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double) : True TE spectrum, with bounds (0:rlmax)
!*    :A/B [l] (double)  : Any power spectra, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :R_p [l] (double) : R_L for even-type estimators (e.g. phi), with bounds (0:lmax)
!*    :R_m [l] (double) : R_L for odd-type estimators (e.g. curl), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: tCTT, fCTE, A, B
  double precision, intent(out), dimension(0:lmax) :: R_p, R_m
  !internal
  integer :: i, l, rL(2)
  character(2) :: SG_type(4)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(4,rlmin:rlmax) :: W0, W1
  double precision, dimension(4,2,lmax) :: SG

  rL  = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    W0(1,l) = A(l)
    W1(1,l) = B(l)*tCTT(l)*fCTE(l)
    W0(2,l) = A(l)*fCTE(l)
    W1(2,l) = B(l)*tCTT(l)
    W0(3,l) = B(l)*fCTE(l)
    W1(3,l) = A(l)*tCTT(l)
    W0(4,l) = B(l)
    W1(4,l) = A(l)*tCTT(l)*fCTE(l)
  end do

  lk2 = 1d0
  SG_type = (/'S0','Gc','G0','Sc'/)
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    do i = 1, 4
      call kernels_lens(rL,W0(i,:),W1(i,:),SG(i,:,:),SG_type(i))
    end do
  case('amp')
    do i = 1, 4
      call kernels_tau(rL,W0(i,:),W1(i,:),SG(i,1,:),SG_type(i))
    end do
  end select

  R_p = 0d0
  R_m = 0d0
  do l = 1, lmax
    R_p(l) = sum(SG(:,1,l))/lk2(l)
    R_m(l) = sum(SG(:,2,l))/lk2(l)
  end do

end subroutine rfunc_ttte


subroutine rfunc_ttee(est,lmax,rlmin,rlmax,tCTT,fCEE,A,B,R_p,R_m,lfac)
!*  Correlation between TT and EE weights, R_L^{TT,EE}[A,B]
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :tCTT [l] (double) : Theory TT/EE spectrum, with bounds (0:rlmax)
!*    :fCEE [l] (double) : True TT/EE spectrum, with bounds (0:rlmax)
!*    :A/B [l] (double)  : Any power spectra, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :R_p [l] (double) : R_L for even-type estimators (e.g. phi), with bounds (0:lmax)
!*    :R_m [l] (double) : R_L for odd-type estimators (e.g. curl), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: tCTT, fCEE, A, B
  double precision, intent(out), dimension(0:lmax) :: R_p, R_m
  !internal
  integer :: i, l, rL(2)
  character(2) :: SG_type(4)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(4,rlmin:rlmax) :: W0, W1
  double precision, dimension(4,2,lmax) :: SG

  do l = rlmin, rlmax
    W0(1,l) = A(l)
    W1(1,l) = B(l)*tCTT(l)*fCEE(l)
    W0(2,l) = A(l)*fCEE(l)
    W1(2,l) = B(l)*tCTT(l)
    W0(3,l) = B(l)*fCEE(l)
    W1(3,l) = A(l)*tCTT(l)
    W0(4,l) = B(l)
    W1(4,l) = A(l)*tCTT(l)*fCEE(l)
  end do

  lk2 = 1d0
  rL  = (/rlmin,rlmax/)
  SG_type = (/'Sc','Gc','Gc','Sc'/)
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    do i = 1, 4
      call kernels_lens(rL,W0(i,:),W1(i,:),SG(i,:,:),SG_type(i))
    end do
  case('amp')
    do i = 1, 4
      call kernels_tau(rL,W0(i,:),W1(i,:),SG(i,1,:),SG_type(i))
    end do
  end select

  R_p = 0d0
  R_m = 0d0
  do l = 1, lmax
    R_p(l) = sum(SG(:,1,l))/lk2(l)
    R_m(l) = sum(SG(:,2,l))/lk2(l)
  end do

end subroutine rfunc_ttee


subroutine rfunc_teet(est,lmax,rlmin,rlmax,tCTE,fCTE,A,B,R_p,R_m,lfac)
!*  Correlation between TE and ET weights, R_L^{TE,ET}[A,B]
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :tCTE [l] (double) : Theory TE spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double) : True TE spectrum, with bounds (0:rlmax)
!*    :A/B [l] (double)  : Any power spectra, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :R_p [l] (double) : R_L for even-type estimators (e.g. phi), with bounds (0:lmax)
!*    :R_m [l] (double) : R_L for odd-type estimators (e.g. curl), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: tCTE, fCTE, A, B
  double precision, intent(out), dimension(0:lmax) :: R_p, R_m
  !internal
  integer :: i, l, rL(2)
  character(2) :: SG_type(4)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(4,rlmin:rlmax) :: W0, W1
  double precision, dimension(4,2,lmax) :: SG

  do l = rlmin, rlmax
    W0(1,l) = A(l)
    W1(1,l) = B(l)*tCTE(l)*fCTE(l)
    W0(2,l) = A(l)*fCTE(l)
    W1(2,l) = B(l)*tCTE(l)
    W0(3,l) = B(l)*fCTE(l)
    W1(3,l) = A(l)*tCTE(l)
    W0(4,l) = B(l)
    W1(4,l) = A(l)*tCTE(l)*fCTE(l)
  end do

  lk2 = 1d0
  rL  = (/rlmin,rlmax/)
  SG_type = (/'Sc','G0','Gp','Sc'/)
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    do i = 1, 4
      call kernels_lens(rL,W0(i,:),W1(i,:),SG(i,:,:),SG_type(i))
    end do
  case('amp')
    do i = 1, 4
      call kernels_tau(rL,W0(i,:),W1(i,:),SG(i,1,:),SG_type(i))
    end do
  end select

  R_p = 0d0
  R_m = 0d0
  do l = 1, lmax
    R_p(l) = sum(SG(:,1,l))/lk2(l)
    R_m(l) = sum(SG(:,2,l))/lk2(l)
  end do

end subroutine rfunc_teet


subroutine rfunc_teee(est,lmax,rlmin,rlmax,tCTE,fCEE,A,B,R_p,R_m,lfac)
!*  Correlation between TE and EE weights, R_L^{TE,EE}[A,B]
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :tCTE [l] (double) : Theory TE spectrum, with bounds (0:rlmax)
!*    :fCEE [l] (double) : True EE spectrum, with bounds (0:rlmax)
!*    :A/B [l] (double)  : Any power spectra, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :R_p [l] (double) : R_L for even-type estimators (e.g. phi), with bounds (0:lmax)
!*    :R_m [l] (double) : R_L for odd-type estimators (e.g. curl), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: tCTE, fCEE, A, B
  double precision, intent(out), dimension(0:lmax) :: R_p, R_m
  !internal
  integer :: i, l, rL(2)
  character(2) :: SG_type(4)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(4,rlmin:rlmax) :: W0, W1
  double precision, dimension(4,2,lmax) :: SG

  do l = rlmin, rlmax
    W0(1,l) = A(l)
    W1(1,l) = B(l)*tCTE(l)*fCEE(l)
    W0(2,l) = A(l)*fCEE(l)
    W1(2,l) = B(l)*tCTE(l)
    W0(3,l) = B(l)*fCEE(l)
    W1(3,l) = A(l)*tCTE(l)
    W0(4,l) = B(l)
    W1(4,l) = A(l)*tCTE(l)*fCEE(l)
  end do

  lk2 = 1d0
  rL  = (/rlmin,rlmax/)
  SG_type = (/'Sc','Gc','Gp','Sp'/)
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    do i = 1, 4
      call kernels_lens(rL,W0(i,:),W1(i,:),SG(i,:,:),SG_type(i))
    end do
  case('amp')
    do i = 1, 4
      call kernels_tau(rL,W0(i,:),W1(i,:),SG(i,1,:),SG_type(i))
    end do
  end select

  R_p = 0d0
  R_m = 0d0
  do l = 1, lmax
    R_p(l) = sum(SG(:,1,l))/lk2(l)
    R_m(l) = sum(SG(:,2,l))/lk2(l)
  end do

end subroutine rfunc_teee


subroutine rfunc_tbeb(est,lmax,rlmin,rlmax,tCTE,fCEE,fCBB,A,B,R_p,R_m,lfac)
!*  Correlation between TB and EB weights, R_L^{TE,EE}[A,B]
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :tCTE [l] (double) : Theory TE spectrum, with bounds (0:rlmax)
!*    :fCEE [l] (double) : True EE spectrum, with bounds (0:rlmax)
!*    :fCEB [l] (double) : True BB spectrum, with bounds (0:rlmax)
!*    :A/B [l] (double)  : Any power spectra, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :R_p [l] (double) : R_L for even-type estimators (e.g. phi), with bounds (0:lmax)
!*    :R_m [l] (double) : R_L for odd-type estimators (e.g. curl), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: tCTE, fCEE, fCBB, A, B
  double precision, intent(out), dimension(0:lmax) :: R_p, R_m
  !internal
  integer :: i, l, rL(2)
  character(2) :: SG_type(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(2,rlmin:rlmax) :: W0, W1
  double precision, dimension(2,2,lmax) :: SG

  do l = rlmin, rlmax
    W0(1,l) = B(l)*fCBB(l)
    W1(1,l) = A(l)*tCTE(l)
    W0(2,l) = B(l)
    W1(2,l) = A(l)*tCTE(l)*fCEE(l)
  end do

  lk2 = 1d0
  rL  = (/rlmin,rlmax/)
  SG_type = (/'Gm','Sm'/)
  select case(est)
  case('lens')
    call get_lfac(lmax,lfac,lk2)
    do i = 1, 2
      call kernels_lens(rL,W0(i,:),W1(i,:),SG(i,:,:),SG_type(i))
    end do
  case('amp')
    do i = 1, 2
      call kernels_tau(rL,W0(i,:),W1(i,:),SG(i,1,:),SG_type(i))
    end do
  end select

  R_p = 0d0
  R_m = 0d0
  do l = 1, lmax
    R_p(l) = sum(SG(:,1,l))/lk2(l)
    R_m(l) = sum(SG(:,2,l))/lk2(l)
  end do

end subroutine rfunc_tbeb


subroutine quad_ttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TT and TE quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fCTT [l] (double) : Theory TT spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double) : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double) : Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT, fCTE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l
  double precision, dimension(0:rlmax) :: A, B

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_qttte): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (norm_qttte): observed clee is zero'
  end do

  A = 0d0
  B = 0d0
  do l = rlmin, rlmax
    A(l) = 1d0/OCT(l)
    B(l) = OCTE(l)/(OCT(l)*OCE(l))
  end do

  call rfunc_ttte(est,lmax,rlmin,rlmax,fCTT,fCTE,A,B,Ig,Ic,lfac)

end subroutine quad_ttte


subroutine quad_ttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TT and EE quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fCTT [l] (double) : Theory TT spectrum, with bounds (0:rlmax)
!*    :fCEE [l] (double) : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double) : Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT, fCEE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l
  double precision, dimension(0:rlmax) :: A, B

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qttee): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qttee): observed clee is zero'
  end do

  do l = rlmin, rlmax
    A(l) = OCTE(l)/(OCT(l)*OCE(l))
    B(l) = OCTE(l)/(OCT(l)*OCE(l))
  end do

  call rfunc_ttee(est,lmax,rlmin,rlmax,fCTT,fCEE,A,B,Ig,Ic,lfac)

end subroutine quad_ttee


subroutine quad_teee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TE and EE quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,rot,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fCEE [l] (double) : Theory EE spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double) : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double) : Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCTE,OCT,OCE,OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l
  double precision, dimension(0:rlmax) :: A, B

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qteee): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qteee): observed clee is zero'
  end do

  do l = rlmin, rlmax
    A(l) = OCTE(l)/(OCT(l)*OCE(l))
    B(l) = 1d0/OCE(l)
  end do
  
  call rfunc_teee(est,lmax,rlmin,rlmax,fCTE,fCEE,A,B,Ig,Ic,lfac)

end subroutine quad_teee


subroutine quad_tbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TB and EB quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type (lens,amp,rot,src)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fCEE [l] (double) : Theory EE spectrum, with bounds (0:rlmax)
!*    :fCBB [l] (double) : Theory BB spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double)  : Observed BB spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double) : Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCBB,fCTE,OCT,OCE,OCTE,OCB
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l
  double precision, dimension(0:rlmax) :: A, B

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qtbeb): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qtbeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (qtbeb): observed clbb is zero'
  end do

  do l = rlmin, rlmax
    A(l) = OCTE(l)/(OCT(l)*OCE(l))
    B(l) = 1d0/OCB(l)
  end do

  call rfunc_tbeb(est,lmax,rlmin,rlmax,fCTE,fCEE,fCBB,A,B,Ig,Ic,lfac)

end subroutine quad_tbeb


subroutine quad_mv(lmax,QDO,Al,Il,MV,Nl)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :lmax (int):    Maximum multipole of the output power spectra
!*    :QDO[6] (bool): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. Currently, BB is always False
!*    :Al [5,l] (double): Normalizations of each estimator (TT, TE, EE, TB, EB). 
!*    :Il [4,l] (double): Correlation between different estimators (TTxTE, TTxEE, TExEE, TBxEB).
!*
!*  Returns:
!*    :MV [l] (double):   Normalization of the MV estimator, with bounds (0:lmax)
!*    :Nl [6,l] (double): Weights for each estimator (TT, TE, EE, TB, EB, BB=0), with bounds (0:lmax)
!*
  implicit none
  ![input]
  integer, intent(in) :: lmax
  logical, intent(in), dimension(6) :: QDO
  double precision, intent(in), dimension(5,0:lmax) :: Al
  double precision, intent(in), dimension(4,0:lmax) :: Il
  double precision, intent(out), dimension(0:lmax) :: MV
  double precision, intent(out), dimension(6,0:lmax) :: Nl
  !internal
  integer :: qn = 5, QTT = 1, QTE = 2, QTB = 4, QEE = 3, QEB = 5!, QBB = 6
  integer :: X, Y, qmax, i, id(6), l
  double precision :: d1, d2, M1(5,5)
  double precision, allocatable :: M(:,:)

  id = 0

  qmax = qn
 
  MV = 0d0
  Nl = 0d0

  do l = 2, lmax

    ! noise covariance
    allocate(M(qmax,qmax));  M = 0d0

    if (QDO(QTT).and.QDO(QTE)) M(1,2) = Il(1,l)*Al(QTT,l)*Al(QTE,l)
    if (QDO(QTT).and.QDO(QEE)) M(1,3) = Il(2,l)*Al(QTT,l)*Al(QEE,l)
    if (QDO(QTE).and.QDO(QEE)) M(2,3) = Il(3,l)*Al(QTE,l)*Al(QEE,l)
    if (QDO(QTB).and.QDO(QEB)) M(4,5) = Il(4,l)*Al(QTB,l)*Al(QEB,l)
    do X = 1, qn
      M(X,X) = 1d10 ! some large value for not used estimator
      if (QDO(X)) M(X,X) = Al(X,l)
      do Y = X + 1, qn
        if(QDO(X).and.QDO(Y)) M(Y,X) = M(X,Y)
      end do
    end do

    ! inverting M with explict expression (symmetric)
    d1 = M(1,1)*(M(2,2)*M(3,3)-M(2,3)**2) + M(1,2)*(M(2,3)*M(3,1)-M(2,1)*M(3,3)) + M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
    d2 = M(4,4)*M(5,5) - M(4,5)**2
    M1 = 0d0
    ! components
    M1(1,1) = M(2,2)*M(3,3) - M(2,3)**2
    M1(1,2) = M(1,3)*M(3,2) - M(1,2)*M(3,3)
    M1(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
    M1(2,2) = M(1,1)*M(3,3) - M(1,3)**2
    M1(2,3) = M(1,3)*M(2,1) - M(1,1)*M(2,3)
    M1(3,3) = M(1,1)*M(2,2) - M(1,2)**2
    M1(4,4) = M(5,5)
    M1(4,5) = -M(4,5)
    M1(5,5) = M(4,4)
    ! symmetric
    M1(2,1) = M1(1,2)
    M1(3,1) = M1(1,3)
    M1(3,2) = M1(2,3)
    M1(5,4) = M1(4,5)
    ! final
    M(1:3,1:3) = M1(1:3,1:3)/d1
    M(4:5,4:5) = M1(4:5,4:5)/d2

    MV(l)   = 1d0/sum(M)
    Nl(:,l) = sum(M,dim=2)

    deallocate(M)

  end do

end subroutine quad_mv


subroutine quad_gmv(est,lmax,rlmin,rlmax,tC,fC,OC,Ag,Ac,lfac,th_vary)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :est (str)             : Estimator type (lens,amp,rot,src)
!*    :lmax (int)            : Maximum multipole of the output power spectra
!*    :rlmin/rlmax (int)     : Minimum/Maximum multipole of CMB for reconstruction
!*    :tC/fC/OC [l] (double) : Theory/True/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax) 
!*
!*  Args(optional):
!*    :lfac (str)          : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*    :th_vary (bool)      : Vary theory or not, default = False
!*
!*  Returns:
!*    :Ag [6,l] (double)  : Normalization of the TT, TE, EE, TB, EB, and GMV estimators for lensing potential, with bounds (6,0:lmax)
!*    :Ac [6,l] (double)  : Same as Ag but for curl mode
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  logical, intent(in) :: th_vary
  integer, intent(in) :: rlmin, rlmax, lmax
  double precision, intent(in), dimension(4,0:rlmax) :: tC, fC, OC
  double precision, intent(out), dimension(6,0:lmax) :: Ag, Ac
  !internal
  integer :: l, TT = 1, EE = 2, BB = 3, TE = 4
  double precision :: rho_sq
  double precision, dimension(0:rlmax) :: A, B
  double precision, dimension(4,0:rlmax) :: tOC
  double precision, dimension(6,2,0:lmax) :: Al
  double precision, dimension(5,2,0:lmax) :: Rl, Sl

  write(*,*) 'norm qGMV'

  ! Modified CMB power spectra
  tOC = 0d0
  do l = rlmin, rlmax
    rho_sq    = OC(TE,l)**2 / (OC(TT,l)*OC(EE,l))
    tOC(TT,l) = OC(TT,l) * ( 1d0 - rho_sq )
    tOC(TE,l) = OC(TE,l) * ( 1d0 - 1d0/rho_sq )
    tOC(EE,l) = OC(EE,l) * ( 1d0 - rho_sq )
    tOC(BB,l) = OC(BB,l)
  end do

  ! Each estimator normalization
  Al = 0d0
  call quad_tt(est,lmax,rlmin,rlmax,tC(TT,:),fC(TT,:),tOC(TT,:),Al(1,:,:),lfac)
  call quad_te(est,lmax,rlmin,rlmax,tC(TE,:),fC(TE,:),tOC(TT,:),tOC(EE,:),Al(2,:,:),lfac)
  call quad_ee(est,lmax,rlmin,rlmax,tC(EE,:),fC(EE,:),tOC(EE,:),Al(3,:,:),lfac)
  call quad_tb(est,lmax,rlmin,rlmax,tC(TE,:),fC(TE,:),tOC(TT,:),tOC(BB,:),Al(4,:,:),lfac)
  call quad_eb(est,lmax,rlmin,rlmax,tC(EE,:),fC(EE,:),tOC(EE,:),tOC(BB,:),tC(BB,:),fC(BB,:),Al(5,:,:),lfac)
  Ag = Al(:,1,:)
  Ac = Al(:,2,:)

  ! Correlations
  A = 0d0
  B = 0d0
  do l = rlmin, rlmax
    A(l) = 1d0/tOC(TT,l)
    B(l) = 1d0/tOC(TE,l)
  end do
  call rfunc_ttte(est,lmax,rlmin,rlmax,tC(TT,:),fC(TE,:),A,B,Rl(1,1,:),Rl(1,2,:),lfac)

  if (th_vary) then
    call rfunc_ttte(est,lmax,rlmin,rlmax,fC(TT,:),tC(TE,:),A,B,Sl(1,1,:),Sl(1,2,:),lfac)
  end if
  
  do l = rlmin, rlmax
    A(l) = 1d0/tOC(TE,l)
    B(l) = 1d0/tOC(TE,l)
  end do
  call rfunc_ttee(est,lmax,rlmin,rlmax,tC(TT,:),fC(EE,:),A,B,Rl(2,1,:),Rl(2,2,:),lfac)

  if (th_vary) then
    call rfunc_ttee(est,lmax,rlmin,rlmax,fC(TT,:),tC(EE,:),A,B,Sl(2,1,:),Sl(2,2,:),lfac)
  end if

  do l = rlmin, rlmax
    A(l) = 1d0/tOC(TE,l)
    B(l) = 1d0/tOC(TE,l)
  end do
  call rfunc_teet(est,lmax,rlmin,rlmax,tC(TE,:),fC(TE,:),A,B,Rl(3,1,:),Rl(3,2,:),lfac)

  do l = rlmin, rlmax
    A(l) = 1d0/tOC(TE,l)
    B(l) = 1d0/tOC(EE,l)
  end do
  call rfunc_teee(est,lmax,rlmin,rlmax,tC(TE,:),fC(EE,:),A,B,Rl(4,1,:),Rl(4,2,:),lfac)

  if (th_vary) then
    call rfunc_teee(est,lmax,rlmin,rlmax,fC(TE,:),tC(EE,:),A,B,Sl(4,1,:),Sl(4,2,:),lfac)
  end if

  do l = rlmin, rlmax
    A(l) = 1d0/tOC(TE,l)
    B(l) = 1d0/tOC(BB,l)
  end do
  call rfunc_tbeb(est,lmax,rlmin,rlmax,tC(TE,:),fC(EE,:),fC(BB,:),A,B,Rl(5,1,:),Rl(5,2,:),lfac)

  if (th_vary) then
    call rfunc_tbeb(est,lmax,rlmin,rlmax,fC(TE,:),tC(EE,:),tC(BB,:),A,B,Sl(5,1,:),Sl(5,2,:),lfac)
  end if

  if (th_vary) then
    call quad_gmv_sum(lmax,1d0,Al(1:5,1,:),Rl(:,1,:),Sl(:,1,:),Ag(6,0:lmax))
    call quad_gmv_sum(lmax,-1d0,Al(1:5,2,:),Rl(:,2,:),Sl(:,2,:),Ac(6,0:lmax))
  else
    call quad_gmv_sum(lmax,1d0,Al(1:5,1,:),Rl(:,1,:),Rl(:,1,:),Ag(6,0:lmax))
    call quad_gmv_sum(lmax,-1d0,Al(1:5,2,:),Rl(:,2,:),Rl(:,2,:),Ac(6,0:lmax))
  end if
  

end subroutine quad_gmv


subroutine quad_gmv_sum(lmax,p,Al,Rl,Sl,MV)
!*  Compute GMV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :lmax (int):        Maximum multipole of the output power spectra
!*    :Al [5,l] (double): Normalizations of each estimator (TT, TE, EE, TB, EB). 
!*    :Rl [5,l] (double): Correlations, R_L^{TT,TE}, R_L^{TT,EE}, R_L^{TE,ET}, R_L^{TE,EE}, and R_L^{TB,EB}.
!*    :Sl [5,l] (double): Counterpart of R_L if tC/=fC or R_L if tC=fC.
!*    :p (double)       : parity of the estimator, e.g. phi for 1, curl for -1
!*
!*  Returns:
!*    :MV [l] (double):   Normalization of the MV estimator, with bounds (0:lmax)
!*
  implicit none
  ![input]
  integer, intent(in) :: lmax
  double precision, intent(in) :: p
  double precision, intent(in), dimension(5,0:lmax) :: Al
  double precision, intent(in), dimension(5,0:lmax) :: Rl, Sl
  double precision, intent(out), dimension(0:lmax) :: MV
  !internal
  integer :: l

  MV = 0d0
  do l = 2, lmax
    MV(l) = sum(1d0/Al(:,l)) + Rl(1,l) + Sl(1,l) + (0.5d0)*(Rl(2,l)+Sl(2,l)) + p*Rl(3,l) + Rl(4,l) + Sl(4,l) + Rl(5,l) + Sl(5,l)
    MV(l) = 1d0/MV(l)
  end do

end subroutine quad_gmv_sum


subroutine quad_all(est,QDO,lmax,rlmin,rlmax,fC,fwC,OC,Ag,Ac,Nlg,Nlc,lfac)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :est (str)          : Estimator type (lens,amp,rot,src)
!*    :QDO[6] (bool)      : Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. 
!*    :lmax (int)         : Maximum multipole of the output power spectra
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC/OC [l] (double) : Theory/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax) 
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ag [6,l] (double)  : Normalization of the TT, TE, EE, TB, EB, and MV estimators for lensing potential, with bounds (6,0:lmax)
!*    :Ac [6,l] (double)  : Same as Ag but for curl mode
!*    :Nlg [6,l] (double) : Weights for TT, TE, EE, TB, EB, and BB (=0) estimators for lensing potential, with bounds (6,0:lmax)
!*    :Nlc [6,l] (double) : Same as Nlg but for curl mode
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  logical, intent(in), dimension(6) :: QDO
  integer, intent(in) :: rlmin, rlmax, lmax
  double precision, intent(in), dimension(4,0:rlmax) :: fC, fwC,OC
  double precision, intent(out), dimension(6,0:lmax) :: Ag, Ac, Nlg, Nlc
  !internal
  character(1) :: gt
  integer :: TT = 1, EE = 2, BB = 3, TE = 4
  double precision, dimension(:,:), allocatable :: Ilg, Ilc
  double precision, dimension(6,2,0:lmax) :: Al

  !//// interface ////!
  Ag  = 0d0
  Ac  = 0d0
  Nlg = 0d0
  Nlc = 0d0
  if (QDO(1))  call quad_tt(est,lmax,rlmin,rlmax,fC(TT,:),fwC(TT,:),OC(TT,:),Al(1,:,:),lfac)
  if (QDO(2))  call quad_te(est,lmax,rlmin,rlmax,fC(TE,:),fwC(TE,:),OC(TT,:),OC(EE,:),Al(2,:,:),lfac)
  if (QDO(3))  call quad_ee(est,lmax,rlmin,rlmax,fC(EE,:),fwC(EE,:),OC(EE,:),Al(3,:,:),lfac)
  if (QDO(4))  call quad_tb(est,lmax,rlmin,rlmax,fC(TE,:),fwC(TE,:),OC(TT,:),OC(BB,:),Al(4,:,:),lfac)
  if (QDO(5))  call quad_eb(est,lmax,rlmin,rlmax,fC(EE,:),fwC(EE,:),OC(EE,:),OC(BB,:),fC(BB,:),fwC(BB,:),Al(5,:,:),lfac)
  Ag = Al(:,1,:)
  Ac = Al(:,2,:)

  allocate(Ilg(4,0:lmax),Ilc(4,0:lmax))
  if (QDO(1).and.QDO(2))  call quad_ttte(est,lmax,rlmin,rlmax,fC(TT,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(1,:),Ilc(1,:),lfac)
  if (QDO(1).and.QDO(3))  call quad_ttee(est,lmax,rlmin,rlmax,fC(TT,:),fC(EE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(2,:),Ilc(2,:),lfac)
  if (QDO(2).and.QDO(3))  call quad_teee(est,lmax,rlmin,rlmax,fC(EE,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(3,:),Ilc(3,:),lfac)
  if (QDO(4).and.QDO(5))  call quad_tbeb(est,lmax,rlmin,rlmax,fC(EE,:),fC(BB,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(BB,:),OC(TE,:), &
       Ilg(4,:),Ilc(4,:),lfac)
  call quad_mv(lmax,QDO,Ag(1:5,0:lmax),Ilg,Ag(6,0:lmax),Nlg)
  call quad_mv(lmax,QDO,Ac(1:5,0:lmax),Ilc,Ac(6,0:lmax),Nlc)
  deallocate(Ilg,Ilc)

end subroutine quad_all


subroutine quad_eb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,fCE,OCE,OCB,Cpp,Ag,Ac,iter,conv)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization
!*    :elmax (int)      : Maximum multipole of input EE spectra, CE and OCE
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :dlmin/dlmax (int): Minimum/Maximum multipole of E mode and lensing potential for delensing
!*    :CE [l] (double)  : Theory EE angular power spectrum, with bounds (0:elmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:elmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*    :Cpp [l] (double) : Theory lensing potential spectrum, with bounds (0:dlmax)
!*
!*  Args(optional):
!*    :iter (int)    : number of iteration, default to 1 (no iteration)
!*    :conv (double) : a parameter for convergence the iteration, default to 0.001
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, elmax, rlmin, rlmax, dlmin, dlmax
  double precision, intent(in), dimension(0:elmax) :: CE,fCE, OCE
  double precision, intent(in), dimension(0:rlmax) :: OCB
  double precision, intent(in), dimension(0:dlmax) :: Cpp
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  integer, intent(in) :: iter
  double precision, intent(in) :: conv
  !internal
  integer :: i, n, l
  double precision :: ratio
  double precision :: Al(2,0:dlmax), rCBB(0:rlmax), BB(0:rlmax)
  !opt4py :: iter = 1
  !opt4py :: conv = 0.001

  if (elmax<dlmax.or.elmax<rlmax) stop 'error (qeb_iter): does not support elmax<dlmax or elmax<rlmax'

  !initial values
  ratio = 1d0
  rCBB  = OCB
  Al   = 0d0
  BB   = 0d0

  do n = 1, iter !loop for iteration 

    !lensing reconstruction with EB
    call quad_eb('lens',dlmax,rlmin,rlmax,CE(0:rlmax),fCE(0:rlmax),OCE(0:rlmax),rCBB,BB(0:rlmax),BB(0:rlmax),Al(:,0:dlmax),'')

    !convergence check using gradient mode
    if (n>=2) then
      ratio = (sum(Ag)/sum(Al(1,:))-1d0)/dble(dlmax)
      write(*,*) n, ratio
    end if
    Ag = Al(1,:)
    Ac = Al(2,:)

    if (abs(ratio) < conv) exit

    !delensing with EB-estimator
    call clbb_est((/rlmin,rlmax/),(/dlmin,dlmax/),CE(1:dlmax),Cpp(1:dlmax),OCE(1:dlmax)-CE(1:dlmax),Al(1,1:dlmax),rCBB(1:rlmax))
    rCBB = OCB - rCBB !delensed B-mode

    if(n==iter) write(*,*) 'not well converged'

  end do

end subroutine quad_eb_iter


! //////////////////// !
! Response Function !
! //////////////////// !

subroutine quad_xtt(est,lmax,rlmin,rlmax,fC,OCT,Rxy,lfac)
!*  Unnormalized response for the symmetric temperature quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator combination (lensamp,lenssrc,ampsrc)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)        : Convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Rxy [l] (double)   : Unnormalized response with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Rxy
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2, W3
  double precision, dimension(2,lmax) :: SG

  write(*,*) 'Response (TT)'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_xtt): observed cltt is zero'
  end do

  W1 = 1d0 / OCT(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2 / OCT(rlmin:rlmax)
  W3 = fC(rlmin:rlmax) / OCT(rlmin:rlmax)

  lk2 = 1d0
  SG  = 0d0
  select case(est)
  case('lensamp')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lenstau(rL,W1,W2,SG(1,:),'S0')
    call kernels_lenstau(rL,W3,W3,SG(2,:),'G0')
  case('lenssrc')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lenstau(rL,W1,W3,SG(1,:),'S0')
    call kernels_lenstau(rL,W1,W3,SG(2,:),'G0')
    SG = SG*0.5d0
  case('ampsrc')
    call Kernels_tau(rL,W1,W3,SG(1,:),'S0')
    call Kernels_tau(rL,W1,W3,SG(2,:),'G0')
    SG = SG*0.5d0
  end select

  Rxy = 0d0
  do l = 1, lmax
    Rxy(l) = sum(SG(:,l))/dsqrt(lk2(l))
  end do

end subroutine quad_xtt


subroutine quad_xtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,Rxy,lfac)
!*  Unnormalized response for the asymmetric temperature quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator combination (lensamp,lenssrc,amplens,ampsrc,srclens,srcamp)
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :glmin/glmax (int) : Minimum/Maximum multipole of gradient leg
!*    :llmin/llmax (int) : Minimum/Maximum multipole of C-inverse leg
!*    :rlmax (int)       : Minimum/Maximum multipole of TT
!*    :TT [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCTG [l] (double) : Observed TT spectrum for gradient leg, with bounds (0:glmax)
!*    :OCTL [l] (double) : Observed TT spectrum for C-inverse leg, with bounds (0:llmax)
!*
!*  Args(optional):
!*    :lfac (str)        : Convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Rl [l] (double)   : Unnormalized response with bounds (0:lmax)
!*
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmax, glmin, glmax, llmin, llmax
  double precision, intent(in), dimension(0:rlmax) :: TT
  double precision, intent(in), dimension(0:glmax) :: OCTG
  double precision, intent(in), dimension(0:llmax) :: OCTL
  double precision, intent(out), dimension(0:lmax) :: Rxy
  !internal
  integer :: rL(2), l
  double precision, dimension(lmax) :: lk2
  double precision, dimension(:,:), allocatable :: WL, WG
  double precision, dimension(2,lmax) :: SG

  write(*,*) 'Response (TT,asym)'

  if (max(llmax,glmax)/=rlmax)  stop 'error (qtt): max(glmax,llmax) should be lmax of TT'

  do l = glmin, glmax
    if (OCTG(l)==0d0) stop 'error (qtt): observed cltt is zero for gradient leg'
  end do

  do l = llmin, llmax
    if (OCTL(l)==0d0) stop 'error (qtt): observed cltt is zero for C-inverse leg'
  end do

  rL = (/min(glmin,llmin),rlmax/)

  !gradient-leg and C-inverse leg
  allocate(WL(2,rL(1):rL(2)),WG(2,rL(1):rL(2))); WL=0d0; WG=0d0

  select case(est)
  case('lensamp','amplens')
    do l = glmin, glmax
      WG(1,l) = TT(l)**2 / OCTG(l)
      WG(2,l) = TT(l) / OCTG(l)
    end do
    do l = llmin, llmax
      WL(1,l) = 1d0 / OCTL(l)
      WL(2,l) = TT(l) / OCTL(l)
    end do
  case('lenssrc','ampsrc')
    do l = glmin, glmax
      WG(1,l) = TT(l) / OCTG(l)
      WG(2,l) = TT(l) / OCTG(l)
    end do
    do l = llmin, llmax
      WL(1,l) = 1d0 / OCTL(l)
      WL(2,l) = 1d0 / OCTL(l)
    end do
  case('srclens','srcamp')
    do l = glmin, glmax
      WG(1,l) = TT(l) / OCTG(l)
      WG(2,l) = 1d0 / OCTG(l)
    end do
    do l = llmin, llmax
      WL(1,l) = 1d0 / OCTL(l)
      WL(2,l) = TT(l) / OCTL(l)
    end do
  end select

  lk2 = 1d0
  SG = 0d0
  select case(est)
  case('lensamp')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lenstau(rL,WL(1,:),WG(1,:),SG(1,:),'S0')
    call kernels_lenstau(rL,WL(2,:),WG(2,:),SG(2,:),'G0')
  case('amplens')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lenstau(rL,WL(1,:),WG(1,:),SG(1,:),'S0')
    call kernels_lenstau(rL,WG(2,:),WL(2,:),SG(2,:),'G0')
  case('lenssrc')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lenstau(rL,WL(1,:),WG(1,:),SG(1,:),'S0')
    call kernels_lenstau(rL,WL(2,:),WG(2,:),SG(2,:),'G0')
    SG = SG*0.5d0
  case('srclens')
    call get_lfac(lmax,lfac,lk2)
    call kernels_lenstau(rL,WL(1,:),WG(1,:),SG(1,:),'S0')
    call kernels_lenstau(rL,WG(2,:),WL(2,:),SG(2,:),'G0')
    SG = SG*0.5d0
  case('ampsrc')
    call Kernels_tau(rL,WL(1,:),WG(1,:),SG(1,:),'S0')
    call Kernels_tau(rL,WL(2,:),WG(2,:),SG(2,:),'G0')
    SG = SG*0.5d0
  case('srcamp')
    call Kernels_tau(rL,WL(1,:),WG(1,:),SG(1,:),'S0')
    call Kernels_tau(rL,WG(2,:),WL(2,:),SG(2,:),'G0')
    SG = SG*0.5d0
  end select

  Rxy = 0d0
  do l = 1, lmax
    Rxy(l) = sum(SG(:,l))/dsqrt(lk2(l))
  end do

end subroutine quad_xtt_asym


end module norm_quad



