!////////////////////////////////////////////////!
! Normalization of (real) quadratic estimator reconstruction
!////////////////////////////////////////////////!

module norm_general
  use norm_quad
  implicit none

contains


subroutine qtt(est,lmax,rlmin,rlmax,TT,OCT,Al,gtype,N)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the temperature quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :TT [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:N) :: TT, OCT
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_qtt(est,lmax,rlmin,rlmax,TT,OCT,Al,gtype)

end subroutine qtt


subroutine qtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCT0,OCT1,Al,gtype,NTT,N0,N1)
!*  Normalization of reconstructed fields from the temperature quadratic estimator (asymmetric case)
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :glmin/glmax (int) : Minimum/Maximum multipole of gradient leg
!*    :llmin/llmax (int) : Minimum/Maximum multipole of C-inverse leg
!*    :rlmax (int)       : Minimum/Maximum multipole of TT
!*    :TT [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT0 [l] (double) : Observed TT spectrum for gradient leg, with bounds (0:glmax)
!*    :OCT1 [l] (double) : Observed TT spectrum for C-inverse leg, with bounds (0:llmax)
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
  integer :: NTT, N0, N1 ! this argument is removed by f2py since it appears in the size of an input array argument
  integer, intent(in) :: lmax, rlmax, glmin, glmax, llmin, llmax
  double precision, intent(in), dimension(0:NTT) :: TT
  double precision, intent(in), dimension(0:N0) :: OCT0
  double precision, intent(in), dimension(0:N1) :: OCT1
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_tt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCT0,OCT1,Al,gtype)

end subroutine qtt_asym


subroutine qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,Al,gtype,N)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the TE quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :TE [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*

  implicit none
  !I/O
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:N) :: TE, OCT, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_te(est,lmax,rlmin,rlmax,TE,OCT,OCE,Al,gtype)

end subroutine qte


subroutine qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,Al,gtype,N)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the TB quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :TE [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:N) :: TE, OCT, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_tb(est,lmax,rlmin,rlmax,TE,OCT,OCB,Al,gtype)

end subroutine qtb


subroutine qee(est,lmax,rlmin,rlmax,EE,OCE,Al,gtype,N)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the E-mode quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in) , dimension(0:N) :: EE, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_ee(est,lmax,rlmin,rlmax,EE,OCE,Al,gtype)

end subroutine qee


subroutine qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,Al,gtype,N)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in) , dimension(0:N) :: EE, OCE, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  double precision, allocatable :: BB(:)

  allocate(BB(0:rlmax))
  BB = 0d0
  call quad_eb(est,lmax,rlmin,rlmax,EE,OCE,OCB,BB,Al,gtype)
  deallocate(BB)

end subroutine qeb


subroutine qbb(est,lmax,rlmin,rlmax,BB,OCB,Al,gtype,N)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the B-mode quadratic estimator
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :BB [l] (double)  : Theory BB spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:N) :: BB, OCB
  double precision, intent(out), dimension(0:lmax) :: Al
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_bb(est,lmax,rlmin,rlmax,BB,OCB,Al,gtype)

end subroutine qbb


subroutine qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype,N)
!*  Correlation between unnormalized TT and TE quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCTT [l] (double): Theory TT spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double): Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:N) :: fCTT, fCTE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_ttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype)

end subroutine qttte


subroutine qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic,gtype,N)
!*  Correlation between unnormalized TT and EE quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCTT [l] (double): Theory TT spectrum, with bounds (0:rlmax)
!*    :fCEE [l] (double): Theory EE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:N) :: fCTT, fCEE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_ttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic,gtype)

end subroutine qttee


subroutine qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype,N)
!*  Correlation between unnormalized TE and EE quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCEE [l] (double): Theory EE spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double): Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:N) :: fCEE,fCTE,OCT,OCE,OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_teee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype)

end subroutine qteee


subroutine qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic,gtype,N)
!*  Correlation between unnormalized TB and EB quadratic estimators
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCEE [l] (double): Theory EE spectrum, with bounds (0:rlmax)
!*    :fCBB [l] (double): Theory BB spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:N) :: fCEE,fCBB,fCTE,OCT,OCE,OCTE,OCB
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_tbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic,gtype)

end subroutine qtbeb


subroutine qall(est,QDO,lmax,rlmin,rlmax,fC,OC,Ag,Ac,Nlg,Nlc,gtype,N)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :est (str)         : Estimator type
!*    :QDO[6] (bool): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. 
!*    :lmax (int):    Maximum multipole of the output power spectra
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC/OC [l] (double): Theory/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax) 
!*
!*  Args(optional):
!*    :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
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
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  logical, intent(in), dimension(6) :: QDO
  integer, intent(in) :: rlmin, rlmax, lmax
  double precision, intent(in), dimension(4,0:N) :: fC, OC
  double precision, intent(out), dimension(6,0:lmax) :: Ag, Ac, Nlg, Nlc
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_all(est,QDO,lmax,rlmin,rlmax,fC,OC,Ag,Ac,Nlg,Nlc,gtype)

end subroutine qall


subroutine qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,Ag,Ac,iter,conv,Ne,Nr,Nd)
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
  integer :: Ne, Nr, Nd ! these are removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:Ne) :: CE, OCE
  double precision, intent(in), dimension(0:Nr) :: OCB
  double precision, intent(in), dimension(0:Nd) :: Cpp
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  integer, intent(in), optional :: iter
  double precision, intent(in), optional :: conv
  !f2py integer :: iter = 1
  !f2py double precision :: conv = 1e-6

  call quad_eb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,Ag,Ac,iter,conv)

end subroutine qeb_iter


subroutine xtt(xest,lmax,rlmin,rlmax,fC,OCT,Ag,gtype,N)
!*  Unnormalized response between two estimators from the temperature quadratic estimator
!*
!*  Args:
!*    :xest (str)        : Combination of two estimator types
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double)   : Unnormalized response, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: xest
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:N) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  call quad_x_tt(xest,lmax,rlmin,rlmax,fC,OCT,Ag,gtype)

end subroutine xtt


end module norm_general


