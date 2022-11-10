!////////////////////////////////////////////////!
! Normalization of quadratic lens reconstruction
!////////////////////////////////////////////////!

module norm_lens
  use norm_quad
  implicit none

contains


subroutine qtt(lmax,rlmin,rlmax,TT,fTT,OCT,Ag,Ac,gtype,temp_arg)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :TT [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :fTT [l] (double)   : Theory TT spectrum used for the weights, with bounds (0:rlmax)

!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double)   : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double)   : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:temp_arg) :: TT,fTT, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''
  !internal
  double precision, dimension(2,0:lmax) :: Al

  call quad_tt('lens',lmax,rlmin,rlmax,TT,fTT,OCT,Al,gtype)
  Ag = Al(1,:)
  Ac = Al(2,:)

end subroutine qtt


subroutine qte(lmax,rlmin,rlmax,TE,fTE,OCT,OCE,Ag,Ac,gtype,temp_arg)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the TE quadratic estimator
!*
!*  Args:
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
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*

  implicit none
  !I/O
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:temp_arg) :: TE,fTE, OCT, OCE
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''
  !internal
  double precision, dimension(2,0:lmax) :: Al

  call quad_te('lens',lmax,rlmin,rlmax,TE,fTE,OCT,OCE,Al,gtype)
  Ag = Al(1,:)
  Ac = Al(2,:)

end subroutine qte

subroutine qtb(lmax,rlmin,rlmax,TE,fTE,OCT,OCB,Ag,Ac,gtype,temp_arg)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the TB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :TE [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :fTE [l] (double)  : Theory TE spectrum for weogjts, with bounds (0:rlmax)

!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:temp_arg) :: TE,fTE, OCT, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''
  !internal
  double precision, dimension(2,0:lmax) :: Al

  call quad_tb('lens',lmax,rlmin,rlmax,TE,fTE,OCT,OCB,Al,gtype)
  Ag = Al(1,:)
  Ac = Al(2,:)

end subroutine qtb


subroutine qee(lmax,rlmin,rlmax,EE,fEE,OCE,Ag,Ac,gtype,temp_arg)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the E-mode quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in) , dimension(0:temp_arg) :: EE,fEE, OCE
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''
  !internal
  double precision, dimension(2,0:lmax) :: Al

  call quad_ee('lens',lmax,rlmin,rlmax,EE,fEE,OCE,Al,gtype)
  Ag = Al(1,:)
  Ac = Al(2,:)

end subroutine qee


subroutine qeb(lmax,rlmin,rlmax,EE,fEE,OCE,OCB,Ag,Ac,gtype,temp_arg)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
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
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in) , dimension(0:temp_arg) :: EE,fEE, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''
  !internal
  double precision, allocatable :: BB(:)
  double precision, dimension(2,0:lmax) :: Al

  allocate(BB(0:rlmax))
  BB = 0d0
  call quad_eb('lens',lmax,rlmin,rlmax,EE,fEE,OCE,OCB,BB,BB,Al,gtype)
  Ag = Al(1,:)
  Ac = Al(2,:)
  deallocate(BB)

end subroutine qeb


subroutine qbb(lmax,rlmin,rlmax,BB,fBB,OCB,Ag,Ac,gtype,temp_arg)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the B-mode quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :BB [l] (double)  : Theory BB spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: BB,fBB, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''
  !internal
  double precision, dimension(2,0:lmax) :: Al

  call quad_bb('lens',lmax,rlmin,rlmax,BB,fBB,OCB,Al,gtype)
  Ag = Al(1,:)
  Ac = Al(2,:)

end subroutine qbb


subroutine qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype,temp_arg)
!*  Correlation between unnormalized TT and TE quadratic estimators
!*
!*  Args:
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
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: fCTT, fCTE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''

  call quad_ttte('lens',lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype)

end subroutine qttte


subroutine qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic,gtype,temp_arg)
!*  Correlation between unnormalized TT and EE quadratic estimators
!*
!*  Args:
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
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: fCTT, fCEE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''

  call quad_ttee('lens',lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic,gtype)

end subroutine qttee


subroutine qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype,temp_arg)
!*  Correlation between unnormalized TE and EE quadratic estimators
!*
!*  Args:
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
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: fCEE, fCTE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''

  call quad_teee('lens',lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype)

end subroutine qteee


subroutine qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic,gtype,temp_arg)
!*  Correlation between unnormalized TB and EB quadratic estimators
!*
!*  Args:
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
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: fCEE, fCBB, fCTE, OCT, OCE, OCTE, OCB
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''

  call quad_tbeb('lens',lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic,gtype)

end subroutine qtbeb


subroutine qall(QDO,lmax,rlmin,rlmax,fC,fwC,OC,Ag,Ac,Nlg,Nlc,gtype,temp_arg)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
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
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  logical, intent(in), dimension(6) :: QDO
  integer, intent(in) :: rlmin, rlmax, lmax
  double precision, intent(in), dimension(4,0:temp_arg) :: fC, fwC,OC
  double precision, intent(out), dimension(6,0:lmax) :: Ag, Ac, Nlg, Nlc
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''

  call quad_all('lens',QDO,lmax,rlmin,rlmax,fC,fwC,OC,Ag,Ac,Nlg,Nlc,gtype)

end subroutine qall


subroutine qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,fCE,OCE,OCB,Cpp,Ag,Ac,iter,conv,temp_arg_E,temp_arg_B,temp_arg_p)
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
  integer :: temp_arg_E, temp_arg_B, temp_arg_p ! these are removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg_E) :: CE,fCE, OCE
  double precision, intent(in), dimension(0:temp_arg_B) :: OCB
  double precision, intent(in), dimension(0:temp_arg_p) :: Cpp
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  integer, intent(in) :: iter
  double precision, intent(in) :: conv
  !opt4py :: iter = 1
  !opt4py :: conv = 1e-6

  call quad_eb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,fCE,OCE,OCB,Cpp,Ag,Ac,iter,conv)

end subroutine qeb_iter


subroutine ttt(lmax,rlmin,rlmax,fC,OCT,Ag,gtype,temp_arg)
!*  Unnormalized response between lensing potential and amplitude modulation from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double)   : Cross normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''

  call quad_xtt('lensamp',lmax,rlmin,rlmax,fC,OCT,Ag,gtype)

end subroutine ttt


subroutine stt(lmax,rlmin,rlmax,fC,OCT,Ag,gtype,temp_arg)
!*  Unnormalized response between lensing potential and poisson sources/inhomogeneous noise with the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double)   : Cross normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: temp_arg ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:temp_arg) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag
  !optional
  character(1), intent(in) :: gtype
  !opt4py :: gtype = ''

  call quad_xtt('lenssrc',lmax,rlmin,rlmax,fC,OCT,Ag,gtype)

end subroutine stt



end module norm_lens


