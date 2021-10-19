!////////////////////////////////////////////////////!
! * Normalization of quadratic src reconstruction
!////////////////////////////////////////////////////!

module norm_src
  use norm_quad
  implicit none

contains


subroutine qtt(lmax,rlmin,rlmax,OCT,As,N)
!*  Normalization of reconstructed src field from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :As [l] (double) : src field normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer :: N ! this argument is removed by f2py since it appears in the size of an input array argument
  double precision, intent(in), dimension(0:N) :: OCT
  double precision, intent(out), dimension(0:lmax) :: As
  !internal
  double precision, dimension(2,0:lmax) :: Al
  double precision, dimension(0:rlmax) :: TT

  TT = 0.5d0
  call quad_qtt('src',lmax,rlmin,rlmax,TT,OCT,Al,'')
  As = Al(1,:)

end subroutine qtt


end module norm_src

