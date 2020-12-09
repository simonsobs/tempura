import libtempura

def qtt(lmax,rlmin,rlmax,OCT):
  """
  Normalization of reconstructed src field from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Returns:
    :As [*l*] (*double*): src field normalization, with bounds (0:lmax)

  Usage:
    :As = tempura.norm_src.qtt(lmax,rlmin,rlmax,OCT):
  """
  return libtempura.norm_src.qtt(lmax,rlmin,rlmax,OCT)

