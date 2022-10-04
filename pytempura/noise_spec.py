from . import _libtempura

def qtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,gtype=''):
  """
  Noise spectrum of reconstructed fields from the temperature quadratic estimator (asymmetric case)

  Args:
    :lmax (*int*): Maximum multipole of output noise spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole for reconstruction
    :wx0/1 [*l*] (*double*): 1/(CXX+NXX) (0:rlmax)
    :wxy0/1 [*l*] (*double*): CXY/(CXX+NXX) (0:rlmax)

  Args(optional):
    :gtype (*str*): Multiplying square of L(L+1)/2, i.e., convergence ('k') or lensing potential ('', default)

  Returns:
    :Nl [*2,l*] (*double*): Noise spectrum (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Nl = _libtempura.noise_spec.qtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,gtype):
  """
  return _libtempura.noise_spec.qtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,gtype)

def xtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,gtype=''):
  """
  Noise spectrum of reconstructed fields from the temperature quadratic estimator (asymmetric case)

  Args:
    :est (*str*): Estimator combination (lensamp,lenssrc,amplens,ampsrc,srclens,srcamp)
    :lmax (*int*): Maximum multipole of output noise spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole for reconstruction
    :wx0/1 [*l*] (*double*): 1/(CXX+NXX) (0:rlmax)
    :wxy0/1 [*l*] (*double*): CXY/(CXX+NXX) (0:rlmax)

  Args(optional):
    :gtype (*str*): Multiplying square of L(L+1)/2, i.e., convergence ('k') or lensing potential ('', default)

  Returns:
    :Nl [*l*] (*double*): Noise cross spectrum with bounds (0:lmax)

  Usage:
    :Nl = _libtempura.noise_spec.xtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,gtype):
  """
  return _libtempura.noise_spec.xtt_asym(est,lmax,rlmin,rlmax,wx0,wxy0,wx1,wxy1,a0a1,b0b1,a0b1,a1b0,gtype)

