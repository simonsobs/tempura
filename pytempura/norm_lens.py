from . import _libtempura 

def qtt(lmax,rlmin,rlmax,fC,fCw,OCT,gtype= ''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qtt(lmax,rlmin,rlmax,fC,OCT,gtype):
  """
  return _libtempura.norm_lens.qtt(lmax,rlmin,rlmax,fC,fCw,OCT,gtype)

def qtt_generic(lmax,rlmin,rlmax,fC,fCw,OCT,gtype= ''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qtt(lmax,rlmin,rlmax,fC,OCT,gtype):
  """
  return _libtempura.norm_lens.qtt_generic(lmax,rlmin,rlmax,fC,OCT,gtype)

def qte(lmax,rlmin,rlmax,fC,fCw,OCT,OCE,gtype= ''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the TE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qte(lmax,rlmin,rlmax,fC,OCT,OCE,gtype):
  """
  return _libtempura.norm_lens.qte(lmax,rlmin,rlmax,fC,fCw,OCT,OCE,gtype)

def qtb(lmax,rlmin,rlmax,fC,fCw,OCT,OCB,gtype= ''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the TB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qtb(lmax,rlmin,rlmax,fC,OCT,OCB,gtype):
  """
  return _libtempura.norm_lens.qtb(lmax,rlmin,rlmax,fC,fCw,OCT,OCB,gtype)

def qee(lmax,rlmin,rlmax,fC,fCw,OCE,gtype= ''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the E-mode quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qee(lmax,rlmin,rlmax,fC,OCE,gtype):
  """
  return _libtempura.norm_lens.qee(lmax,rlmin,rlmax,fC,fCw,OCE,gtype)

def qeb(lmax,rlmin,rlmax,fC,fCw,OCE,OCB,gtype= ''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qeb(lmax,rlmin,rlmax,fC,OCE,OCB,gtype):
  """
  return _libtempura.norm_lens.qeb(lmax,rlmin,rlmax,fC,fCw,OCE,OCB,gtype)

def qbb(lmax,rlmin,rlmax,fC,fCw,OCB,gtype= ''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the B-mode quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qbb(lmax,rlmin,rlmax,fC,OCB,gtype):
  """
  return _libtempura.norm_lens.qbb(lmax,rlmin,rlmax,fC,fCw,OCB,gtype)

def qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype= ''):
  """
  Correlation between unnormalized TT and TE quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCTT [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :fCTE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = pytempura.norm_lens.qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype):
  """
  return _libtempura.norm_lens.qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype)

def qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype= ''):
  """
  Correlation between unnormalized TT and EE quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCTT [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :fCEE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = pytempura.norm_lens.qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype):
  """
  return _libtempura.norm_lens.qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype)

def qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype= ''):
  """
  Correlation between unnormalized TE and EE quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCEE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :fCTE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = pytempura.norm_lens.qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype):
  """
  return _libtempura.norm_lens.qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype)

def qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype= ''):
  """
  Correlation between unnormalized TB and EB quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCEE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :fCBB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = pytempura.norm_lens.qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype):
  """
  return _libtempura.norm_lens.qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype)

def qmv(lmax,QDO,Al,Il):
  """
  Compute MV estimator normalization. Currently BB is ignored. 

  Args:
    :lmax (*int*): Maximum multipole of the output power spectra
    :QDO[*6*] (*bool*): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. Currently, BB is always False
    :Al [*5,l*] (*double*): Normalizations of each estimator (TT, TE, EE, TB, EB).
    :Il [*4,l*] (*double*): Correlation between different estimators (TTxTE, TTxEE, TExEE, TBxEB).

  Returns:
    :MV [*l*] (*double*): Normalization of the MV estimator, with bounds (0:lmax)
    :Nl [*6,l*] (*double*): Weights for each estimator (TT, TE, EE, TB, EB, BB=0), with bounds (0:lmax)

  Usage:
    :MV,Nl = pytempura.norm_lens.qmv(lmax,QDO,Al,Il):
  """
  return _libtempura.norm_lens.qmv(lmax,QDO,Al,Il)

def qall(QDO,lmax,rlmin,rlmax,fC,fC1,OC,gtype= ''):
  """
  Compute MV estimator normalization. Currently BB is ignored. 

  Args:
    :QDO[*6*] (*bool*): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB.
    :lmax (*int*): Maximum multipole of the output power spectra
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC/OC [*l*] (*double*): Theory/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*6,l*] (*double*): Normalization of the TT, TE, EE, TB, EB, and MV estimators for lensing potential, with bounds (6,0:lmax)
    :Ac [*6,l*] (*double*): Same as Ag but for curl mode
    :Nlg [*6,l*] (*double*): Weights for TT, TE, EE, TB, EB, and BB (=0) estimators for lensing potential, with bounds (6,0:lmax)
    :Nlc [*6,l*] (*double*): Same as Nlg but for curl mode

  Usage:
    :Ag,Ac,Nlg,Nlc = pytempura.norm_lens.qall(QDO,lmax,rlmin,rlmax,fC,OC,gtype):
  """
  return _libtempura.norm_lens.qall(QDO,lmax,rlmin,rlmax,fC,fC1,OC,gtype)

def qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter= 1,conv= 1e-6):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization
    :elmax (*int*): Maximum multipole of input EE spectra, CE and OCE
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :dlmin/dlmax (*int*): Minimum/Maximum multipole of E mode and lensing potential for delensing
    :CE [*l*] (*double*): Theory EE angular power spectrum, with bounds (0:elmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:elmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)
    :Cpp [*l*] (*double*): Theory lensing potential spectrum, with bounds (0:dlmax)

  Args(optional):
    :iter (*int*): number of iteration, default to 1 (no iteration)
    :conv (*double*): a parameter for convergence the iteration, default to 0.001

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = pytempura.norm_lens.qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter,conv):
  """
  return _libtempura.norm_lens.qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter,conv)

def ttt(lmax,rlmin,rlmax,fC,OCT,gtype= ''):
  """
  Unnormalized response between lensing potential and amplitude modulation from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): Cross normalization, with bounds (0:lmax)

  Usage:
    :Ag = pytempura.norm_lens.ttt(lmax,rlmin,rlmax,fC,OCT,gtype):
  """
  return _libtempura.norm_lens.ttt(lmax,rlmin,rlmax,fC,OCT,gtype)

def stt(lmax,rlmin,rlmax,fC,OCT,gtype= ''):
  """
  Unnormalized response between lensing potential and poisson sources/inhomogeneous noise with the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): Cross normalization, with bounds (0:lmax)

  Usage:
    :Ag = pytempura.norm_lens.stt(lmax,rlmin,rlmax,fC,OCT,gtype):
  """
  return _libtempura.norm_lens.stt(lmax,rlmin,rlmax,fC,OCT,gtype)

