from . import _libtempura

def qtt(est,lmax,rlmin,rlmax,TT,OCT,gtype=''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the temperature quadratic estimator

  Args:
    :est (*str*): Estimator type (lens,amp,src)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TT [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = _libtempura.norm_general.qtt(est,lmax,rlmin,rlmax,TT,OCT,gtype):
  """
  return _libtempura.norm_general.qtt(est,lmax,rlmin,rlmax,TT,OCT,gtype)

def qtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,gtype=''):
  """
  Normalization of reconstructed fields from the temperature quadratic estimator (asymmetric case)

  Args:
    :est (*str*): Estimator type (lens,amp,src)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :glmin/glmax (*int*): Minimum/Maximum multipole of gradient leg
    :llmin/llmax (*int*): Minimum/Maximum multipole of C-inverse leg
    :rlmax (*int*): Minimum/Maximum multipole of TT
    :TT[*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCTG[*l*] (*double*): Observed TT spectrum for gradient leg, with bounds (0:glmax)
    :OCTL[*l*] (*double*): Observed TT spectrum for C-inverse leg, with bounds (0:llmax)

  Args(optional):
    :gtype (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = _libtempura.norm_general.qtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,gtype):
  """
  return _libtempura.norm_general.qtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,gtype)

def qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,gtype=''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the TE quadratic estimator

  Args:
    :est (*str*): Estimator type (lens,amp,rot,src)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = _libtempura.norm_general.qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,gtype):
  """
  return _libtempura.norm_general.qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,gtype)

def qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,gtype=''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the TB quadratic estimator

  Args:
    :est (*str*): Estimator type (lens,amp,rot,src)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = _libtempura.norm_general.qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,gtype):
  """
  return _libtempura.norm_general.qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,gtype)

def qee(est,lmax,rlmin,rlmax,EE,OCE,gtype=''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the E-mode quadratic estimator

  Args:
    :est (*str*): Estimator type (lens,amp,rot,src)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = _libtempura.norm_general.qee(est,lmax,rlmin,rlmax,EE,OCE,gtype):
  """
  return _libtempura.norm_general.qee(est,lmax,rlmin,rlmax,EE,OCE,gtype)

def qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,gtype=''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :est (*str*): Estimator type (lens,amp,rot,src)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = _libtempura.norm_general.qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,gtype):
  """
  return _libtempura.norm_general.qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,gtype)

def qbb(est,lmax,rlmin,rlmax,BB,OCB,gtype=''):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the B-mode quadratic estimator

  Args:
    :est (*str*): Estimator type (lens,amp,rot,src)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :BB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = _libtempura.norm_general.qbb(est,lmax,rlmin,rlmax,BB,OCB,gtype):
  """
  return _libtempura.norm_general.qbb(est,lmax,rlmin,rlmax,BB,OCB,gtype)

def qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype=''):
  """
  Correlation between unnormalized TT and TE quadratic estimators

  Args:
    :est (*str*): Estimator type (lens,amp,src)
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
    :Ig,Ic = _libtempura.norm_general.qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype):
  """
  return _libtempura.norm_general.qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype)

def qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype=''):
  """
  Correlation between unnormalized TT and EE quadratic estimators

  Args:
    :est (*str*): Estimator type (lens,amp,src)
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
    :Ig,Ic = _libtempura.norm_general.qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype):
  """
  return _libtempura.norm_general.qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype)

def qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype=''):
  """
  Correlation between unnormalized TE and EE quadratic estimators

  Args:
    :est (*str*): Estimator type (lens,amp,rot,src)
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
    :Ig,Ic = _libtempura.norm_general.qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype):
  """
  return _libtempura.norm_general.qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype)

def qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype=''):
  """
  Correlation between unnormalized TB and EB quadratic estimators

  Args:
    :est (*str*): Estimator type (lens,amp,rot,src)
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
    :Ig,Ic = _libtempura.norm_general.qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype):
  """
  return _libtempura.norm_general.qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype)

def qall(est,QDO,lmax,rlmin,rlmax,fC,OC,gtype=''):
  """
  Compute MV estimator normalization. Currently BB is ignored. 

  Args:
    :est (*str*): Estimator type
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
    :Ag,Ac,Nlg,Nlc = _libtempura.norm_general.qall(est,QDO,lmax,rlmin,rlmax,fC,OC,gtype):
  """
  return _libtempura.norm_general.qall(est,QDO,lmax,rlmin,rlmax,fC,OC,gtype)

def qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter=1,conv=1e-6):
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
    :Ag,Ac = _libtempura.norm_general.qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter,conv):
  """
  return _libtempura.norm_general.qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter,conv)

def xtt(est,lmax,rlmin,rlmax,fC,OCT,gtype=''):
  """
  Unnormalized response for the symmetric temperature quadratic estimator

  Args:
    :est (*str*): Estimator combination (lensamp,lenssrc,ampsrc)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Rxy [*l*] (*double*): Unnormalized response, with bounds (0:lmax)

  Usage:
    :Rxy = _libtempura.norm_general.xtt(est,lmax,rlmin,rlmax,fC,OCT,gtype):
  """
  return _libtempura.norm_general.xtt(est,lmax,rlmin,rlmax,fC,OCT,gtype)

def xtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,gtype=''):
  """
  Unnormalized response for the asymmetric temperature quadratic estimator

  Args:
    :est (*str*): Estimator combination (lensamp,lenssrc,amplens,ampsrc,srclens,srcamp)
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :glmin/glmax (*int*): Minimum/Maximum multipole of gradient leg
    :llmin/llmax (*int*): Minimum/Maximum multipole of C-inverse leg
    :rlmax (*int*): Minimum/Maximum multipole of TT
    :TT[*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCTG[*l*] (*double*): Observed TT spectrum for gradient leg, with bounds (0:glmax)
    :OCTL[*l*] (*double*): Observed TT spectrum for C-inverse leg, with bounds (0:llmax)

  Args(optional):
    :gtype (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Rxy [*l*] (*double*): Unnormalized response, with bounds (0:lmax)

  Usage:
    :Rxy = _libtempura.norm_general.xtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,gtype):
  """
  return _libtempura.norm_general.xtt_asym(est,lmax,glmin,glmax,llmin,llmax,rlmax,TT,OCTG,OCTL,gtype)

