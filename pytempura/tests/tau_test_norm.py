# Written by Darby Kramer, fixed Sept 2022. This is modeled mostly from tempura/pytempura/tests/test_norm.py 
# except this is an executable file and not a callable function.

# This file should be moved to and run from the directory containing tempura and orphics. 
# You will have to make the imports work if you want to run it from within the test folder. 

import numpy as np
from matplotlib.pyplot import *
import orphics.cosmology as cosmology

# import pytempura
import tempura.pytempura as cs

#import norm file
from tempura.pytempura.norm import get_norms

def test(noiselevel=20):
    
    noiselevel = noiselevel # TT noise in µK' , set to 0 if no noise is desired. Can be whitenoise level or specific noise curve.

    Tcmb  = 2.726e6    # CMB temperature
    Lmax  = 3000       # maximum multipole of output normalization
    lmax  = 3000
    rlmin, rlmax = 100, 3000  # CMB multipole range for reconstruction
    L = np.linspace(0,Lmax,Lmax+1)
    l = L.copy()
    Lfac = (L*(L+1.))**2/(2*np.pi) # this is not used because of the next line.
    Lfac_tau = 1/(2*np.pi) # Implemented this factor to get a flat tau tt normalization plot.

    thloc = "falafel/data/cosmo2017_10K_acc3"
    theory = cosmology.loadTheorySpectraFromCAMB(thloc,get_dimensionless=False)
    
    if type(noiselevel) is np.ndarray:
        assert len(noiselevel) == len(theory.uCl('TT',L))
        noisecl = noiselevel
    
    else:
        print("No array detected: constructing Noise Cl")
        noisecl = []
        
        for i in range(lmax+1):
            noisecl.append(((noiselevel * np.pi / (180 * 60))**2)) # conversion from arcminutes to radians^2

    cl = np.zeros((4,lmax+1)) # TT, EE, BB, TE
    cl[:,2:] = [theory.uCl('TT',L[2:]),theory.uCl('EE',L[2:]),\
                    theory.uCl('BB',L[2:]),theory.uCl('TE',L[2:])] 

    noise_tt = noisecl
    noise_ee, noise_bb = 2*noise_tt, 2*noise_tt

    ocl1 = [np.add(cl[0],noise_tt)*1., (cl[1])*1., (cl[2])*1., (cl[3])*1.]  # observed Cl, noise added to tt lcl
    ocl2 = [(cl[0])*1., (cl[1])*1., (cl[2])*1., (cl[3])*1.]  # observed Cl, NO noise added to tt lcl

    # Construction of dictionaries necessary for the norm functions.

    ocls1 = {'TT' : ocl1[0], 'EE': ocl1[1], 'BB': ocl1[2], 'TE': ocl1[3]} 
    ocls2 = {'TT' : ocl2[0], 'EE': ocl2[1], 'BB': ocl2[2], 'TE': ocl2[3]} 

    ucls = {'TT' : cl[0], 'EE': cl[1], 'BB': cl[2], 'TE': cl[3]} 


    # Norm Calculator. Cannot call tau and lens in the same function call. Returns Aest as a dictionary of arrays

    Aest1 = get_norms(['tt'],ucls,ocls1,rlmin,rlmax,coupling=["tau"]) # only tt and eb are implemented for tau.
    Aest2 = get_norms(['tt'],ucls,ocls2,rlmin,rlmax,coupling=["tau"])

    # plot spectra

    clf()
    figure(figsize=(10,8))
    title(r"Normalization vs. $\ell$")
    xlim(0,3000)
    yscale('log')

    for est in Aest1:
        plot(L,Lfac_tau*Aest1[est], label='Noise, '+est, color = 'k', linewidth=5)

    for est in Aest2: 
        plot(L,Lfac_tau*Aest2[est], label='No Noise, '+est, color='r') 
        
#     plot(L, noisecl, label=r"$N_{\ell}^{TT}$", color='b')

    ylabel(r"$A_{\tau}$ or $N_{\ell}^{TT}$")
    xlabel(r"L or $\ell$")
    legend()
    savefig('tautest.png') #save figure


