# Written by Darby Kramer, March 2021. This is modeled mostly from tempura/pytempura/tests/test_norm.py 
# except this is an executable file and not a callable function.

# This file should be moved to and run from the directory containing the whole tempura folder. 
# You will have to change the imports if you want to run it from within the test folder. 
# Sorry if there is a more preferred way to do this, I am still getting the hang of these imports.

import numpy as np
from matplotlib.pyplot import *
import pdb

# import pytempura
import tempura.pytempura as cs

#import norm file
from tempura.pytempura.norm import get_norms

Tcmb  = 2.726e6    # CMB temperature
Lmax  = 3000       # maximum multipole of output normalization
lmax  = 3000
rlmin, rlmax = 100, 3000  # CMB multipole range for reconstruction
L = np.linspace(0,Lmax,Lmax+1)
l = L.copy()
Lfac = (L*(L+1.))**2/(2*np.pi) # this is not used because of the next line.
Lfac_tau = 1/(2*np.pi) # Did this to get a flat tau tt normalization plot. Advised by Alex Van Engelen

lcl = np.zeros((4,lmax+1)) # TT, EE, BB, TE
lcl[:,2:] = np.loadtxt('tempura/pytempura/tests/data/cosmo2017_10K_acc3_lensedCls.dat',unpack=True,usecols=(1,2,3,4))[:,:lmax-1] 


lcl *= 2.*np.pi / (l**2+l+1e-30) / Tcmb**2

noise_tt = (1 * np.pi/180/60)**2 # conversion from arcminutes to radians

noise_ee, noise_bb = 2*noise_tt, 2*noise_tt


ocl = [(lcl[0]+noise_tt)*1., (lcl[1])*1., (lcl[2])*1., (lcl[3])*1.]  # observed Cl, noise added to tt lcl

# Construction of dictionaries necessary for the norm functions.

ocls = {'TT' : ocl[0], 'EE': ocl[1], 'BB': ocl[2], 'TE': ocl[3]} 

lcls = {'TT' : lcl[0], 'EE': lcl[1], 'BB': lcl[2], 'TE': lcl[3]} 


# Norm Calculator. Cannot call tau and lens in the same function call. Returns Aest as a dictionary of arrays

Aest = get_norms(['tt','eb'],lcls,ocls,rlmin,rlmax,coupling=["tau"]) # only tt and eb are implemented for tau.


# plot spectra
clf()
title("Normalization vs. ell")
xlim(2,1000)
xscale('log')
yscale('log')
for est in Aest: # plots both tt and eb 
    plot(L,Lfac_tau*Aest[est], label='L = 1/2pi, '+est) 
    
ylabel("At")
xlabel("l")
legend()
savefig('normtau.png') #save figure
