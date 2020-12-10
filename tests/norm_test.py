import numpy as np, sys
import matplotlib.pyplot as plt
import tempura as cs


# First define parameters

Tcmb  = 2.726e6    # CMB temperature
Lmax  = 3000       # maximum multipole of output normalization
lmax  = 3000
rlmin, rlmax = 100, 3000  # CMB multipole range for reconstruction
L = np.linspace(0,Lmax,Lmax+1)
l = L.copy()
Lfac = (L*(L+1.))**2/(2*np.pi)


# Load arrays of CMB unlensed and lensed Cls. Unlensed Cls are not used for now. The Cls should not be multiplied by any factors and should not have units.  

lcl = np.zeros((4,lmax+1)) # TT, EE, BB, TE
lcl[:,2:] = np.loadtxt('data/cosmo2017_10K_acc3_lensedCls.dat',unpack=True,usecols=(1,2,3,4))[:,:lmax-1] 
lcl *= 2.*np.pi / (l**2+l+1e-30) / Tcmb**2

Cpp = np.zeros((lmax+1,))
Cpp[2:] = np.loadtxt('data/cosmo2017_10K_acc3_lenspotentialCls.dat',unpack=True,usecols=(5))[:lmax-1] 
Ckk = Cpp * 2.*np.pi / 4.
Cpp[1:] *= (2.*np.pi  / (l[1:]*(l[1:]+1))**2. )


ocl = lcl*1. # observed Cl (here, no CMB noise)


# ### Compute normalization

# QDO below specifies which normalizations will be computed: TT, TE, EE, TB, EB, BB (currently BB is ignored even if you set True)

# MV is also automatically computed from these specified estimators


QDO = [True,True,True,True,True,True] # this means that TT, TE, EE, TB and EB are used for MV estimator


# "Wg" ("Wc") below is the optimal weight for constructing the MV estimator; $$\phi^{MV} = A^{MV}\sum_q W_q \phi^q$$ where $W$ is Wg (Wc) and $q = TT, TE, \dots$. 
# BB is not output, and the array has 6 normalizations (TT, TE, EE, TB, EB and MV)

Ag, Ac, Wg, Wc = cs.norm_lens.qall(QDO,Lmax,rlmin,rlmax,lcl,ocl)

# Iterative EB
ebAg,ebAc = cs.norm_lens.qeb_iter(Lmax,rlmax,rlmin,rlmax,rlmin,rlmax,lcl[1],lcl[1],lcl[2],Cpp,iter= 5,conv= 1e-3)


# plot spectra (gradient)
plt.xlim(2,Lmax)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$L$')
plt.ylabel('$C_L$')
plt.plot(L,Ckk,lw=2,color='k')
for i in range(6):
    plt.plot(L,Lfac*Ag[i],label='TT TE EE TB EB MV'.split(' ')[i])
plt.plot(L,Lfac*ebAg,ls='--',label='Iterative EB')
plt.legend()
plt.savefig('test.png')



