import numpy as np
import camb
from camb import model, initialpower
import mcfit
from scipy import integrate
from optparse import OptionParser
import warnings
warnings.filterwarnings('ignore')

#-----------------------------------------------------------------------------
def Pklin(redshift,cosmology):
  z    = redshift
  H_null,Ombh2,Omch2,nps=cosmology
  pars = camb.CAMBparams()
  pars.set_cosmology(H0=H_null,ombh2=Ombh2,omch2=Omch2)
  pars.set_dark_energy()
  pars.InitPower.set_params(ns=nps)
  pars.set_matter_power(redshifts=[0.0,z],kmax=2.0)

  #Linear spectra--------
  pars.NonLinear=model.NonLinear_none
  results       =camb.get_results(pars)
  kh,znon,pk    =results.get_matter_power_spectrum(minkh=1e-4,maxkh=1000.0,npoints=1024)
  xi            =mcfit.P2xi(kh,l=0)
  #Calculate corr from PS------
  nxx   = 50
  Rmin  = -2.0
  Rmax  = 2.0
  rr    = np.logspace(Rmin,Rmax,nxx)
  rx,corrfunc   =xi(pk,extrap=True)

  #print np.shape(kh),np.shape(pk)  
  return {'k':kh,'pklin':pk[1,:],'r':rx,'corr':corrfunc[0,:]}
#---------------------------------------------------------------------------------------

def galaxybias(logM):
  Mnl = 8.73*10e+12
  Mh  = 10.0**logM
  xx  = Mh/Mnl
  b0  = 0.53+0.39*xx**0.45+(0.13/(40.0*xx+1.0))\
       +5.0*0.0004*xx**1.5
  bias= b0+\
        np.log10(xx)*(0.4*(omega_m-0.3+nps-1)+\
        0.3*(sigma8-0.9+h-0.7)+0.8*alphas)
  return bias

#--------------------------------------------------------------------------------------------

z  = 0.1
cosmology  = [100,0.022,0.122,0.965]
results    = Pklin(z,cosmology)
k          = results['k']
pklin      = results['pklin']
r          = results['r']
corr       = results['corr']
print np.shape(r),np.shape(corr)  
import matplotlib.pyplot as plt
plt.plot(r,corr,'k-')
plt.xscale('log')
plt.yscale('log')
plt.show()

