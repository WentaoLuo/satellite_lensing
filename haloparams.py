import numpy as np

h       = 0.732
w       = -1.0
omega_m = 0.31
omega_l = 0.69
omega_k = 1.0-omega_m-omega_l
rho_crt0= 2.78e11                # M_sun Mpc^-3 *h*h 
rho_bar0= rho_crt0*omega_m       # M_sun Mpc^-3 *h*h
pi      = np.pi
ns      = 0.95
alphas  = -0.04
sigma8  = 0.815

def haloparams(logM,hostsub,zl):
   efunc  = 1.0/np.sqrt(omega_m*(1.0+zl)**3+\
             omega_l*(1.0+zl)**(3*(1.0+w))+\
             omega_k*(1.0+zl)**2)
   rhoc      = rho_crt0/efunc/efunc
   omegmz    = omega_m*(1.0+zl)**3*efunc**2
   ov        = 1.0/omegmz-1.0
   dv        = 18.8*pi*pi*(1.0+0.4093*ov**0.9052)
   rhom      = rhoc*omegmz
   if hostsub ==0:
       con  = 4.67*(10.0**(logM-14)*h)**(-0.11) # Neto et al 2007 
       r200 = (10.0**logM*3.0/200./4.0/rhom/pi)**(1./3.)
       rs   = r200/con
       rtrc = r200
   if hostsub ==1:
       con  = 4.67*(10.0**(logM-14)*h)**(-0.11) # Neto et al 2007 
       r200 = (10.0**logM*3.0/200./4.0/rhom/pi)**(1./3.)
       rs   = r200/con

   delta= (200./3.0)*(con**3)\
          /(np.log(1.0+con)-con/(1.0+con))

   amp  = 2.0*rs*delta*rhoc*10e-12
   res  = np.array([amp,rs,r200])

   return res

