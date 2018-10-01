# This is a libs 
# for galaxy-galaxy lensing modelling.

import numpy as np

# Basic Parameters---------------------------
omega_m = 0.28
rho_crit= 9.47e-27
ckm     = 3.24078e-20
ckg     = 5.0e-31
rhoc    = rho_crit*ckg*10e+9/ckm/ckm/ckm
rhom    = rhoc*omega_m
pi      = np.pi


class ESD(object):
  def __init__(self,logM=None,Ms=None,c=None,Roff=None,Rsig=None):
     self.logM = logM
     self.Ms   = Ms
     self.c    = c
     self.Roff = Roff
     self.Rsig = Rsig

     self.r200 = (10.0**self.logM*3.0/200./rhom/pi)**(1./3.)
     self.rs   = self.r200/self.c
     self.delta= (200./3.0)*(self.c**3)\
                 /(np.log(1.0+self.c)-self.c/(1.0+self.c))
     self.amp  = 2.0*self.rs*self.delta*rhoc*1e-13
    
  # Stellar contribution of ESD----------------
  def stellar(self,Rp):
     sdens = self.Ms/(np.pi*Rp*Rp)
   
     return sdens

  # functions needs for repeatedly calculation-------
  def funcs(self,Rp):
     x   = Rp/self.rs
     x1  = x*x-1.0
     x2  = 1.0/np.sqrt(np.abs(1.0-x*x))
     x3  = np.sqrt(np.abs(1.0-x*x))
     x4  = np.log((1.0+x3)/(x))
     s1  = Rp*0.0 
     s2  = Rp*0.0 

     ixa = x>0. 
     ixb = x<1.0
     ix1 = ixa&ixb
     s1[ix1] = 1.0/x1[ix1]*(1.0-x2[ix1]*x4[ix1])
     s2[ix1] = 2.0/(x1[ix1]+1.0)*(np.log(0.5*x[ix1])\
               +x2[ix1]*x4[ix1])

     ix2 = x==1.0
     s1[ix2] = 1.0/3.0
     s2[ix2] = 2.0+2.0*np.log(0.5)

     ix3 = x>1.0
     s1[ix3] = 1.0/x1[ix3]*(1.0-x2[ix3]*np.arctan(x3[ix3]))
     s2[ix3] = 2.0/(x1[ix3]+1.0)*(np.log(0.5*x[ix3])+\
              x2[ix3]*np.arctan(x3[ix3]))

     res = {'funcf':s1,'funcg':s2}
     return res
  
  # NFW without Roff---------------------------------------
  def NFWcen(self,Rp):
     functions = self.funcs(Rp) 
     funcf     = functions['funcf']
     funcg     = functions['funcg']
     res       = self.amp*(funcf)
     
     return res
  # NFW with off center effect Roff, but not averaged over
  # directions.
  def NFWRoff(self,Rp,Rf,theta):
     cosa  = np.cos(theta*pi/180.)

     x     = np.sqrt(Rp*Rp+Rf*Rf+2.0*Rp*Rf*cosa)
     func  = self.funcs(x)
     res   = self.amp*func['funcf']

     return res
  # NFW with Roff and also averaged over directions-----
  def avrNFWRoff(self,Rp,Rf):
     nstep = 50
     vmax  = 360.0
     vmin  = 0.0
     step  = ((vmax-vmin)/nstep)
     temp  = np.zeros(len(Rp))

     for i in range(nstep):
        temp1 = self.NFWRoff(Rp,Rf,vmin+step*i)      
        temp2 = self.NFWRoff(Rp,Rf,vmin+step*(i+1)) 
        temp  = temp+0.5*step*np.pi*(temp2+temp1)/180.
  
     res = temp/2.0/np.pi          

     return res

  # Inner average density Sigma(<R)----------------------
  def innerNFWRoff(self,Rp,Rof):
     nsteps = 50
     esd    = np.zeros(len(Rp))
     for i in range(len(Rp)):
          Rf = np.linspace(0.0,Rp[i],nsteps)
          dR = Rp[i]/nsteps
          tem= 0.0
          for j in range(nsteps):
             rr = np.array([0.001,Rf[j]])
             tep= self.avrNFWRoff(rr,Rof)*Rf[j]*dR
             #print tep
             tem = tem + tep[1]
          esd[i]=2.0*tem/Rp[i]/Rp[i]

     return esd
 
