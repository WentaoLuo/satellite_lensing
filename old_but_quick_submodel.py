#!/home/wtluo/anaconda/bin/python2.7 
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
import scipy.optimize as opt
from scipy import interpolate
import sys
#----Basic Parameters---------------------------
h       = 1.0
w       = -1.0
omega_m = 0.315
omega_l = 0.72
omega_k = 1.0-omega_m-omega_l
rho_crt0= 2.78e11                # M_sun Mpc^-3 *h*h 
rho_bar0= rho_crt0*omega_m       # M_sun Mpc^-3 *h*h
pi      = np.pi
ns      = 0.95
alphas  = -0.04
sigma8  = 0.815
#-----------------------------------------------------------
fname,Rc1,Rc2 = np.loadtxt('r_c-range.tsv',\
                dtype=np.str,usecols=(0,1,2),\
		                unpack=True)
nx  = len(fname)
ny  = 3077
rs  = 1
rc1 = np.zeros(nx)
rc2 = np.zeros(nx)
tabs= np.zeros((nx,ny,2))
for i in range(nx):
  rc1[i] = float(Rc1[i])
  rc2[i] = float(Rc2[i])
  fsub   = 'mean/'+fname[i]
  Rps,dsig = np.loadtxt(fsub,dtype=np.str,usecols=(0,1),unpack=True)
  for j in range(len(Rps)):
     tabs[i,j,0] = Rps[j]
     tabs[i,j,1] = dsig[j]

tmx  = np.linspace(1,ny,ny)
tmx  = tmx.astype(np.int)
rtmp = tabs[0,:,0]

#-------------------------------------------------------------------------------------
def haloparams(logM,con,hostsub,frac,zl):
   efunc     = 1.0/np.sqrt(omega_m*(1.0+zl)**3+\
               omega_l*(1.0+zl)**(3*(1.0+w))+\
               omega_k*(1.0+zl)**2)
   rhoc      = rho_crt0/efunc/efunc
   omegmz    = omega_m*(1.0+zl)**3*efunc**2
   ov        = 1.0/omegmz-1.0
   dv        = 18.8*pi*pi*(1.0+0.4093*ov**0.9052)
   rhom      = rhoc*omegmz
   if hostsub ==0:
     r200 = (10.0**logM*3.0/200./4.0/rhom/pi)**(1./3.)
     rs   = r200/con
     rtrc = r200
   if hostsub ==1:
     r200 = (10.0**logM*3.0/200./4.0/rhom/pi)**(1./3.)
     rtrc = frac*r200
     rs   = rtrc/con

   delta= (200./3.0)*(con**3)\
          /(np.log(1.0+con)-con/(1.0+con))

   amp  = 2.0*rs*delta*rhoc*10e-14
   res  = np.array([amp,rs,r200,rtrc])

   return res

def nfwfuncs(rsub,Rp): 
   x   = Rp/rsub
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

   res = s2-s1
   #print res
   return res

#------------------
def subhalo(theta,logMh,ch,Roff,zl,Rp):
  #logMh, host halomass; ch, host halo concentration from Ryuma's fitting mass-richness relation
  #logMsub,sub halo mass; ratio: between truncated radius and r200 of subhalo
  #Roff, separation between subhalo and host halo
  #logMinter intervening subhalo mass
  logMsub,ratio,logMinter = theta
  con     = 10.0
  hostsub = 0
  amp,rs,r200,rtc= haloparams(logMsub,con,hostsub,ratio,zl)
  esd1    = amp*nfwfuncs(rs,Rp)
  esd2    = amp*nfwfuncs(rs,Rp)
  idx     = Rp>=rtc
  #esd2[idx] = 0.0
  return {'notruncation':esd1,'truncation':esd2}
def hostRsep(logMh,ch,Roff,zl,Rp):
  hostsub = 0
  ratio   = 1
  amp,rs,r200,rtrc = haloparams(logMh,ch,hostsub,ratio,zl)
  rr      = Roff/rs
  res     = np.zeros(len(Rp))
  xx      = np.linspace(1,nx,nx)-1
  idx     = np.abs(rc2-rr)==np.min(np.abs(rc2-rr))
  inx     = int(xx[idx])
  summ    = np.zeros(ny)
  tmp     = 34.0*amp*tabs[inx,:,1]*(rc2[inx]-rc1[inx])/rc2[inx]
  res     = np.interp(Rp,rs*rtmp,tmp)
  return res
def intervener():
  return 0
#-------------------
def lnlike():
  return 0
def lnprior():
  return 0
def lnprob():
  return 0
#--------------------------------------------------------------------------------------
def main():
  fnn      = '../results/esds/sat_0.7_0.8_esd.dat'
  results  = np.loadtxt(fnn,unpack=True,comments='#')
  esd      = results[5,:] 
  rp       = results[7,:]
  error    = results[11,:]

  plt.errorbar(rp,esd,yerr=error,fmt='k.',ms=15,elinewidth=2.5)
  plt.plot([0.75,0.75],[-5,150],'r-',linewidth=3)
  plt.xlabel('R Mpc/h')
  plt.title('Rp(0.7,0.8)')
  plt.ylabel(r'$ESD(hM_{}/pc^2)$')
  plt.xlim(0.02,5.0)
  plt.xscale('log')
  plt.ylim(-5,100.)
  plt.legend()
  plt.show()
  filename = 'camira_sat_Pmem_all'
  Ra_bcg,Dec_bcg,Ra_sat,Dec_sat,zcl,Rsep,xpos,ypos,Pmem,rich,\
       Ms_bcg,Ms_sat=np.loadtxt(filename,unpack=True,comments='#')
  #--selection----------------------
  ixa = Rsep >=0.7
  ixb = Rsep <0.8
  idx = ixa&ixb
  rac = Ra_bcg[idx]  
  decc= Dec_bcg[idx]  
  rab = np.unique(rac)
  nxx = len(rab)
  xitv= np.array([])
  yitv= np.array([])
  Ritv= np.array([])
  for i in range(nxx):
     iax = Ra_bcg==rab[i]
     ibx = Rsep < 0.7
     ixx = iax&ibx
     xitv=np.append(xitv,[xpos[ixx]])
     yitv=np.append(yitv,[ypos[ixx]])
     Ritv=np.append(Ritv,[Rsep[ixx]])

  #plt.plot(xpos[idx],ypos[idx],'r.')
  #plt.plot(xitv,yitv,'b.')
  #plt.show()
  logMh = 14.2
  ch    = 5.0
  logMs = 12.5
  logMit= 14.2
  ratio = 1.0
  Roff  = 0.1
  zl    = 0.1
  npts  = 200
  Rpair = 0.5
  fpair = 0.9
  Rp    = np.logspace(-2,1,npts)
  msesd = 1.5*(np.mean(10.0**Ms_sat[idx]))/1.0e+12/pi/Rp/Rp
  print Rp[3],msesd[3]
  Rsp   = Rsep[idx]
  esdmh0  = np.zeros(npts) 
  esdmh1  = np.zeros(npts) 
  esdmh2  = np.zeros(npts) 
  esdmh3  = np.zeros(npts) 
  esdms  = np.zeros(npts) 
  esditv= np.zeros(npts) 
  Rsig1  = np.random.normal(loc=0.0,scale=0.08,size=len(Rsp))
  Rsig2  = np.random.normal(loc=0.0,scale=0.2,size=len(Rsp))
  Rsig3  = np.random.normal(loc=0.0,scale=0.25,size=len(Rsp))
  for i in range(len(Rsp)):
     tmp0   = hostRsep(logMh,ch,Rsp[i],zl,Rp)
     tmp1   = hostRsep(logMh,ch,Rsp[i]+Rsig1[i],zl,Rp)
     tmp2   = hostRsep(logMh,ch,Rsp[i]+Rsig2[i],zl,Rp)
     tmp3   = hostRsep(logMh,ch,Rsp[i]+Rsig3[i],zl,Rp)
     #tmp2   = hostRsep(logMs,ch,0.2,zl,Rp)
     esdmh0  = esdmh0+tmp0
     esdmh1  = esdmh1+tmp1
     esdmh2  = esdmh2+tmp2
     esdmh3  = esdmh3+tmp3
     #esdms  = esdms+tmp2
  for j in range(len(Rsp)):

     #rss   = np.mean(Rsp)-Rpair-Rsig[j]
     rss   = np.mean(Rsp)-Rpair-Rsig1[j]
     tmp   = hostRsep(logMit,15,rss,zl,Rp)
     esditv  = esditv+tmp
  
  substr= subhalo([logMs,ratio,0.0],logMh,ch,Roff,zl,Rp)
  subesd= substr['notruncation']
  subtrc= substr['truncation']
  esdms  = esdms/float(len(Rsp))
  esdmh0  = esdmh0/float(len(Rsp))
  esdmh1  = esdmh1/float(len(Rsp))
  esdmh2  = esdmh2/float(len(Rsp))
  esdmh3  = esdmh3/float(len(Rsp))
  esditv= fpair*esditv/float(len(Ritv))
  plt.plot(Rp,msesd,'--',color='cyan',linewidth=2,label='Ms_sat')
  #plt.plot(Rp,esdms,'r--',linewidth=2,label='Host halo')
  plt.plot(Rp,esdmh0,'r--',linewidth=2,label='Host halo Rsig=0.0')
  plt.plot(Rp,esdmh1,'r--',linewidth=2,label='Host halo Rsig=0.08')
  plt.plot(Rp,esdmh2,'r--',linewidth=2,label='Host halo Rsig=0.2')
  plt.plot(Rp,esdmh3,'r--',linewidth=2,label='Host halo Rsig=0.25')
  plt.plot(Rp,subtrc,'b--',linewidth=2,label='Subhalo')
  #plt.plot(Rp,esditv,'g--',linewidth=2,label='Intervening structrue')
  #plt.plot(Rp,esdm+subtrc+msesd+subesd,'k-',linewidth=3,label='All')
  plt.plot(Rp,esdmh0+subtrc+msesd,'k-',linewidth=3,label='All Rsig=0.0')
  #plt.plot(Rp,esdmh1+subtrc+msesd,'k-',linewidth=3,label='All Rsig=0.08')
  #plt.plot(Rp,esdmh2+subtrc+msesd,'k-',linewidth=3,label='All Rsig=0.2')
  #plt.plot(Rp,esdmh3+subtrc+msesd,'k-',linewidth=3,label='All Rsig=0.25')
  #plt.plot(Rp,esdm+subesd+esditv,'y-',linewidth=3)
  plt.errorbar(rp,esd,yerr=error,fmt='k.',ms=25,elinewidth=3.5)
  plt.xlim(0.02,2.0)
  plt.ylim(-20.0,150.0)
  plt.xscale('log')
  plt.xlabel('R Mpc/h')
  plt.ylabel(r'$ESD(hM_{}/pc^2)$')
  plt.legend()
  plt.savefig('modeling.eps')
  plt.show()
  #------------------------------------------------------------------ 
  #logMh = 14.0
  #ch    = 5.0
  #logMs = 10.0
  #ratio = 1.0
  #Roff  = 0.1
  #zl    = 0.1
  
  Rp    = np.linspace(0.01,5,200)
  substr= subhalo([logMs,ratio,0.0],logMh,ch,Roff,zl,Rp)
  subesd= substr['notruncation']
  subtrc= substr['truncation']
  esd1  = hostRsep(logMh,ch,0.1,zl,Rp)
  esd11 = hostRsep(logMh,ch,0.14,zl,Rp)
  esd12 = hostRsep(logMh,ch,0.17,zl,Rp)
  esd2  = hostRsep(logMh,ch,0.2,zl,Rp)
  esd21 = hostRsep(logMh,ch,0.24,zl,Rp)
  esd22 = hostRsep(logMh,ch,0.27,zl,Rp)
  esd3  = hostRsep(logMh,ch,0.3,zl,Rp)
  esd31 = hostRsep(logMh,ch,0.34,zl,Rp)
  esd32 = hostRsep(logMh,ch,0.37,zl,Rp)
  esd4  = hostRsep(logMh,ch,0.4,zl,Rp)
  esd41 = hostRsep(logMh,ch,0.44,zl,Rp)
  esd42 = hostRsep(logMh,ch,0.47,zl,Rp)
  esd5  = hostRsep(logMh,ch,0.5,zl,Rp)
  esd51 = hostRsep(logMh,ch,0.54,zl,Rp)
  esd52 = hostRsep(logMh,ch,0.57,zl,Rp)
  plt.plot(Rp,esd1)
  plt.plot(Rp,esd11)
  plt.plot(Rp,esd12)
  plt.plot(Rp,esd2)
  plt.plot(Rp,esd21)
  plt.plot(Rp,esd22)
  plt.plot(Rp,esd3)
  plt.plot(Rp,esd31)
  plt.plot(Rp,esd32)
  plt.plot(Rp,esd4)
  plt.plot(Rp,esd41)
  plt.plot(Rp,esd42)
  plt.plot(Rp,esd5)
  plt.plot(Rp,esd51)
  plt.plot(Rp,esd52)
  plt.plot(Rp,subesd,'k--',linewidth=4)
  plt.plot(Rp,subtrc,'b-.',linewidth=4)
  plt.plot(Rp,(esd1+esd11+esd12+esd2+esd21+esd22+esd3+esd31+esd32+esd4+esd41+esd42+esd5+esd51+esd52)/15.0,'k.',linewidth=4)
  plt.plot(Rp,subesd+(esd1+esd11+esd12+esd2+esd21+esd22+esd3+esd31+esd32+esd4+esd41+esd42+esd5+esd51+esd52)/15.0,'r-',linewidth=4)
  plt.plot(Rp,subtrc+(esd1+esd11+esd12+esd2+esd21+esd22+esd3+esd31+esd32+esd4+esd41+esd42+esd5+esd51+esd52)/15.0,'g-',linewidth=4)
  plt.xscale('log')
  plt.xlim(0.01,2.0)
  plt.xlabel('R Mpc/h')
  plt.ylabel(r'$ESD(hM_{}/pc^2)$')

  #plt.yscale('log')
  plt.show()
  
  #----------------------------------------------------
if __name__=='__main__':
  main()



