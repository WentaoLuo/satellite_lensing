#!/home/wtluo/anaconda/bin/python2.7 
import numpy as np
import matplotlib.pyplot as plt
import gglens 
import emcee
import corner
# Part I model---------------------------

#---------------------------------------------------------
def model(theta,Rp):
  logMh,c,Roff,logMs,finterloper=theta
  shost  = gglens.ESD(logMh,0.0,c,Roff,0.0)
  ssub   = gglens.ESD(logMs,0.0,15,0,0)
  signal1= shost.innerNFWRoff(Rp,Roff)-shost.avrNFWRoff(Rp,Roff)
  signal2= ssub.NFWcen(Rp)
  #  
  result = signal1+signal2+finterloper*signal2

  return result

#-------------------------------------------------------------
def lnprior(theta):
  logM,con,Roff,Ms,finterloper = theta
  if 12.0<logM<20. and 1.0<con<17.0 and \
	  0.0<Roff<1.0 and 8.0<Ms<logM-0.3 and  \
	  0.0<=finterloper<=1.0:
       return 0.0
  return -np.inf
#---------------------------------------------
#def lnlike(theta,Rp,esd,covar):
def lnlike(theta,Rp,esd,err):
  logM,con,Roff,logMs,finterloper = theta
    
  esdmodel = model(theta,Rp)
  #cov  = np.dot(np.linalg.inv(covar),(esdmodel-esd))
  #chi2 = np.dot((esdmodel-esd).T,cov)
  chi2  = ((esdmodel-esd)**2/err/err)
  diff  = -0.5*(chi2)
  return diff.sum()

#-----------------------------------------------------
#def lnprob(theta,Rp,esd,covar):
def lnprob(theta,Rp,esd,err):
  lp = lnprior(theta)
  if not np.isfinite(lp):
        return -np.inf

  #return lp+lnlike(theta,Rp,esd,covar)
  return lp+lnlike(theta,Rp,esd,err)

#--------------------------------------------------------
def main():

  Rmax = 1.5
  Rmin = 0.01
  Nbin = 9 
  rbin = np.zeros(Nbin+1)
  r    = np.zeros(Nbin)
  xtmp = (np.log10(Rmax)-np.log10(Rmin))/Nbin
  for i in range(Nbin):
    ytmp1 = np.log10(0.01)+float(i)*xtmp
    ytmp2 = np.log10(0.01)+float(i+1)*xtmp
    rbin[i] = 10.0**ytmp1
    rbin[i+1] = 10.0**ytmp2
    r[i] =(rbin[i])*1./2.+(rbin[i+1])*1.0/2.0


  data = np.loadtxt('shear_sat_bin1_v3',unpack=True)
  Rp   = data[0][:]
  esd  = data[1][:]
  err  = data[2][:]
  #print esd
  covar    = np.loadtxt('cmatrix_bin1',unpack=True)
  Rr   = np.linspace(0.05,3.0,50)
  logMh= 13.6
  c    = 6.0
  Roff = 0.3
  logMs= 12.2
  finterloper =0.2
  zinterloper = np.random.normal(loc=0,scale=0.01)+zl
  pars = [logMh,c,Roff,logMs,finterloper]

  # finterloper is the fraction of false detection of satellites, which are
  # actually central galaxies from another system but with similar halo mass
  # or stellar mass along with other properties alike the ones that are real
  # satellites. So far, we simply assume a simple average halo mass identical
  # to real ones and normal redshift distribution with scatter 0.01 around
  # the clusters. And we did not take the redshift as an issue for the modelling
  # because we do not think it effects the results so much, for one thing, it is  # close to zl, for another the signal would be simply NFW.
  ndim,nwalkers = 5,100
  pos = [pars+1e-4*np.random.randn(ndim) for i in range(nwalkers)]
  # please adjust the threads according to the computer.
  sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(Rp,esd,err),threads=10)
  sampler.run_mcmc(pos,2000)

  burnin = 100
  samples=sampler.chain[:,burnin:,:].reshape((-1,ndim))
  Mh,con,Roff,Mst = map(lambda v: (v[1],v[2]-v[1],v[1]-v[0]),zip(*np.percentile(samples,[16,50,84],axis=0)))
  print 'logM: ',Mh
  print 'c: ',con
  print 'Rsig: ',Roff
  print 'Msub: ',Mst

  fig = corner.corner(samples,labels=["logM","c","Roff","Msub"],\
        truths=[Mh[0],con[0],Roff[0],Mst[0]],color="b",\
	plot_datapoints=False,plot_density=True)
  plt.savefig('submcmc_3.eps')

if __name__=='__main__':
  main()
