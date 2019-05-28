import numpy as np
import haloparams as halos
import matplotlib.pyplot as plt

#---read interpolation tables------------------------
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
   Rps,dsig = np.loadtxt(fsub,dtype=np.str,\
              usecols=(0,1),unpack=True)
   for j in range(len(Rps)):
      tabs[i,j,0] = Rps[j]
      tabs[i,j,1] = dsig[j]

      tmx  = np.linspace(1,ny,ny)
      tmx  = tmx.astype(np.int)
      rtmp = tabs[0,:,0]
#------------------------------------------------------------------

def richmass(rich):

  return 0

def offcenter(amp,rs,Roff,Rp,sigr):
  nrp     = 300
  npi     = 100
  rprob   = np.linspace(0.0,1.0,nrp)
  protb   = (rprob/sigr/sigr)*np.exp(-0.5*(rprob*rprob)/sigr/sigr)
  totprb  = protb.sum()
  summ    = 0.0
  probp   = 2.0*pi/npi
  for i in range(nrp):
    for j in range(npi):
      rxx     = np.sqrt(Roff*Roff+rprob[i]*rprob[i]\
              -2.0*np.cos(float(j)*probp)*Roff*rprob[i])
      rr      = rxx/rs
      xx      = np.linspace(1,nx,nx)-1

      idx     = np.abs(rc2-rr)==np.min(np.abs(rc2-rr))
      inx     = int(xx[idx])
      tmp     = float(inx)*amp*tabs[inx,:,1]*(rc2[inx]-rc1[inx])/rc2[inx]
      summ    = summ+(protb[i]/totprb)*float(1.0/npi)*tmp
  res     = np.interp(Rp,rs*rtmp,summ)
  return res
 
def hostRsep(Roff,rich,Pm,zl,Rp):
  hostsub = 0
  Pcen    = 0.68

  #fo
  #amp,rs,r200 = haloparams(logMh,hostsub,zl)

  return 0 
