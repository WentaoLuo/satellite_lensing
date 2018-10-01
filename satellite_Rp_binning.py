#!/home/wtluo/anaconda/bin/python2.7

import numpy as np
import sys
import matplotlib.pyplot as plt
import mycosmology as cosmos
from scipy.spatial import KDTree

#---------------------------------------------
pi = np.pi
vc = 2.9970e5
G  = 6.754e-11
ckm= 3.240779e-17
ckg= 5.027854e-31
fac= (vc*vc*ckg)/(4.0*pi*G*ckm)

#-----------Functions for shear--------------------------
def angulardis2(z1,z2):
   dis   = np.zeros(len(z2))
   for i in range(len(z2)):
       dls  = cosmos.Da2(z1,z2[i])
       dis[i]  = dls
   return dis

def readsource(filename,iredshift):
  ra,dec,e1,e2,res,wt,m,c1,c2,erms,z1,\
  z2,z3,z4,z5,z6=np.loadtxt(filename,unpack=True)
  if iredshift ==1:
      zs = z1
  if iredshift ==2:
      zs = z2
  if iredshift ==3:
      zs = z3
  if iredshift ==4:
      zs = z4
  if iredshift ==5:
      zs = z5
  if iredshift ==6:
      zs = z6
  dis    = np.zeros(len(ra))
  for i in range(len(ra)):
     dis[i]=cosmos.Da(zs[i])
  
  ix       = res>1.0/3.0
  shapes   = {"ra":ra[ix],"dec":dec[ix],"z":zs[ix],\
              "dis":dis[ix],"e1":e1,"e2":e2,\
              "erms":erms,"res":res,"m":m,"c1":c1,"c2":c2}
  return shapes

#-----------Functions for clusters------------------------
def separation(ras,decs,ral,decl):
  a1=np.cos(pi/2.-decl*pi/180.)
  a2=np.cos(pi/2.-decs*pi/180.)
  a3=np.sin(pi/2.-decl*pi/180.)
  a4=np.sin(pi/2.-decs*pi/180.)
  a5=np.cos((ras-ral)*pi/180.)

  thet=np.arccos(a1*a2+a3*a4*a5)

  return thet

def phi(x,y,xc,yc,pm):
  n  =np.size(pm)
  qxx=0.0
  qyy=0.0
  qxy=0.0
  T  =0.0
  for i in range(n):
     if ~np.isnan(x[i]):
       qxx=qxx+(x[i]-xc)*(x[i]-xc)*pm[i]
       qyy=qyy+(y[i]-yc)*(y[i]-yc)*pm[i]
       qxy=qxy+(x[i]-xc)*(y[i]-yc)*pm[i]
       T  =T+pm[i]

       qxx=qxx/T
       qyy=qyy/T
       qxy=qxy/T

       theta=0.5*np.arctan(qxy/(qxx-qyy))
       e1   =(qxx-qyy)/(qxx+qyy)
       e2   =(qxy)/(qxx+qyy)
  return [e1,e2]
#------------------------------------------------------------------------------

def main():
  fclusters   = sys.argv[1] 
  fsatellites = sys.argv[2] 
  #print fclusters,fsatellites
  ra,dec,z,rich,Ms,Izspec=np.loadtxt(fclusters,unpack=True,skiprows=1)
  
  xx1,xx2,xx3,xx4,ras,decs,mss,Pmem=np.loadtxt(fsatellites,unpack=True,skiprows=1)

  for i in range(len(ra)):
    if rich[i]>=10.0 and rich[i]<30000.0 and z[i]>=0.55 and z[i]<1.10:
      rax = ra[i]
      decx= dec[i]
      zx  = z[i]
      dl  = cosmos.Da(zx)
      ix  = xx1==rax
      iy  = Pmem>=0.8
      idx = ix&iy
      rr1 = ras[idx]
      rr2 = decs[idx]
      for j in range(len(xx1[idx])): 
          sep = dl*separation(rax,decx,rr1[j],rr2[j])*(1.0+zx)
          if sep >=0.1 and sep<=0.2:
             print rr1[j],rr2[j],zx,Pmem[j]
    

   
if __name__=='__main__':
  main()
