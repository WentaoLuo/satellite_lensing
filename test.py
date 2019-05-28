#!/home/wtluo/anaconda/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import emcee
import hostmodule as host
#import satemodule as sate
#import starmodule as star

#-----------sub routines---------------------------
def model():
  return 0 
def lnlike():
  return 0

#-----------main------------------------------------
def main():
   import sys
   rplow    = sys.argv[1]
   rphig    = sys.argv[2]
   filename = 'camira_sat_Pmem_all'
   Ra_bcg,Dec_bcg,Ra_sat,Dec_sat,zcl,Rsep,xpos,ypos,Pmem,rich,\
   Ms_bcg,Ms_sat=np.loadtxt(filename,unpack=True,comments='#')

   fname = '../results/esds/sat_'+str(rplow)+'_'+str(rphig)+'_esd.dat'
   data  = np.loadtxt(fname,unpack=True)
   Rp    = data[7,:]
   esd   = data[5,:]
   err   = data[11,:]

   #--selection-----------------------------------------------
   ixa  = Rsep >= float(rplow)
   ixb  = Rsep <= float(rphig)
   idx  = ixa&ixb
   Roff = Rsep[idx]
   rhns = rich[idx]
   Ms   = Ms_sat[idx]
   Pm   = Pmem[idx]
   zl   = zcl[idx]
   print zl  
   #hostesd = host(Roff,rhns,Pm,zl,Rp)
if __name__=='__main__':
  main()
