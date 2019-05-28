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
   ixc  = rich >=15
   idx  = ixa&ixb&ixc
   Roff = Rsep[idx]
   rhns = rich[idx]
   Ms   = Ms_sat[idx]
   Pm   = Pmem[idx]
   zl   = zcl[idx]
   zlunq= np.unique(zl)
   print len(zl)
   print len(zlunq)
   #------test rich mass integration--------------------------
   """
   This test is designed to test the accuracy of the integration while
   using a very simple shit.
   """
   """
   testPm1   = host.richmass(15,100)
   testPm2   = host.richmass(25,100)
   testPm3   = host.richmass(55,100)
   testPm4   = host.richmass(150,100)
   testPm5   = host.richmass(500,100)
   plt.plot(testPm1['Mrange'],testPm1['Prob'],'r-',linewidth=3,label=r'$\lambda$=15')
   plt.plot(testPm2['Mrange'],testPm2['Prob'],'g-',linewidth=3,label=r'$\lambda$=25')
   plt.plot(testPm3['Mrange'],testPm3['Prob'],'b-',linewidth=3,label=r'$\lambda$=55')
   plt.plot(testPm4['Mrange'],testPm4['Prob'],'y-',linewidth=3,label=r'$\lambda$=150')
   plt.plot(testPm5['Mrange'],testPm5['Prob'],'c-',linewidth=3,label=r'$\lambda$=500')
   plt.xlabel(r'logMh',fontsize=15)
   plt.xlim(12,17)
   plt.ylabel(r'P(logM|$\lambda$)',fontsize=15)
   plt.legend()
   #plt.savefig('massrich.eps')
   plt.show()
   """
   hostesd = host.hostRsep(Roff,rhns,Pm,zl,Rp)
   plt.plot(Rp,hostesd,'k-',linewidth=3)
   plt.errorbar(Rp,esd,yerr=err,fmt='k.',ms=10,elinewidth=3)
   plt.xscale('log')
   plt.xlabel(r'Rp Mpc/h',fontsize=15)
   plt.ylabel(r'ESD',fontsize=15)
   plt.show()
if __name__=='__main__':
  main()
