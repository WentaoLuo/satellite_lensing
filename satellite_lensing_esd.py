#!/home/wtluo/anaconda/bin/python2.7

#----README--------------------------------------
# This is a small code to convert 2D dark matter
# particle distribution to 2D ESD from simulation.
# This demo only plots the profile given a sample
# halo from illustris simulation snap128-halo50.
# You can also do the stacking bby modify this version.
# Have fun!
#-------------------------------------------------
import numpy as np
import illustris_python as illus
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3d

pi =np.pi

def dens2esd(x,y,rb,nbin,mass):
  dis = np.sqrt(x*x+y*y) 
  esd = np.zeros(nbin)
  for i in range(nbin):
      ixa = dis >= rb[i]
      ixb = dis <= rb[i+1]
      iy  = dis <= (rb[i]/2.0+rb[i+1]/2.0)
      ix  = ixa&ixb
      denin = mass*len(x[iy])/2.0/pi/(rb[i]/2.0+rb[i+1]/2.0)**2
      annul = 2.0*pi*(rb[i+1]*rb[i+1]-rb[i]*rb[i])
      denat = mass*len(x[ix])/annul
      esd[i]= denin-denat
  return esd/10e+12  # 10e+12 factor comes from Mpc^2 to pc^2
def main():
 SnapID = 128
 HaloID = 51 

 mass   = 10e+9  # particle mass 
 baseDir= '../'
 dm     = illus.snapshot.loadHalo(baseDir,SnapID,HaloID,'dm') 
 npp    = dm["count"]
 x1     = dm["Coordinates"][:,0] # kpc/h
 x2     = dm["Coordinates"][:,1] # kpc/h
 x3     = dm["Coordinates"][:,2] # kpc/h
# for i in range(len(x1)):
#    print x1[i]/1000.,x2[i]/1000.,x3[i]/1000.
 xc1    = np.mean(x1)
 xc2    = np.mean(x2)
 xc3    = np.mean(x3)

 xs1    = 67.18
 ys1    = 6.05
 xs2    = 66.94
 ys2    = 5.20
 xs3    = 66.81
 ys3    = 4.93
#---plots to check if everything is in control-------- 
 #print xc1,xc2,xc3
 #fig    = plt.figure()
 #axs    = p3d.Axes3D(fig)
 #axs.plot((x1-xc1)/1000.,(x2-xc2)/1000.,(x3-xc3)/1000.,'k.')
 #axs.set_xlim3d([-0.5,0.5])
 #axs.set_xlabel('X kpc/h')
 #axs.set_ylim3d([-0.5,0.5])
 #axs.set_xlabel('Y kpc/h')
 #axs.set_zlim3d([-0.5,0.5])
 #axs.set_xlabel('Z kpc/h')
 #plt.plot((x1-xc1)/1000.,(x2-xc2)/1000.,'k.')
 #plt.xlim([-2,2])
 #plt.xlabel('X kpc/h')
 #plt.ylim([-2,2])
 #plt.ylabel('Y kpc/h')
 plt.show()
#--------------------------------------------------------
#---starts to calculate the ESD--------------------------
 Rmax = 5.0
 Rmin = 0.02
 Nbin = 20
 rbin = np.zeros(Nbin+1)
 Rp   = np.zeros(Nbin)
 xtmp = (np.log10(Rmax)-np.log10(Rmin))/Nbin
 for i in range(Nbin):
     ytmp1 = np.log10(0.02)+float(i)*xtmp
     ytmp2 = np.log10(0.02)+float(i+1)*xtmp
     rbin[i] = 10.0**ytmp1
     rbin[i+1] = 10.0**ytmp2
     Rp[i] =(rbin[i]+rbin[i+1])/2.0


 ESD = dens2esd((x1-xs3*1000.0)/1000.0,(x2-ys3*1000.0)/1000.0,rbin,Nbin,mass)
#---plot check if it make sense----------------------------
 plt.plot(Rp,ESD,'r-',linewidth=2.5)
 plt.xlabel('Rp kpc')
 plt.ylabel('ESD ')
 plt.xlim(0.025,5.)
 plt.ylim(-20.,150)
 plt.xscale('log')
 #plt.yscale('log')
 plt.show()
#---END----------------------------------------------------
if __name__=='__main__':
   main()
