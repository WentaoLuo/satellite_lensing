import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import sys
import pandas

def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


# First get a 2d histogram of Galaxies, QSOs in the Mstar versus z bins
# Figure out the ratio of objects in each bin, and select a galaxy sample
# that mimics the QSO Mstar-z distribution.

#fig=plt.Figure((6,6))
zmin = 0.2 
zmax = 1.0
mstmin = 8.0
mstmax = 12.5
bins=21

galname = 'galaxy.dat'
df = pandas.read_csv(galname, delim_whitespace=True, header=None, names=(["ra", "dec", "z", "mst"]))
ra, dec, z, mst = df.ra.values, df.dec.values, df.z.values, df.mst.values 

idx = (z>zmin) & (z<zmax) & (mst<mstmax) & (mst>mstmin)
z = z[idx]
mst = mst[idx]
ra = ra[idx]
dec = dec[idx]

qsoname = '../mst_zdist.dat'
#qsoname = 'mst_z_i.dat'
#qsoname = 'mst_z_ii.dat'
dfqso = pandas.read_csv(qsoname, delim_whitespace=True, header=None, names=(["zqso", "mstqso"]))
zqso, mstqso = dfqso.zqso.values, dfqso.mstqso.values
mstqso = np.log10(mstqso)

idx = (zqso>zmin) & (zqso<zmax) & (mstqso<mstmax) & (mstqso>mstmin)
zqso = zqso[idx]
mstqso = mstqso[idx]

zarr = np.linspace(zmin,zmax,bins)
mstarr = np.linspace(mstmin,mstmax,bins)
histgal, xedges, yedges = np.histogram2d(z, mst, bins=[zarr, mstarr], normed="density")
histqso, xedges, yedges = np.histogram2d(zqso, mstqso, bins=[zarr, mstarr], normed="density")

ibingal = ((z-zmin)*(bins-1)/(zmax-zmin)).astype("int")
jbingal = ((mst-mstmin)*(bins-1)/(mstmax-mstmin)).astype("int")

weight = (histqso[ibingal,jbingal])/(histgal[ibingal,jbingal])
idx = np.isnan(weight) #| (np.random.random(size=weight.size)>0.10)
#idx = np.isnan(weight | (np.random.random(size=weight.size)>0.10))
#print weight
histgal_weight, xedges, yedges = np.histogram2d(z[~idx], mst[~idx], bins=[zarr, mstarr], normed="density", weights=weight[~idx])

#----New plots---------------------
fig  = plt.figure(figsize=(12,4))
grid = ImageGrid(fig,111,           # as in plt.subplot(111)
                 nrows_ncols=(1,3),
		 axes_pad=0.15,
		 share_all=True,
		 cbar_location="right",
		 cbar_mode="single",
		 cbar_size="7%",
		 cbar_pad=0.15,
		 )
#for ax in grid:
for i in range(3):
   if i ==0:
     image = histqso.T
     title = 'QSO'
   if i ==1:
     image = histgal.T
     title = 'Control'
   if i ==2:
     image = histgal_weight.T
     title = 'Weights'
   #im = grid[i].imshow(image,interpolation='nearest',origin='lower',extent=[zmin,zmax,mstmin,mstmax])
   im = grid[i].imshow(image,interpolation='nearest',origin='lower')
   #forceAspect(grid[0],aspect=0.1)
   grid[i].set_title(title)
   grid[i].set_xlim(0.2,1.0)
   grid[i].set_ylim(8.0,12.5)
#im1 = grid[0].imshow(histqso.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
#im1 = grid[0].imshow(histqso.T, origin='lower', interpolation="nearest", extent=[zmin,zmax,mstmin,mstmax])
#grid[0].set_xlim(0.2,1.0)
#grid[0].set_ylim(8.0,12.5)
"""
im2 = grid[1].imshow(histgal.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
grid[1].set_xlim(0.2,1.0)
grid[1].set_ylim(8.0,12.5)
# Colorbar
im3 = grid[2].imshow(histgal_weight.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
grid[2].set_xlim(0.2,1.0)
grid[2].set_ylim(8.0,12.5)

grid[2].cax.colorbar(im)
grid[2].cax.toggle_label(True)

fig.text(0.5,0.03,r'$\mathbf{Redshift z}$',ha='center',size=15)
fig.text(0.08,0.5,r'$\mathbf{log(M_{star})}$',va='center',size=15,rotation='vertical')
plt.show()
"""
#-END new plot------------------------------------------------------
"""
fig=pl.figure(figsize=(30,8))
ax = pl.subplot(131)
#ax.imshow(histgal.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
ax.imshow(histqso.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
ax.set_title(r'$\mathbf{QSO}$')
#pl.imshow(histgal_weight.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
#pl.xlabel(r'$\mathbf{Redshift}$',fontsize=20)
#pl.ylabel(r'$\mathbf{M_{stellar}}$',fontsize=20)

ax = pl.subplot(132)
ax.imshow(histgal.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
ax.set_title(r'$\mathbf{Control}$')
ax.set_yticks(())

ax = pl.subplot(133)
im = ax.imshow(histgal_weight.T, origin="lower", interpolation="nearest", aspect="auto",extent=[zmin,zmax,mstmin,mstmax])
ax.set_title(r'$\mathbf{Weights}$')
ax.set_yticks(())
#pl.colorbar(histgal)

fig.text(0.5,0.03,r'$\mathbf{Redshift z}$',ha='center',size=20)
fig.text(0.07,0.5,r'$\mathbf{log(M_{star})}$',va='center',size=20,rotation='vertical')
fig.colorbar(im)
pl.subplots_adjust(wspace = 0.05, hspace = 0.0 )

pl.show()


"""
idx= idx&(weight==0.0)
dt = {}
dt["0ra"] = ra[~idx]
dt["1dec"] = dec[~idx]
dt["2z"] = z[~idx]
dt["3weight"] = weight[~idx]
dt["4mst"] = mst[~idx]

df = pandas.DataFrame(data=dt)
df.to_csv("galaxy_control_reweighted_mst_full.dat", index=False, sep=" ", header=False)

#cbar = pl.colorbar()
#cbar.set_label(r'$\mathbf{Density}$',fontsize=20)
#pl.savefig("galwt_map.eps")

