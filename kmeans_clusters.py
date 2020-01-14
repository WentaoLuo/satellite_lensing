#!/home/wtluo/anaconda3/bin/python3.7

import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt


x,y,z = np.loadtxt('particles.dat',unpack=True)
pos   = np.array([x,y]).T
kmeans= KMeans(n_clusters=8,random_state=0,max_iter=2000).fit(pos)
cntrs = kmeans.cluster_centers_
xc    = cntrs[:,0]
yc    = cntrs[:,1]

print(cntrs) 
idx   = np.random.randint(low=0,high=84189,size=20000)
plt.plot(x[idx],y[idx],'k.')
plt.plot(xc,yc,'r*',ms=15)
plt.show()

