#!/home/wtluo/anaconda/bin/python2.7

import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt


def main():
  basePath = '..'
  snapNum  = 84
  fields = ['SubhaloMass','SubhaloSFRinRad']
  subhalos = il.groupcat.loadSubhalos(basePath,snapNum,fields=fields)
  print(subhalos.keys())
  mass_msun = subhalos['SubhaloMass'] * 1e10 / 0.704
  plt.plot(mass_msun,subhalos['SubhaloSFRinRad'],'.')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('Total Mass [$M_\odot$]')
  plt.ylabel('Star Formation Rate [$M_\odot / yr$]')
  #secPath = ''
  GroupFirstSub = il.groupcat.loadHalos(basePath,snapNum,fields=['GroupFirstSub'])
  #ptNumDm       = il.snapshot.partTypeNum('dm') 
  #ptNumGas      = il.snapshot.partTypeNum('gas') 
  #ptNumStars    = il.snapshot.partTypeNum('stars') 
  partype        = 'dm'
  haloid         = GroupFirstSub[0]
  dmpos          = il.snapshot.loadHalo(basePath,snapNum,haloid,partype,fields=None)

  x1 = dmpos['Coordinates'][:,0] # ckpc/h
  x2 = dmpos['Coordinates'][:,1] # ckpc/h
  x3 = dmpos['Coordinates'][:,2] # ckpc/h

  ptnl_p = dmpos['Potential']
  idx_ptnl_p_min = ptnl_p == ptnl_p.min()

  xcp1 = x1[idx_ptnl_p_min][0]
  xcp2 = x2[idx_ptnl_p_min][0]
  xcp3 = x3[idx_ptnl_p_min][0]

  xin1 =x1-xcp1
  xin2 =x2-xcp2
  xin3 =x3-xcp3

  plt.plot(xin1,xin2,'k.')
  plt.show()

  #dmsub = il.snapshot.loadHalo(basePath,snapNum,haloid,partype,fields=None)
  """
  for i in range(len(GroupFirstSub)):
    all_fields = il.groupcat.loadSingle(basePath,84,subhaloID=GroupFirstSub[i])
    stars_mass = all_fields['SubhaloMassInHalfRadType'][ptNumStars]
    gas_mass   = all_fields['SubhaloMassInHalfRadType'][ptNumGas]
    dm_mass    = all_fields['SubhaloMassInHalfRadType'][ptNumDm]
    frac = gas_mass / (gas_mass + stars_mass)
    print(GroupFirstSub[i], stars_mass,gas_mass,dm_mass)
  """
if __name__=='__main__':
  main()
