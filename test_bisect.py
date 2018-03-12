
"""
 driver of the whole pipe line
 for finding current sheet and null pts
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import athena4_read as ath
import scipy.ndimage.measurements as measurements
import scipy.ndimage.morphology as morphology
import cPickle as pickle
import time
import multiprocessing as mp
from itertools import product
import sys


def loadData(fname='Unstra.0080.vtk'):
  t,x,y,z,data=ath.vtk(fname)
  bx = data['cell_centered_B'][...,0]
  by = data['cell_centered_B'][...,1]
  bz = data['cell_centered_B'][...,2]
  # ---
  def curl(vx,vy,vz):
    [dzvx,dyvx,dxvx] = np.gradient(vx)
    [dzvy,dyvy,dxvy] = np.gradient(vy)
    [dzvz,dyvz,dxvz] = np.gradient(vz)
    cx = dyvz-dzvy
    cy = dzvx-dxvz
    cz = dxvy-dyvx
    # No need to del the reference by one manually
    # allow python to perform its own garbage collection
    # after the function return cxyz
    #del dzvx
    #del dzvy
    #del dzvz
    return cx,cy,cz
  # ---
  dx = dy = dz = x[1]-x[0]
  jx,jy,jz = curl(bx,by,bz)/dx
  j2 = jx**2+jy**2+jz**2
  return t,x,y,z,bx,by,bz,j2


def estimateJth(j2):
  """
  estimate the j_threshold => 25% of the total 
  ohmic dissipation having j>j_threshold;
  return jth (j_threshold^2)
  """
  j2min=np.min(j2)
  j2max=np.max(j2)
  j2th = 0.5*(j2min+j2max)
  diss4 = np.sum(j2)*0.25
  diss  = np.sum(j2[j2>=j2th])
  epsilon = 0.01
  NMAX  = 100
  count = 0
  while ((np.fabs(1.0-diss/diss4) > epsilon) or (count >=NMAX)): 
    count += 1 
    if (diss < diss4): # need to lower jth
      j2max = j2th
      j2th = 0.5*(j2min+j2th)
    else:              # need to raise jth
      j2min = j2th
      j2th = 0.5*(j2th+j2max)
    diss = np.sum(j2[j2>=j2th])
  
  print 'after ',count,' iterations: '
  print 'find threshold current sheet magnitude'
  print 'jth^2 = ',j2th,' 0.25*qdiss_tot = ',diss4,' actual qdiss = ',diss
  return j2th


t,x,y,z,bx,by,bz,j2 = loadData()
tstart = time.time()
jth = estimateJth(j2)
tend = time.time()
print 'time cost for the estimate jth: ',tend-tstart,' seconds'
