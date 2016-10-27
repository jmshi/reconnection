"""
map Athena++ output data to var and grid.
"""
import numpy as np
import athena_read as ath
import athena_read  as ath
from shwave import cn4, cn5
import sys
import os
#=======================================================================================

def map(trunk,key,shcoord='xy'):
  """
  input athena dataset and keyword of the variable
  return the variable data
  input:  trunk = ath.vtk(filename)
  output: 
  """
  if key == 'time': 
    t = trunk[0]
    return t
  elif key == 'grid':
    if shcoord == 'xz':
      x = trunk[1]
      y = trunk[3]
      z = trunk[2]
    else:
      x = trunk[1]
      y = trunk[2]
      z = trunk[3]
    return (x,y,z)
  elif key == 'velocity':  
    vel = trunk[4]['vel']
    return vel
  elif key == 'density':
    rho = trunk[4]['rho']
    return rho
  elif key == 'magnetic':
    bcc = trunk[4]['cc-B']
    return bcc


def rmghost(datai,ng=2):
  ndim = len(np.shape(datai))
  if ndim == 1: 
    datao = datai[ng:-ng]
  elif ndim == 3:
    ni = np.shape(datai)[2]
    nj = np.shape(datai)[1]
    nk = np.shape(datai)[0]
    if ng != 0:
      iss = ng; ie = -ng
      js = ng; je = -ng
      ks = ng; ke = -ng
      if ni < 2*ng: iss = 0; ie=ni
      if nj < 2*ng: js = 0; je=nj
      if nk < 2*ng: ks = 0; ke=nk
      datao = datai[ks:ke,js:je,iss:ie]
    else:
      datao = datai[:,:,:]
  else:
    print 'wrong dim: has to be 3-dim for data or 1-dim for grid'
    return

  return datao

