#!/usr/bin/env python
"""
scripts to manipulate shwave data
"""

#
import numpy as np
import matplotlib.pyplot as plt
import athena_read as ath
from astropy.io import ascii
from astropy.table import Table, Column
import sys
import os

def cn4(number):
  snumber = str(number)
  nlen = len(snumber)
  if nlen < 2:
    snumber = '000'+snumber
  elif nlen < 3:
    snumber = '00'+snumber
  elif nlen < 4:
    snumber = '0'+snumber
  else:
    snumber = snumber
  return snumber

def cn5(number):
  snumber = str(number)
  nlen = len(snumber)
  if nlen < 2:
    snumber = '0000'+snumber
  elif nlen < 3:
    snumber = '000'+snumber
  elif nlen < 4:
    snumber = '00'+snumber
  elif nlen < 5:
    snumber = '0'+snumber
  else:
    snumber = snumber
  return snumber


def vx_amp(targ,ts,te):
  """ read filename.vtk from frame=ts to frame=te """
  """ and plot the vx amplitude over time         """
  nt = te-ts+1
  time = np.empty([nt])
  vx   = np.empty([nt])
  for i in range(te-ts+1):
    filename = targ+'/sst.block0.out2.'+cn5(ts+i)+'.vtk'
    trunk = ath.vtk(filename)
    time[i] = trunk[0]
    absvx = trunk[4]['vel'][:,:,:,0]
    vx[i] = np.amax(np.sqrt(np.square(absvx)*trunk[4]['rho']))
  torb = 2.0*np.pi/1e-3
  vx  = vx*1e3
  time = time/torb
  print type(vx)
  ascii.write(Table([time,vx],names=['#time(orbits)','|v_x|(c_s)']),targ+'/absvx.tab')
  return 



if __name__=='__main__':
    if len( sys.argv ) < 4:
        print "Please specify input file dir" \
            + " and ts, te"
        exit(  );
    #

    fdir = sys.argv[ 1 ];
    ts = int(float(sys.argv[ 2 ]));
    te = int(float(sys.argv[ 3 ]));

    vx_amp(fdir,ts,te);
    #shwave.vx_amp('../develop/bin/debug',0,400)
#

