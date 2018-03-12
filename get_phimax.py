import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import my_athena_read as myath
import os.path
import sys
import cPickle as pickle
from reconnect_2d import get_path
from reconnect_2d import loadpickle
from reconnect_2d import stream_max
import reconnect_2d



if __name__=='__main__':
  if len(sys.argv)>=5:
    ts=int(sys.argv[2])
    te=int(sys.argv[3])
    nz=int(sys.argv[4])
  else:
    ts=0;te=100;nz=1

  targ = sys.argv[1]
  datapath = get_path(targ)
  tt,xx,zz,phimax=stream_max(datapath,ts,te,nz,dump=True)

  
