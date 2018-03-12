from reconnect_2d import get_path
from reconnect_2d import loadpickle
from reconnect_2d import stream_max
from reconnect_2d import plot_Jy
import reconnect_2d
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path
import os,sys


if __name__=='__main__':
  """ use this script to generate series of snapshots of Jy and vx
    and then generate animation with ffmpeg
    pay attention to the jmin and jmax, and zmin/zmax,and aspect ratio
    see plot_Jy() in reconnect_2d.py for reference
  """
  if len(sys.argv)>=4:
    ts=int(sys.argv[2])
    te=int(sys.argv[3])
  else:
    ts,te = 0,100

  targ = sys.argv[1]
  fps = 10
  if targ == 's1e4z' or targ == 's1e4z.ns' or targ == 's3e3z':
    jmin=-210;jmax=70;zmin=-0.05;zmax=0.05;aspect=2
  elif targ == 's3e4z':
    jmin=-300;jmax=100;zmin=-0.05;zmax=0.05;aspect=2
  elif targ == 's1e5z':
    jmin=-600;jmax=200;zmin=-0.025;zmax=0.025;aspect=4
  elif targ == 's3e5z':
    jmin=-1800;jmax=600;zmin=-0.025;zmax=0.025;aspect=4
  elif targ == 's1e6z':
    jmin=-1800;jmax=600;zmin=-0.025;zmax=0.025;aspect=4
  else:
    jmin=-1800;jmax=600;zmin=-0.025;zmax=0.025;aspect=4

  datapath = get_path(targ)
  for iframe in np.arange(ts,te+1):
    fname = datapath+'/tmp/Jy-'+targ+'.'+str(iframe).zfill(5)+'.png'
    plot_Jy(targ,iframe=iframe,jmin=jmin,jmax=jmax,zmin=zmin,zmax=zmax,aspect=aspect,filename=fname)
  
  moviename = datapath+'/tmp/Jy-'+targ+'.mp4'
  os.system('rm '+moviename)
  os.system('ffmpeg -r '+str(fps)+' -i '+ datapath+'/tmp/Jy-'+targ+'.%05d.png  '+moviename)
  
