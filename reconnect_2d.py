import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import my_athena_read as myath
import os.path
import sys
import cPickle as pickle

def get_path(targname='s1e4'):
  direname='/tigress/jiming/current_sheet/bin/'
  targname=targname+'/'
  return direname+targname

def loaddata(datapath,iframe,level=4,subsample=True,quantities=['Bcc1','Bcc2','vel1'],basename='CS.out2',zrange=0.05,appdname='athdf'):
  fname=datapath+'/'+basename+'.'+str(iframe).zfill(5)+'.'+appdname 
  print "load data from "+fname
  time,data = myath.athdf(fname,level=level,subsample=subsample,quantities=quantities)
  x = data['x1f'];z = data['x2f']
  index = np.where(np.fabs(z)<=zrange)
  zz = z[index]
  xx = x
  bx = data['Bcc1'][0,index,:]
  bx = bx[0,...]
  bz = data['Bcc2'][0,index,:]
  bz = bz[0,...]
  vx = data['vel1'][0,index,:]
  vx = vx[0,...]
  return time,xx,zz,bx,bz,vx

def pickledata(time,xx,zz,bx,bz,vx,datapath,iframe):
  tmppath = datapath+'/tmp/'
  if not os.path.exists(tmppath):
    os.makedirs(tmppath)
  pname = tmppath+'cs.'+str(iframe).zfill(5)+'.p'
  pickle.dump([time,xx,zz,bx,bz,vx], open(pname, "wb" ),2)

def loadpickle(datapath,iframe):
  tmppath = datapath+'/tmp/'
  pname = tmppath+'cs.'+str(iframe).zfill(5)+'.p'
  time,xx,zz,bx,bz,vx  = pickle.load( open(pname, "rb" ) )
  return float(time),np.array(xx),np.array(zz),\
	np.array(bx),np.array(bz),np.array(vx)

def stream_max(datapath,ts=0,te=50,nz=11,dump=False):
  nt = te-ts+1
  phimax=np.zeros([nt,nz])
  tt =np.zeros(nt)
  for iframe in np.arange(ts,te+1):
    time,xx,zz,bx,bz,vx = loadpickle(datapath,iframe)
    tt[iframe-ts]=time
    iz0 = len(zz)/2; dx = xx[1]-xx[0]
    for iz in np.arange(0,nz,1):
      if(zz[0]==0.0):
        phimax[iframe-ts,iz]=np.max(np.cumsum(bz[iz,:]))*dx
      else:
        phimax[iframe-ts,iz]=np.max(np.cumsum(bz[iz0+iz-nz/2,:]))*dx
  
  if dump:
    tmppath = datapath+'/tmp/'
    if not os.path.exists(tmppath):
      os.makedirs(tmppath)
    pname = tmppath+'phimax.ts='+str(ts)+'-te='+str(te)+'.p'
    pickle.dump([tt,xx,zz,phimax], open(pname, "wb" ),2)
    print "phimax was stored for "+datapath
    
 
  return tt,xx,zz,phimax

def plot_Jy(targ,probid='cs',iframe=0,jmin=None,jmax=None,vmin=-1,vmax=1,aspect=1,cmap='jet',zmin=None,zmax=None,filename=None):
  """ load pickled data and generated Jy and vx plot"""
  datapath=get_path(targ)
  pad,fraction = 0.01,0.02;
  tt,xx,zz,bx,bz,vx = loadpickle(datapath,iframe)
  dx = xx[1]-xx[0];dz = zz[1]-zz[0]
  extent = [np.min(xx),np.max(xx),np.min(zz),np.max(zz)]
  j = np.gradient(bx,axis=0)/dz -np.gradient(bz,axis=1)/dx
  
  fig = plt.figure()
  plt.subplot(1,2,1) 
  plt.imshow((j),cmap=cmap,origin='lower',extent=extent,aspect=aspect,vmin=jmin,vmax=jmax)
  plt.title('t='+str("{:4.1f}".format(tt))+', '+'$J_y$=['+str("{:5.2e}".format(np.min(j)))+ ', '+\
                str("{:5.2e}".format(np.max(j)))+']')
  plt.colorbar(pad=pad,fraction=fraction)
  if zmin and zmax:
    plt.ylim([zmin,zmax])
  plt.subplot(1,2,2) 
  plt.imshow((vx),cmap=cmap,origin='lower',extent=extent,aspect=aspect,vmin=vmin,vmax=vmax)
  plt.title('$v_x$=['+str("{:5.2e}".format(np.min(vx)))+ ', '+\
                str("{:5.2e}".format(np.max(vx)))+']')
  plt.colorbar(pad=pad,fraction=fraction)
  if zmin and zmax:
    plt.ylim([zmin,zmax])
  #plt.show()
  if filename:
    fig.savefig(filename, bbox_inches='tight')


if __name__=='__main__':
  if len(sys.argv)>=6:
    ts=int(sys.argv[2])
    te=int(sys.argv[3])
    lev=int(sys.argv[4])
    zrange=float(sys.argv[5])
  else:
    ts=0;te=100;lev=4;zrange=0.05
    print "use ts=0 te=100 for the calc" 
    print "plz specify ts/te if other than that" 
  targ = sys.argv[1]
  datapath = get_path(targ)
  for iframe in np.arange(ts,te+1):
    time,xx,zz,bx,bz,vx = loaddata(datapath,iframe,level=lev,zrange=zrange)
    pickledata(time,xx,zz,bx,bz,vx,datapath,iframe)


