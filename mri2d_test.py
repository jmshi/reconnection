import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import athena_read as ath
import sys
from scipy import integrate

def read_all_hdf5(fname,level=0,subsample=False):
  time,data=ath.athdf(fname,level=level,subsample=subsample)
  bx = data['Bcc1'][0]
#time,data=ath.athdf(fname,quantities=['Bcc2'])
  by = data['Bcc2'][0]
#time,data=ath.athdf(fname,quantities=['Bcc3'])
  bz = data['Bcc3'][0]
#time,data=ath.athdf(fname,quantities=['vel1'])
  vx = data['vel1'][0]
#time,data=ath.athdf(fname,quantities=['vel2'])
  vy = data['vel2'][0]
#time,data=ath.athdf(fname,quantities=['vel3'])
  vz = data['vel3'][0]
  rho = data['rho'][0]
  T = data['press'][0]/rho
  x  = data['x1f'];y  = data['x2f']; z  = data['x3f']
  x = (x + 0.5*(x[1]-x[0]))[:-1]
  y = (y + 0.5*(y[1]-y[0]))[:-1]
  z = (z + 0.5*(z[1]-z[0]))[:-1]
  return time,x,y,z,rho,T,bx,by,bz,vx,vy,vz

def read_all_vtk(fname):
  time,x,y,z,data=ath.vtk(fname)
  bx = data['Bcc'][0,...,0]
#time,data=ath.athdf(fname,quantities=['Bcc2'])
  by = data['Bcc'][0,...,1]
#time,data=ath.athdf(fname,quantities=['Bcc3'])
  bz = data['Bcc'][0,...,2]
#time,data=ath.athdf(fname,quantities=['vel1'])
  vx = data['vel'][0,...,0]
#time,data=ath.athdf(fname,quantities=['vel2'])
  vy = data['vel'][0,...,1]
#time,data=ath.athdf(fname,quantities=['vel3'])
  vz = data['vel'][0,...,2]
  rho = data['rho'][0]
  T = data['press'][0]/rho
  x = (x + 0.5*(x[1]-x[0]))[:-1]
  y = (y + 0.5*(y[1]-y[0]))[:-1]
  if len(z) > 1:
    z = (z + 0.5*(z[1]-z[0]))[:-1]
  else:
    z = z
  return time,x,y,z,rho,T,bx,by,bz,vx,vy,vz

def read_bz_hdf5(fname,level=0,subsample=False):
  time,data=ath.athdf(fname,level=level,subsample=subsample,quantities=['Bcc3'])
  bz = data['Bcc3'][0]
  x  = data['x1f'];y  = data['x2f']; z  = data['x3f']
  x = (x + 0.5*(x[1]-x[0]))[:-1]
  y = (y + 0.5*(y[1]-y[0]))[:-1]
  z = (z + 0.5*(z[1]-z[0]))[:-1]
  return time,x,y,z,bz

def plot_all_hdf5(time,x,y,z,rho,T,bx,by,bz,vx,vy,vz,cmap='RdBu_r',to_png=False,pname=None):
    matplotlib.rcParams['figure.figsize'] = (6,12)
    nrow,ncol,count = 4,2,8

    fig,axes = plt.subplots(nrows=nrow,ncols=ncol)#,sharey=True,sharex=True)
    
    for i in np.arange(nrow*ncol):
      if (i==0):
        label=r'$v_r$'
        dmin,dmax = -0.14,0.14
        var = vx
      if (i==1):
        label=r'$B_r$'
        var = bx
        dmin,dmax = -0.18,0.18
        dmin,dmax = -0.3,0.3
      if (i==2):
        label=r'$v_z$'
        var = vy
        dmin,dmax = -0.12,0.12
      if (i==3):
        label=r'$B_z$'
        var = by
        dmin,dmax = -0.1,0.1
        dmin,dmax = -0.3,0.3
      if (i==4):
        label=r'$v_{\phi}$'
        var = vz
        dmin,dmax = -0.4,0.4
        dmin,dmax = -0.6,0.6
      if (i==5):
        label=r'$B_{\phi}$'
        var = bz
        dmin,dmax = -0.22,0.22
        dmin,dmax = -0.6,0.6
      if (i==6):
        label=r'$log \rho$'
        var = np.log10(rho)
        dmin,dmax = -0.012,0.012
      if (i==7):
        label=r'$\delta T$'
        var = T-np.average(T) #np.log10(T/10.0)
        dmin,dmax = -0.05,0.05
        dmin,dmax = -0.1,0.1

      #plt.subplot(nrow,ncol,i+1)    
      #im = plt.imshow(np.transpose(var),vmin=dmin,vmax=dmax,extent=[0,1,-0.5,0.5],cmap='RdBu_r',origin='lower')
#       plt.xlabel('z')
#       plt.ylabel('r')
#       plt.text(1-0.15,0.5-0.1,label,fontsize=15)
#       cb = plt.colorbar(orientation='horizontal',pad=0.02)
#       tick_locator = ticker.MaxNLocator(nbins=5)
#       cb.locator = tick_locator
#       cb.update_ticks()
#       if i==0:
#         plt.text(1.05,0.6,'t='+str(time),fontsize=10)

      idn,jdn = i/ncol,i%ncol
      ax = axes[idn,jdn]
      im = ax.imshow(np.transpose(var),vmin=dmin,vmax=dmax,extent=[0,1,-0.5,0.5],cmap=cmap,origin='lower')
      ax.text(1-0.15,0.5-0.1,label,fontsize=15)
      cb = fig.colorbar(im,ax=ax,orientation='horizontal',pad=0.02)
      tick_locator = ticker.MaxNLocator(nbins=5)
      cb.locator = tick_locator
      cb.update_ticks()
      if i==0:
        ax.set_title('t='+str(time),fontsize=10)
        #ax.text(0.5,0.6,'t='+str(time),fontsize=10)
      if i==7:
        ax.set_xlabel('z')
        ax.set_ylabel('r')
    plt.tight_layout()
    # export to png image
    if to_png:
      if pname == None:
        pname = fname+'.png'
      plt.savefig(pname, format='png', dpi=300)
      plt.close(fig)


def plot_one_panel(time,x,y,z,var,label='default',cmap='RdBu_r',to_png=False,pname=None):
    matplotlib.rcParams['figure.figsize'] = (6,12)
    nrow,ncol,count = 1,1,1

    #fig,ax = plt.subplots(nrows=nrow,ncols=ncol)#,sharey=True,sharex=True)
    dmin,dmax = np.min(var),np.max(var)
    print dmin,dmax
    if dmin < 0.0:
      print np.abs(dmin),np.abs(dmax)
      smax = np.maximum(np.abs(dmin),np.abs(dmax))
      dmin,dmax = -smax,smax

    matplotlib.rcParams['figure.figsize'] = (6,6)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(np.transpose(var),vmin=dmin,vmax=dmax,extent=[0,1,-0.5,0.5],cmap=cmap,origin='lower')
    ax.text(1-0.15,0.5-0.1,label,fontsize=15)
    cb = fig.colorbar(im,ax=ax,orientation='horizontal',pad=0.062,fraction=0.045)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    #cb.ax.xaxis.set_ticks_position('top')
    cb.update_ticks()
    ax.set_title('t='+str(time),fontsize=10)
    ax.set_xlabel('z')
    ax.set_ylabel('r')
    # export to png image
    if to_png:
      if pname == None:
        pname = fname+'.png'
      plt.savefig(pname, format='png', dpi=300)
      plt.close(fig)


def plot_aphi_contour(b_z,x,y,nlev=40,t=0,vmin=5e-5,vmax=2e-3,cmap='RdBu_r',to_png=False,pname=None):

  ndim= np.shape(b_z)
  A_phi = np.zeros(ndim)
  for k in np.arange(ndim[0]):
    yint = b_z[k,:]
    xint = x
    phi = integrate.cumtrapz(yint, xint, initial=0)
    A_phi[k,:]=phi[:]

  matplotlib.rcParams['figure.figsize'] = (6,6)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  cs = ax.contour(x,y,A_phi,nlev,vmin=5e-5,vmax=2e-3,cmap=cmap,antialias=True)
  # export to png image
  if to_png:
    if pname == None:
      pname = 'Aphi0000.png'
    plt.savefig(pname, format='png', dpi=300)
    plt.close(fig)


def curl3d(vx,vy,vz,dx,dy,dz):
  [dzvx,dyvx,dxvx] = np.gradient(vx)
  [dzvy,dyvy,dxvy] = np.gradient(vy)
  [dzvz,dyvz,dxvz] = np.gradient(vz)
  cx = dyvz/dy-dzvy/dz
  cy = dzvx/dz-dxvz/dx
  cz = dxvy/dx-dyvx/dy
  return cx,cy,cz

def curl2d(vx,vy,vz,dx,dz):
  [dzvx,dxvx] = np.gradient(vx)
  [dzvy,dxvy] = np.gradient(vy)
  [dzvz,dxvz] = np.gradient(vz)
  cx = -dzvy/dz
  cy = dzvx/dz-dxvz/dx
  cz = dxvy/dx
  return cx,cy,cz

def plot_jcurrent(time,x,y,bx,by,bz,ext1=[-50,50],ext2=[-150,150],ext3=[-80,80],ext4=None,cmap='RdBu_r',to_png=False,pname=None):

  dx,dz = x[1]-x[0],y[1]-y[0]
  jx,jy,jz = curl2d(bx,by,bz,dx,dz)

  matplotlib.rcParams['figure.figsize'] = (6,6)
  nrow,ncol,count = 2,2,4
  fig,axes = plt.subplots(nrows=nrow,ncols=ncol)#,sharey=True,sharex=True)
  for i in np.arange(4):
    if i==0:
      var = jx
      label=r'$j_r$'
      dmin,dmax = ext1[0],ext1[1]
    if i==1:
      var = jy
      label=r'$j_{\phi}$'
      dmin,dmax = ext2[0],ext2[1]
    if i==2:
      var = jz
      label=r'$j_z$'
      dmin,dmax = ext3[0],ext3[1]
    if i==3:
      var = np.sqrt(jx**2+jy**2+jz**2)
      label=r'$|j|$'
      if ext4 == None:
	ext4=[0,0]
	ext4 = np.sqrt(np.array(ext1)**2+np.array(ext2)**2+np.array(ext3)**2)
      dmin,dmax = 0,np.max(ext4)

    idn,jdn = i/ncol,i%ncol
    ax = axes[idn,jdn]
    im = ax.imshow(np.transpose(var),vmin=dmin,vmax=dmax,extent=[0,1,-0.5,0.5],cmap=cmap,origin='lower')
    ax.text(1-0.15,0.5-0.1,label,fontsize=15)
    cb = fig.colorbar(im,ax=ax,orientation='horizontal',pad=0.02)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    if i==0:
      ax.set_title('t='+str(time),fontsize=10)
      #ax.text(0.5,0.6,'t='+str(time),fontsize=10)
    if i==3:
      ax.set_xlabel('z')
      ax.set_ylabel('r')
    #plt.tight_layout()
    plt.show
    # export to png image
    if to_png:
      if pname == None:
        pname = fname+'.png'
      plt.savefig(pname, format='png', dpi=300)
      plt.close(fig)





def batch_png(fdir,ts=0,te=1,tstride=1):

  for i in np.arange(ts,te+1,tstride):
    fname=fdir+'/mri2d.out2.'+str(i).zfill(5)+'.athdf'
    time,x,y,z,rho,T,bx,by,bz,vx,vy,vz = read_all_hdf5(fname)
    plot_all_hdf5(time,x,y,z,rho,T,bx,by,bz,vx,vy,vz,to_png=True,pname=fdir+'/mri2d_t='+str(i).zfill(4)+'.png')
    plot_aphi_contour(by,x,y,nlev=40,t=time,vmin=5e-5,vmax=2e-3,to_png=True,pname=fdir+'/mri2d.aphi.'+str(i).zfill(5)+'.png')
    plot_jcurrent(time,x,y,bx,by,bz,ext1=[-30,30],ext2=[-80,80],ext3=[-40,40],ext4=[0,80],to_png=True,pname=fdir+'/mri2d.jcurr.'+str(i).zfill(5)+'.png')


if __name__=='__main__':
    if len( sys.argv ) < 4:
        print "Please specify input targ,ts,te,tstride" 
        exit(  )
    targ = sys.argv[1]
    ts,te,tstride = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
    fdir = '/tigress/jiming/athena/bin/for_daniel/data/'+targ
    #batch_png(fdir,ts,te,tstride)

    for i in np.arange(ts,te,tstride):
      fname=fdir+'/mri2d.out2.'+str(i).zfill(5)+'.athdf'
      #time,x,y,z,rho,T,bx,by,bz,vx,vy,vz = read_all_hdf5(fname,level=2,subsample=True)
      time,x,y,z,bz = read_bz_hdf5(fname) #,level=2,subsample=True)
      plot_one_panel(time,x,y,z,bz,label=r'$B_{\phi}$',to_png=True,pname=fdir+'/bphi.'+str(i).zfill(5)+'.png')

#ideal 
#1024x1024_nu=1e-6' 
#512x512_nu=1e-5'
