#import pyfftw
import my_athena_read as ath
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import fftpack
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
import sys
import os
import pandas as pd

def shear_map(x,dy,qomt,datai,flag=1):
  """ 
  depends on the flag, perform the forward (flag=1) 
  or backward (flag=-1) transfer between shearing
  periodic and exact periodic
  """
  if flag !=1 and flag !=-1:
    print "shear_map: incorrect flag,+1 or -1"
    return 0
  else:
    ndim = datai.ndim
    dim  = np.array(datai.shape)# datai[nz,ny,nx]
    sh_data = np.empty([dim[0],dim[1],dim[2]],dtype='float64')
    tp_data = np.empty([dim[0],dim[1]],dtype='float64')
    sh_y = -qomt*x/dy #-qshear*omega*dtn*x/dy
    
    for i in np.arange(0,dim[2]):
      quot = int(np.floor(sh_y[i]))
      res  = sh_y[i]-np.float(quot)
      tp_data[:,:] = datai[:,:,i]
      sh_data[:,:,i] = (1.0-res)*np.roll(tp_data,flag*quot,axis=1)\
                +res*np.roll(tp_data,flag*(quot+1),axis=1)
    
    #print type(sh_data)
    #print sh_data.shape
    return sh_data


def remap(kx,ky,lx,ly,qomt,datai):
  """
  remap the k-space variable back to shearing 
  periodic frame to reflect the time dependent 
  Eulerian wave number
  """
  ndim = datai.ndim
  dim  = np.array(datai.shape)# datai[nz,ny,nx]
  sh_data = np.empty([dim[0],dim[1],dim[2]])
  tp_data = np.empty([dim[0],dim[2]])
  sh_kx = -qomt*ky*lx/ly
  #nquist= np.max(np.fabs(kx))

  for j in np.arange(0,dim[1]):
    quot = int(np.floor(sh_kx[j]))
    res  = sh_kx[j]-float(quot)
    #kx_new = kx[:] + sh_kx[j]
    tp_data[:,:]= datai[:,j,:]
    sh_data[:,j,:] = (1.0-res)*np.roll(tp_data,quot, axis=1) \
                         + res*np.roll(tp_data,quot+1,axis=1)
    #sh_data[:,j,kx_new[:]>nquist] = 0.0

  return sh_data

def onedim_int(kx,ky,kz,lx,ly,lz,datai):
  """
  average the spectra over constant k_i
  where i could be x,y,or z, the direction
  of the guide field,i.e. k_parallel
  """
  nx = kx.shape[0]; ny = ky.shape[0]; nz = kz.shape[0]
  # k_i = kx
  wcnt   = int(nx/2)+1
  kmodx  = np.arange(wcnt) 
  powerx = np.zeros(wcnt)
  for i in xrange(0,wcnt):
    powerx[i] = np.sum(datai[:,:,i])
    if (np.abs(i) != 1) and (np.abs(i) != nx/2):
      powerx[i] += np.sum(datai[:,:,-i])
  # k_i = ky
  wcnt   = int(ny/2)+1
  kmody  = np.arange(wcnt) 
  powery = np.zeros(wcnt)
  for i in xrange(0,wcnt):
    powery[i] = np.sum(datai[:,i,:])
    if (np.abs(i) != 1) and (np.abs(i) != ny/2):
      powery[i] += np.sum(datai[:,-i,:]) 
  # k_i = kz
  wcnt   = int(nz/2)+1
  kmodz  = np.arange(wcnt) 
  powerz = np.zeros(wcnt)
  for i in xrange(0,wcnt):
    powerz[i] = np.sum(datai[i,:,:])
    if (np.abs(i) != 1) and (np.abs(i) != nz/2):
      powerz[i] += np.sum(datai[-i,:,:]) 

  return kmodx,kmody,kmodz,powerx,powery,powerz


def shell_int(kx,ky,kz,lx,ly,lz,datai):
  """
  average the spectra over spherical shells of 
  constant k
  """
  nx = kx.shape[0]; ny = ky.shape[0]; nz = kz.shape[0]
  wcnt = int((min([nx,ny,nz]))/2)+1
  kmod = np.arange(wcnt) 
  power = np.zeros(wcnt)
  k3,k2,k1 = np.meshgrid(kz,ky,kx,indexing='ij')
  kmod3d = (np.sqrt((k1/lx)**2+(k2/ly)**2+(k3/lz)**2)+0.5).astype(int)
  for i in xrange(0,wcnt):
    #power[i] = np.sum(datai[np.where(kmod3d == kmod[i])])
    power[i] = np.average(datai[np.where(kmod3d == kmod[i])])
  return kmod,power


def powerspectra(datai,x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp,int_opt=None):
  """ 
  calling procedure to calculate the power
  spectrum of a given datai in a shearing
  box, and return datao
  """
  # using scipy.fftpack: slow when array size is big
  fft_unwrapped = fftpack.fftn(shear_map(x,dy,qomt,datai*mtp))
  fft_unwrapped = np.real(fft_unwrapped*np.conj(fft_unwrapped))
  fft_unwrapped = remap(kx,ky,lx,ly,qomt,fft_unwrapped)   
  if int_opt == 'shell':
    kmod,power    = shell_int(kx,ky,kz,lx,ly,lz,fft_unwrapped)
    return kmod,power
  elif int_opt == 'onedim':
    kmodx,kmody,kmodz,powerx,powery,powerz = onedim_int(kx,ky,kz,lx,ly,lz,fft_unwrapped)
    return kmodx,kmody,kmodz,powerx,powery,powerz
  else:
    kmod,power    = shell_int(kx,ky,kz,lx,ly,lz,fft_unwrapped)
    kmodx,kmody,kmodz,powerx,powery,powerz = onedim_int(kx,ky,kz,lx,ly,lz,fft_unwrapped)
    return kmod,kmodx,kmody,kmodz,power,powerx,powery,powerz
  
  ### using pyfftw wrapper of FFTW3
  ##nx,ny,nz = len(kx),len(ky),len(kz)
  ##fft_unwrapped = pyfftw.empty_aligned((nz,ny,nx), dtype='complex128')
  ##datai = shear_map(x,dy,qomt,datai*mtp)
  ##fft_unwrapped[:] = datai
  ##fft_object = pyfftw.builders.fft(fft_unwrapped)
  ##power = fft_object()
  ##power = remap(kx,ky,lx,ly,qomt,np.real(power*np.conj(power)))
  ##kmod,power = shell_int(kx,ky,kz,lx,ly,lz,power)


def plot_pspec1d(targname,ts=50,te=100,stride=10):
  """
  plot the energy spectra (B and v) 
  """
  dirname = '/tigress/jiming/reconnect/athena/bin/'
  #targname = 'x2y4z1r128pm0.5re3000' #x2y8z1r128pm0.5re3000/'
  ncells = 256*256*4*128# default 2x8x1 with 128/H
  if targname[0:6] == 'x2y4z1':
    ncells /=2
  if targname[6:9] == 'r64':
    ncells /=8
  else:
    resol = int(targname[7:10])
    ncells = ncells * (resol/128)**3
  if targname[10:12] == 'ry':
    yresol = int(targname[12:14])
    ncells = int(ncells*yresol/128)
  
  fnorm = float(ncells)**2*2.0 # extra factor of 1/2 for B^2/2 and rho v^2/2
  bhist = 0; ahist = 0; cnt = 0
  for i in np.arange(ts,te,stride):
    bhist += np.loadtxt(dirname+targname+'/'+'Unstra-1b.'+str(i).zfill(5)+'.pwr', skiprows=1)
    ahist += np.loadtxt(dirname+targname+'/'+'Unstra-1v.'+str(i).zfill(5)+'.pwr', skiprows=1)
    cnt += 1 
  ahist /= float(cnt)
  bhist /= float(cnt)
  #matplotlib.rcParams['figure.figsize'] = (10,6)
  plt.plot(ahist[:,0],(ahist[:,1]+bhist[:,1])*ahist[:,0]**2/fnorm,'.-',lw=2,label=r'$B^2/2+\rho v^2/2$')
  plt.plot(ahist[:,0],ahist[:,1]*ahist[:,0]**2/fnorm,'.-',lw=2,label=r'$\rho v^2/2$')
  #plt.plot(ahist[:,0],ahist[:,2]*ahist[:,0]**2/fnorm,'.-',lw=2,label=r'$\rho v_x^2$')
  #plt.plot(ahist[:,0],ahist[:,3]*ahist[:,0]**2/fnorm,'.-',lw=2,label=r'$\rho v_y^2$')
  #plt.plot(ahist[:,0],ahist[:,4]*ahist[:,0]**2/fnorm,'.-',lw=2,label=r'$\rho v_z^2$')
  plt.plot(ahist[:,0],bhist[:,1]*bhist[:,0]**2/fnorm,'.-',lw=2,label=r'$B^2/2$')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim([1,100])
  plt.ylim([1e-6,1])
  amp = 2
  plt.plot(np.arange(0,100,0.1),amp*np.arange(0.01,100,0.1)**(-1.5),':',lw=2)
  #plt.plot(np.arange(0,100,0.1),amp*np.arange(0,100,0.1)**(-1),':',lw=2)
  #plt.plot(np.arange(0,100,0.1),amp*np.arange(0,100,0.1)**(-2),':',lw=2)
  plt.legend(fontsize=20,loc=3)



def get_pspec1dall(targname,nkx=128,nky=126,nkz=64,ts=50,te=100,stride=10,noby=False):
  """
  return the energy spectra (B and v) for plots
  """
  dirname = '/tigress/jiming/reconnect/athena/bin/'
  #targname = 'x2y4z1r128pm0.5re3000' #x2y8z1r128pm0.5re3000/'
  ncells = 256*256*4*128# default 2x8x1 with 128/H
  if targname[0:6] == 'x2y4z1':
    ncells /=2
  if targname[6:9] == 'r64':
    ncells /=8
  else:
    resol = int(targname[7:10])
    ncells = ncells * (resol/128)**3
  if targname[10:12] == 'ry':
    yresol = int(targname[12:14])
    ncells = int(ncells*yresol/128)
  
  fnorm = float(ncells)**2*2.0 # extra factor of 1/2 for B^2/2 and rho v^2/2
  bhist = 0; ahist = 0; cnt = 0
  bkxhist = 0; akxhist = 0
  bkyhist = 0; akyhist = 0
  bkzhist = 0; akzhist = 0
  nk = np.min([nkx,nky,nkz])
  for i in np.arange(ts,te,stride):
    fname = dirname+targname+'/'+'Unstra-spec1b.'+str(i).zfill(5)+'.pwr'
    df = pd.read_table(fname,delimiter=' ',skiprows=1,nrows=nk+1,header=None)
    bhist += df.values
    df = pd.read_table(fname,delimiter=' ',skiprows=1+nk+2,nrows=nkx+1,header=None)
    bkxhist += df.values
    df = pd.read_table(fname,delimiter=' ',skiprows=1+nk+nkx+4,nrows=nky+1,header=None)
    bkyhist += df.values
    df = pd.read_table(fname,delimiter=' ',skiprows=1+nk+nkx+nky+6,nrows=nkz+1,header=None)
    bkzhist += df.values
    #bhist += np.loadtxt(dirname+targname+'/'+'Unstra-1b.'+str(i).zfill(5)+'.pwr', skiprows=1)
    fname = dirname+targname+'/'+'Unstra-spec1v.'+str(i).zfill(5)+'.pwr'
    df = pd.read_table(fname,delimiter=' ',skiprows=1,nrows=nk+1,header=None)
    ahist += df.values
    df = pd.read_table(fname,delimiter=' ',skiprows=1+nk+2,nrows=nkx+1,header=None)
    akxhist += df.values
    df = pd.read_table(fname,delimiter=' ',skiprows=1+nk+nkx+4,nrows=nky+1,header=None)
    akyhist += df.values
    df = pd.read_table(fname,delimiter=' ',skiprows=1+nk+nkx+nky+6,nrows=nkz+1,header=None)
    akzhist += df.values
    #ahist += np.loadtxt(dirname+targname+'/'+'Unstra-1v.'+str(i).zfill(5)+'.pwr', skiprows=1)
    cnt += 1 

  for i in np.arange(4):
    if i == 0:
      fn = float(cnt)
    else:
      fn = fnorm*float(cnt)
    ahist[:,i] /= fn;akxhist[:,i] /= fn;akyhist[:,i] /= fn; akzhist[:,i] /=fn 
    bhist[:,i] /= fn;bkxhist[:,i] /= fn;bkyhist[:,i] /= fn; bkzhist[:,i] /=fn 

  kmod,kmodx,kmody,kmodz,pwra,pwrb,pwrax,pwrbx,pwray,pwrby,pwraz,pwrbz = \
  ahist[:,0],akxhist[:,0],akyhist[:,0],akzhist[:,0],\
  ahist[:,1]*ahist[:,0]**2,bhist[:,1]*bhist[:,0]**2,\
  akxhist[:,1]*akxhist[:,0],bkxhist[:,1]*bkxhist[:,0],\
  akyhist[:,1]*akyhist[:,0],bkyhist[:,1]*bkyhist[:,0],\
  akzhist[:,1]*akzhist[:,0],bkzhist[:,1]*bkzhist[:,0]


  return kmod,kmodx,kmody,kmodz,pwra,pwrb,pwrax,pwrbx,pwray,pwrby,pwraz,pwrbz


def get_pspec1d(targname,ts=50,te=100,stride=10,nk=256,noby=False):
  """
  return the energy spectra (B and v) for plots
  """
  dirname = '/tigress/jiming/reconnect/athena/bin/'
  #targname = 'x2y4z1r128pm0.5re3000' #x2y8z1r128pm0.5re3000/'
  ncells = 256*256*4*128# default 2x8x1 with 128/H
  if targname[0:6] == 'x2y4z1':
    ncells /=2
  if targname[6:9] == 'r64':
    ncells /=8
  else:
    resol = int(targname[7:10])
    ncells = ncells * (resol/128)**3
  if targname[10:12] == 'ry':
    yresol = int(targname[12:14])
    ncells = int(ncells*yresol/128)
  
  fnorm = float(ncells)**2*2.0 # extra factor of 1/2 for B^2/2 and rho v^2/2
  bhist = 0; ahist = 0; cnt = 0
  for i in np.arange(ts,te,stride):
    df = pd.read_table(dirname+targname+'/'+'Unstra-spec1b.'+str(i).zfill(5)+'.pwr',delimiter=' ',skiprows=1,nrows=nk+1,header=None)
    bhist += df.values
    df = pd.read_table(dirname+targname+'/'+'Unstra-spec1v.'+str(i).zfill(5)+'.pwr',delimiter=' ',skiprows=1,nrows=nk+1,header=None)
    ahist += df.values
    cnt += 1 
  ahist /= float(cnt)
  bhist /= float(cnt)
  if noby :
    kmod,pwra,pwrb,pwrby = ahist[:,0],ahist[:,1]*ahist[:,0]**2/fnorm,(bhist[:,2]+bhist[:,4])*bhist[:,0]**2/fnorm, bhist[:,3]*bhist[:,0]**2/fnorm
    return kmod,pwra,pwrb,pwrby
  else:
    kmod,pwra,pwrb = ahist[:,0],ahist[:,1]*ahist[:,0]**2/fnorm,bhist[:,1]*bhist[:,0]**2/fnorm
    return kmod,pwra,pwrb


if __name__ == '__main__':
  """
  calling sequence:
  python pspec.py targname [ts] [te] [stride]
  """
  dirname = '/tigress/jiming/reconnect/athena/bin/'
  #dirname = '/tigress/jiming/reconnect/athena.idealMHD/bin/'
  targname = sys.argv[1] #x2y8z1r128pm0.5re3000/'
  ts,te,tstride = 50,100,10
  qshear, omg = 1.5, 1.0
  lx,ly,lz = 2.0,8.0,1.0
  nx,ny,nz = 256,1008,128

  if targname[0:6] == 'x2y4z1':
    ly,ny = 4.0,504
  if targname[0:6] == 'x4y4z1':
    lx,ly,nx,ny = 4.0,4.0,512,504
  if targname[10:14] == 'ry64':
    ny = int(ny/2)

  if len(sys.argv) > 2: # plz specify ts te and tstride
    ts = int(sys.argv[2])
    te = int(sys.argv[3])
    tstride = int(sys.argv[4])

  if len(sys.argv)> 5:
    lx,ly,lz = int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7])
    nx,ny,nz = int(sys.argv[8]),int(sys.argv[9]),int(sys.argv[10])


  for i in np.arange(ts,te,tstride):
    fname = dirname+targname+'/'+'Unstra.out2.'+str(i).zfill(5)+'.athdf'
    time, grid = ath.athdf(fname,quantities=['x1f','x2f','x3f'])
    if (i == ts): # get the grid info
      x = grid['x1f']+0.5*(grid['x1f'][1]-grid['x1f'][0])
      x = x[:-1]
      dx = grid['x1f'][1]-grid['x1f'][0]
      dy = grid['x2f'][1]-grid['x2f'][0]
      dz = grid['x3f'][1]-grid['x3f'][0]
      nx = len(grid['x1f'])-1
      ny = len(grid['x2f'])-1
      nz = len(grid['x3f'])-1
      lx = grid['x1f'][-1]-grid['x1f'][0]
      ly = grid['x2f'][-1]-grid['x2f'][0]
      lz = grid['x3f'][-1]-grid['x3f'][0]
      kx = np.roll(np.arange(nx)-nx/2+1, nx/2+1)
      ky = np.roll(np.arange(ny)-ny/2+1, ny/2+1)
      kz = np.roll(np.arange(nz)-nz/2+1, nz/2+1)
    # calculate the shear amount
    nt = np.rint(time*qshear*omg*lx/ly)
    dtn = time - np.float64(nt)*ly/(qshear*omg*lx)
    qomt = qshear*omg*dtn

    print 'analyze '+fname+': '
    mtp = 1.0 # weight
    time, data = ath.athdf(fname,quantities=['Bcc1'])
    print 'calc pwrspec of Bcc1'
    #kmod,pwr1 = powerspectra(data['Bcc1'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    kmod,kmodx,kmody,kmodz,pwr1,pwr1x,pwr1y,pwr1z = powerspectra(data['Bcc1'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    time, data = ath.athdf(fname,quantities=['Bcc2'])
    print 'calc pwrspec of Bcc2'
    #kmod,pwr2 = powerspectra(data['Bcc2'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    kmod,kmodx,kmody,kmodz,pwr2,pwr2x,pwr2y,pwr2z = powerspectra(data['Bcc2'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    time, data = ath.athdf(fname,quantities=['Bcc3'])
    print 'calc pwrspec of Bcc3'
    #kmod,pwr3 = powerspectra(data['Bcc3'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    kmod,kmodx,kmody,kmodz,pwr3,pwr3x,pwr3y,pwr3z = powerspectra(data['Bcc3'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    ## the convention here is: 1,2,3 for B-component; x,y,z for k-axis
    pwr  = pwr1 + pwr2 + pwr3
    pwrx  = pwr1x + pwr2x + pwr3x 
    pwry  = pwr1y + pwr2y + pwr3y
    pwrz  = pwr1z + pwr2z + pwr3z

    dumpname = dirname+targname+'/'+'Unstra-spec1b.'+str(i).zfill(5)+'.pwr'
    dumpdata = Table([kmod,pwr,pwr1,pwr2,pwr3],names=['kmod','pwr','pwr1','pwr2','pwr3'])
    dumpdata.write(dumpname,format='ascii')

    dumpdata = Table([kmodx,pwrx,pwr1x,pwr2x,pwr3x],names=['kmodx','pwrx','pwr1x','pwr2x','pwr3x'])
    with open(dumpname,mode='a') as f:
      f.seek(0,os.SEEK_END)
      dumpdata.write(f,format='ascii')
      dumpdata = Table([kmody,pwry,pwr1y,pwr2y,pwr3y],names=['kmody','pwry','pwr1y','pwr2y','pwr3y'])
      f.seek(0,os.SEEK_END)
      dumpdata.write(f,format='ascii')
      dumpdata = Table([kmodz,pwrz,pwr1z,pwr2z,pwr3z],names=['kmodz','pwrz','pwr1z','pwr2z','pwr3z'])
      f.seek(0,os.SEEK_END)
      dumpdata.write(f,format='ascii')

    print 'dumped pwrspec to '+dumpname

    time, data = ath.athdf(fname,quantities=['rho'])
    mtp = np.sqrt(data['rho'])
    vshear = qshear*omg*np.resize(x,(nz,ny,nx))
    time, data = ath.athdf(fname,quantities=['vel1'])
    print 'calc pwrspec of vel1'
    #kmod,pwr1 = powerspectra(data['vel1'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    kmod,kmodx,kmody,kmodz,pwr1,pwr1x,pwr1y,pwr1z = powerspectra(data['vel1'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    time, data = ath.athdf(fname,quantities=['vel2'])
    print 'calc pwrspec of vel2'
    #kmod,pwr2 = powerspectra(data['vel2']+vshear,x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    kmod,kmodx,kmody,kmodz,pwr2,pwr2x,pwr2y,pwr2z = powerspectra(data['vel2']+vshear,x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    time, data = ath.athdf(fname,quantities=['vel3'])
    print 'calc pwrspec of vel3'
    #kmod,pwr3 = powerspectra(data['vel3'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    kmod,kmodx,kmody,kmodz,pwr3,pwr3x,pwr3y,pwr3z = powerspectra(data['vel3'],x,dy,kx,ky,kz,lx,ly,lz,qomt,mtp)
    pwr  = pwr1 + pwr2 + pwr3
    pwrx  = pwr1x + pwr2x + pwr3x 
    pwry  = pwr1y + pwr2y + pwr3y
    pwrz  = pwr1z + pwr2z + pwr3z

    dumpname = dirname+targname+'/'+'Unstra-spec1v.'+str(i).zfill(5)+'.pwr'
    dumpdata = Table([kmod,pwr,pwr1,pwr2,pwr3],names=['kmod','pwr','pwr1','pwr2','pwr3'])
    dumpdata.write(dumpname,format='ascii')
    dumpdata = Table([kmodx,pwrx,pwr1x,pwr2x,pwr3x],names=['kmodx','pwrx','pwr1x','pwr2x','pwr3x'])
    with open(dumpname,mode='a') as f:
      f.seek(0,os.SEEK_END)
      dumpdata.write(f,format='ascii')
      dumpdata = Table([kmody,pwry,pwr1y,pwr2y,pwr3y],names=['kmody','pwry','pwr1y','pwr2y','pwr3y'])
      f.seek(0,os.SEEK_END)
      dumpdata.write(f,format='ascii')
      dumpdata = Table([kmodz,pwrz,pwr1z,pwr2z,pwr3z],names=['kmodz','pwrz','pwr1z','pwr2z','pwr3z'])
      f.seek(0,os.SEEK_END)
      dumpdata.write(f,format='ascii')

    print 'dumped pwrspec to '+dumpname



