import numpy as np
import my_athena_read as myath
import os
import sys

def dump_hist2d(targname,ts=0,te=5,resol=64,pres=1.0,nx=256,nz=64,nvar=3):
  direname='/tigress/jiming/reconnect/athena/bin/'
  basename='Unstra.out2.'
  appdname='athdf'
  dumprate = 0.5  # in unit of orbits
  nframe = te-ts+1
  have = np.zeros((nz,nframe,nvar))
  zave = np.zeros((nx,nframe,nvar))
  tave = np.zeros(nframe)
  
  for i in range(ts,te+1):
    fname = direname+targname+'/'+basename+str(i).zfill(5)+'.'+appdname
    if os.path.isfile(fname): 
      time,data = myath.athdf(fname,quantities=['Bcc1','Bcc2','Bcc3','rho'])
      tave[i]=time
      have[:,i,0] = np.average(data['rho'],axis=(1,2)) # rho
      have[:,i,1] = np.average(data['Bcc2'],axis=(1,2)) # By
      have[:,i,2] = np.average(-data['Bcc2']*data['Bcc1'],axis=(1,2))/pres # maxwell stress
      zave[:,i,0] = np.average(data['rho'],axis=(0,1)) # rho
      zave[:,i,1] = np.average(data['Bcc2'],axis=(0,1)) # By
      zave[:,i,2] = np.average(-data['Bcc2']*data['Bcc1'],axis=(0,1))/pres # maxwell stress

  # dump have/zave data
  fmt = '%.15e'
  halfcellwidth=(data['x1f'][-1:]-data['x1v'])[0]
  xx = data['x1f'][:-1]+halfcellwidth
  zz = data['x3f'][:-1]+halfcellwidth
  dst = direname+targname+'/'+'tmp/'
  if not os.path.exists(dst):
    os.makedirs(dst)
  np.savetxt(dst+'xx.dat',xx,fmt=fmt)
  np.savetxt(dst+'zz.dat',zz,fmt=fmt)
  np.savetxt(dst+'tt.dat',tave,fmt=fmt)
  for i in range(nvar):
    np.savetxt(dst+'have'+str(i).zfill(2)+'.dat',have[...,i],fmt=fmt)
    np.savetxt(dst+'zave'+str(i).zfill(2)+'.dat',zave[...,i],fmt=fmt)
  return


if __name__=='__main__':
    if len(sys.argv)>=4:
	ts=int(sys.argv[2])
	te=int(sys.argv[3])
    else:
	ts=0;te=100
        print "use ts=0 te=100 for the calc" 
        print "plz specify ts/te if other than that" 

    resol=64;pres=1.0;nx=4*64;nz=64*1;nvar=3
    targ = sys.argv[ 1 ];
    if(targ[0:6]=='x2y4z1'):
	    Lx,Lz=2,1
    if(targ[0:6]=='x4y4z1'):
	    Lx,Lz=4,1
    if(targ[0:6]=='x2y8z1'):
	    Lx,Lz=2,1
    if(targ[6:9]=='r64'):
	    resol=64
    if(targ[6:10]=='r128'):
	    resol=128
    if(targ[6:10]=='r256'):
	    resol=256
    if(targ[6:10]=='r512'):
	    resol=512
    if(targ[0:8]=='adb.r256'):
	    resol=256
	    Lx,Ly,Lz = 2,4,1
    nx=resol*Lx
    nz=resol*Lz
    
    print 'calculating have/zave for '+targ+' within ts= '+str(ts)+' te= '+str(te)
    dump_hist2d(targ,ts=ts,te=te,resol=resol,nx=nx,nz=nz,nvar=nvar)








