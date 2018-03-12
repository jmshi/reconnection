import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import athena4_read as ath
import athena_read as ath
import cPickle as pickle
import sys


def loadData(fname='Unstra.out2.00008.athdf'):
  """load 3d bfield and calc the current density"""
  #data=ath.athdf(fname,quantities=['B1','B2','B3'])
  time,data=ath.athdf(fname,quantities=['vel1'])
  vx = data['vel1']
  time,data=ath.athdf(fname,quantities=['vel2'])
  vy = data['vel2']
  time,data=ath.athdf(fname,quantities=['vel3'])
  vz = data['vel3']
  x  = data['x1f']
  y  = data['x2f']
  z  = data['x3f']
  # ---
  def curl(vx,vy,vz,dx,dy,dz):
    [dzvx,dyvx,dxvx] = np.gradient(vx)
    [dzvy,dyvy,dxvy] = np.gradient(vy)
    [dzvz,dyvz,dxvz] = np.gradient(vz)
    cx = dyvz/dy-dzvy/dz
    cy = dzvx/dz-dxvz/dx
    cz = dxvy/dx-dyvx/dy
    # No need to del the reference by one manually
    # allow python to perform its own garbage collection
    # after the function return cxyz
    #del dzvx
    #del dzvy
    #del dzvz
    return cx,cy,cz
  # ---
  dx = dz = x[1]-x[0]
  dy = y[1]-y[0]
  jx,jy,jz = curl(vx,vy,vz,dx,dy,dz)
  w2 = jx**2+jy**2+jz**2
  del jx,jy,jz,vx,vy,vz
  return w2


def calc_cdf(j2,fname,fsize=100):
  """
  calc and store the cumulative distribution of j2
  """
  # normalized with jrms^2
  jrms2 = np.average(j2)
  j2 /= jrms2
  j2tot = np.sum(j2)
  jaxis = np.linspace(0,20,num=fsize,endpoint=True)
  jaxis = jaxis**2 # (j/jrms)^2
  voltot = np.count_nonzero(j2+1e-3)
  # dump the cdf data as table
  fhandler=open(fname+'.vcdf.tab','a')
  for i in np.arange(fsize):
    jval = jaxis[i]
    cdf  = np.sum(j2[j2>jval])/j2tot
    fvol  = np.count_nonzero(j2[j2>jval])/float(voltot)
    tmp = np.array([jval,cdf,fvol]).reshape(1,3)
    np.savetxt(fhandler,tmp,fmt='%.10e')
  fhandler.close()



  
if __name__=='__main__':
    if len( sys.argv ) < 4:
        print "Please specify input targ,ts,te,tstride" 
        exit(  )
    targ = sys.argv[1]
    #targ = x2y4z1r64pm1re4000
    ts,te,tstride = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
    fdir = '/tigress/jiming/reconnect/athena/bin/'
    basename = 'Unstra.out2.'
    for frame in np.arange(ts,te,tstride):
      fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.athdf'	    
      w2 = loadData(fname)
      calc_cdf(w2,fname)
 


# end of the script

