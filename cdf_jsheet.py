"""
 driver of the whole pipe line
 for finding current sheet and null pts
"""
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
  time,data=ath.athdf(fname,quantities=['Bcc1'])
  bx = data['Bcc1']
  time,data=ath.athdf(fname,quantities=['Bcc2'])
  by = data['Bcc2']
  time,data=ath.athdf(fname,quantities=['Bcc3'])
  bz = data['Bcc3']
  x  = data['x1f']
  y  = data['x2f']
  z  = data['x3f']
  # refinement
  rfac = 1.0
  ##if bx.shape[0] < 512:
  ##  nz,ny,nx = bx.shape
  ##  rfac = int(512/bx.shape[0])
  ##  bx = np.repeat(bx,rfac,axis=0)
  ##  bx = np.repeat(bx,rfac,axis=1)
  ##  bx = np.repeat(bx,rfac,axis=2)
  ##  by = np.repeat(by,rfac,axis=0)
  ##  by = np.repeat(by,rfac,axis=1)
  ##  by = np.repeat(by,rfac,axis=2)
  ##  bz = np.repeat(bz,rfac,axis=0)
  ##  bz = np.repeat(bz,rfac,axis=1)
  ##  bz = np.repeat(bz,rfac,axis=2)
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
  dx = dz = (x[1]-x[0])/rfac
  dy = (y[1]-y[0])/rfac
  jx,jy,jz = curl(bx,by,bz,dx,dy,dz)
  j2 = jx**2+jy**2+jz**2
  return j2


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
  fhandler=open(fname+'.jcdf.tab','a')
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
      j2 = loadData(fname)
      calc_cdf(j2,fname)
 


# end of the script

