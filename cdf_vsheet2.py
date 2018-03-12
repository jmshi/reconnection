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
  time,data=ath.athdf(fname,quantities=['rho'])
  d = data['rho']
  # ---
  def curl(vx,vy,vz,dx,dy,dz):
    [dzvx,dyvx,dxvx] = np.gradient(vx)
    [dzvy,dyvy,dxvy] = np.gradient(vy)
    [dzvz,dyvz,dxvz] = np.gradient(vz)
    # 1) w^2
    w2 = (dyvz/dy-dzvy/dz)**2
    w2 += (dzvx/dz-dxvz/dx)**2
    w2 += (dxvy/dx-dyvx/dy)**2
    w2 *= d
    # 2) 2*delv:delv
    disp = 2.0*((dxvx/dx)**2 + (dyvy/dy)**2 + (dzvz)**2)
    disp += 4.0*(dxvy*dyvx/dx/dy + dxvz*dzvx/dx/dz + dyvz*dzvy/dy/dz)
    # 3) -2/3 (div v)^2
    disp -=(2.0/3.0)*(dxvx/dx+dyvy/dy+dzvz/dz)**2
    disp *= d
    # No need to del the reference by one manually
    # allow python to perform its own garbage collection
    # after the function return cxyz
    #del dzvx
    #del dzvy
    #del dzvz
    return w2,disp
  # ---
  # 1) get vorticity
  dx = dz = x[1]-x[0]
  dy = y[1]-y[0]
  w2,disp2 = curl(vx,vy,vz,dx,dy,dz)
  #del jx,jy,jz,vx,vy,vz
  # 2) get 
  return w2,disp2


def calc_cdf(w2,k2,fname,fsize=100):
  """
  calc and store the cumulative distribution of j2
  """
  # normalized with jrms^2
  j2 = w2+k2
  jrms2 = np.average(j2)
  j2 /= jrms2
  w2 /= jrms2
  k2 /= jrms2
  j2tot = np.sum(j2)
  jaxis = np.linspace(0,20,num=fsize,endpoint=True)
  jaxis = jaxis**2 # (j/jrms)^2
  voltot = np.count_nonzero(j2+1e-3)
  # dump the cdf data as table
  fhandler=open(fname+'.vcdf2.tab','a')
  for i in np.arange(fsize):
    jval = jaxis[i]
    cdf  = np.sum(j2[j2>jval])/j2tot
    cdfw2  = np.sum(w2[w2>jval])/j2tot
    #cdfk2  = np.sum(k2[k2>jval])/j2tot
    fvol  = np.count_nonzero(j2[j2>jval])/float(voltot)
    fvolw2  = np.count_nonzero(w2[w2>jval])/float(voltot)
    #fvolk2  = np.count_nonzero(k2[k2>jval])/float(voltot)
    tmp = np.array([jval,cdf,cdfw2,fvol,fvolw2]).reshape(1,5)
    np.savetxt(fhandler,tmp,fmt='%.10e')
  fhandler.close()



  
if __name__=='__main__':
    """
    use this script to derive the viscous dissipation
    with both compressible and incompressible contribution
    """
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
      w2,disp2 = loadData(fname)
      calc_cdf(w2,disp2,fname)
 


# end of the script

