"""
 driver of the whole pipe line
 for derive jsheet properties
"""
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
#import athena4_read as ath
import athena_read as ath
import cPickle as pickle
import time
import multiprocessing as mp
from itertools import product
import sys
from sklearn.decomposition import PCA

  
def dimensions_jsheet(num):
  js = jlist_sorted[num]
  locs = np.array(js.keys()) 
  xlocs = np.array(zip(z[locs[:,0]],y[locs[:,1]],x[locs[:,2]]))

  if (len(xlocs[:,0])<nlim):
    #print 'reach < ',nlim
    return [0,0,0,(0,0,0),(0,0,0)]
  else:
    #(1) num of cells in given sheet
    ncell = len(js)
    #(2) find j_max
    maxj = np.sqrt(np.max(np.array(js.values())))
    xmaxj = xlocs[np.argmax(np.array(js.values())),:]
    #(3) compute cell averaged dissipation \epsilon
    eps  = np.average(eta*np.array(js.values()))
    #(4) compute (lambda,xi,l)
    # pca.fit(locs)
    pca = PCA(n_components=3).fit(xlocs)
    i1 = pca.components_[0]; i2= pca.components_[1]; i3 = pca.components_[2]
    nlocs = len(locs[:,0])
    d1 = np.array([xlocs[i,:].dot(i1) for i in range(0,nlocs)])
    d2 = np.array([xlocs[i,:].dot(i2) for i in range(0,nlocs)])
    d3 = np.array([xlocs[i,:].dot(i3) for i in range(0,len(locs[:,0]))])
    # dimensions of the sheet
    L1 = np.max([(max(d1)-min(d1)),dx])
    L2 = np.max([(max(d2)-min(d2)),dx])
    L3 = np.max([(max(d3)-min(d3)),dx])
    return [ncell,maxj,eps,(L1,L2,L3),(i1,i2,i3)]
  
def get_xcenter(x1f,x2f,x3f):
  x = x1f+0.5*(x1f[1]-x1f[0])
  y = x2f+0.5*(x2f[1]-x2f[0])
  z = x3f+0.5*(x3f[1]-x3f[0])
  return x[:-1],y[:-1],z[:-1]
  
def loadData(fname):
  jlist_sorted = pickle.load( open(fname, "rb" ) )
  return jlist_sorted


if __name__=='__main__':
  if len( sys.argv ) < 4:
    print "calling sequence: python get_jsheet_prop.py targ,ts,te,tstride,[eta=2.5e-4],[nlim=8]" 
    exit(  )
  targ = sys.argv[1]
  #targ = x2y4z1r64pm1re4000
  ts,te,tstride = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
  fdir = '/tigress/jiming/reconnect/athena/bin/'
  basename = 'Unstra.out2.'
  t,grid = ath.athdf(fdir+targ+'/'+basename+str(ts).zfill(5)+'.athdf',quantities=['x1f','x2f','x3f'])	    
  x,y,z = get_xcenter(grid['x1f'],grid['x2f'],grid['x3f'])
  dx = x[1]-x[0]
  dy = y[1]-y[0]
  dz = z[1]-z[0]
  dv = dx*dy*dz
  if len(sys.argv)>=6:
    eta = float(sys.argv[5])
  else:
    eta = 2.5e-4
  if len(sys.argv)>=7:
    nlim = int(sys.argv[6])
  else:
    nlim = 8
  print 'nlim = ',nlim
  print 'eta  = ',eta

  ncells=[]; jmax=[];diss=[];size=[];orient=[]
  for frame in np.arange(ts,te,tstride):
    fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.athdf.jlist.p'	    
    jlist_sorted = loadData(fname)
    nlist = len(jlist_sorted)
    pca = PCA(n_components=3)
  
    # only consider sheet with more than nlim grid points
    nlist = np.max([i for i in range(nlist) if len(jlist_sorted[i])>=nlim])
    
    #create a pool and map the target function with multi-arguments
    tstart = time.time() 
    npr = 8
    p = mp.Pool(processes=npr)
    result = p.map(dimensions_jsheet,range(nlist),chunksize=npr)
  
    p.close()
    p.join()
  
    tend = time.time()
    print 'time spent with ',str(npr), ' processors = ', tend-tstart, ' seconds'
  
    ncells0 = [result[i][0] for i in range(nlist)]
    jmax0   = [result[i][1] for i in range(nlist)]
    diss0   = [result[i][2] for i in range(nlist)]
    size0   = [result[i][3] for i in range(nlist)]
    orient0 = [result[i][4] for i in range(nlist)]
    pname = fdir+targ+'/'+'jprop.'+str(frame).zfill(5)+'.p'
    pickle.dump([ncells0,jmax0,diss0,size0,orient0], open(pname, "wb" ),2)

    ncells += ncells0 #[result[i][0] for i in range(nlist)]
    jmax   += jmax0   #[result[i][1] for i in range(nlist)]
    diss   += diss0   #[result[i][2] for i in range(nlist)]
    size   += size0   #[result[i][3] for i in range(nlist)]
    orient += orient0 #[result[i][4] for i in range(nlist)]

  pname = fdir+targ+'/'+'jprop_ts='+str(ts)+'_te='+str(te)+'.p'
  pickle.dump([ncells,jmax,diss,size,orient], open(pname, "wb" ),2)
 


# end of the script

