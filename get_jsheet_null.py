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
import scipy.ndimage.measurements as measurements
import scipy.ndimage.morphology as morphology
from hessian import hessian

  
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
    idx_maxj = np.argmax(np.array(js.values()))
    #(3) compute cell averaged dissipation \epsilon
    eps  = np.average(eta*np.array(js.values()))
    #(4) compute (lambda,xi,l)
    # pca.fit(locs)
    pca = PCA(n_components=3).fit(xlocs)
    i1 = pca.components_[0]; i2= pca.components_[1]; i3 = pca.components_[2]
    nlocs = len(locs[:,0])
    d1 = np.array([xlocs[i,:].dot(i1) for i in range(0,nlocs)])
    d2 = np.array([xlocs[i,:].dot(i2) for i in range(0,nlocs)])
    #d3 = np.array([xlocs[i,:].dot(i3) for i in range(0,len(locs[:,0]))])
    # dimensions of the sheet
    L1 = np.max([(max(d1)-min(d1)),dx])
    L2 = np.max([(max(d2)-min(d2)),dx])

    ### method 5: using Minkowski functionals where 
    ###           thickness is ~ volume/surface_area
    L3 = 0.0
    volume = ncell*dv
    surface = 0.0
    locset = set(js.keys())
    refset0 = locset.copy() # xy,yz faces
    refset1 = locset.copy() # xz face
    for item in locset:
      (iz,iy,ix) = (item[0],item[1],item[2])
      for i in range(-1,2,2):
	refset0.add((iz,iy,ix+i))
	refset0.add((iz+i,iy,ix))
	refset1.add((iz,iy+i,ix))
    nface0 = len(refset0.difference(locset))
    nface1 = len(refset1.difference(locset))
    surface = nface0*dx*dy+nface1*dx*dz

    L3 = volume*3.0/surface
    #print 'num=',num,' jsheet ncell, volume,surface,L3= ',ncell, volume,surface,L3

    ### method 1: 
    ### projection onto i3: overestimate
    #L3 = np.max([(max(d3)-min(d3)),dx])

    ### method 2:
    ### estimate thickness by volume/area: underestimate as sheet might
    ### be sparse.
    ##vol = ncells*dx**3*2
    ##L3 = np.max([vol/(L1*L2),dx])

    #### method 3:
    #### local projection i3 centered on jmax: still overestimate(?)
    #### similar as method4 except use i3 instead of k
    ##i3 = [1,0,0]
    #rad = np.array([np.linalg.norm(np.cross(xlocs[i,:],i3)) for i in range(0,nlocs)])
    #idx_dis = [i for i in range(0,nlocs) if (np.abs(rad[i]-rad[idx_maxj]) <= 2.0*dx)]
    #ndis = len(idx_dis)
    ##d3 = np.array(sorted([xlocs[i,:].dot(i3) for i in idx_dis]))
    #d3 = np.array(sorted([xlocs[i,:].dot(i3) for i in idx_dis]))
    #d_maxj = xlocs[idx_maxj,:].dot(i3)
    #idx_maxj = (np.abs(d3-d_maxj)).argmin()
    #lower,upper = 0,ndis-1
    #for i in range(idx_maxj,1,-1):
    #  if np.abs(d3[i-1]-d3[i]) >dx*1.414:
    #    lower = i
    #    break
    #for i in range(idx_maxj,ndis-2,1):
    #  if np.abs(d3[i+1]-d3[i]) >dx*1.414:
    #    upper = i
    #    break
    #
    ##L3 = np.max([(max(d3)-min(d3)),dx])
    #L3 = np.abs(d3[upper]-d3[lower])
    #L3 = np.max([L3,dx])

    ### method 4:
    ### local projection of the fastest slope (k) centered on jmax
    ### 0) first get estimate of k using 3d hession
    ### 1) first get distance from the unit vector k for all cells
    ### 2) then pick out all the cells which have distance |d-d(jmax)|
    ### smaller than 2 cells;
    ### 3) separate out the physically connected features within the above
    ### cylinder selection;
    ### 4) do the projection among those connected cells and measure the
    ### maximum distance
#    ncol,nstep = 5,2 # nstep=ncol/2
#    jbox = np.zeros([ncol,ncol,ncol])
#    loclist = [(locs[i,0],locs[i,1],locs[i,2]) for i in range(0,nlocs)]
#    knd,jnd,ind = locs[idx_maxj,0],locs[idx_maxj,1],locs[idx_maxj,2]
#    for i in range(0,ncol):
#      for j in range(0,ncol):
#        for k in range(0,ncol):
#          tnd = (knd-nstep+k,jnd-nstep+j,ind-nstep+i)
#          if tnd in loclist:
#            jbox[k,j,i] = js[tnd]
#    #print 'in jsheet # ',num, ' jbox= ',jbox
#    hess_jbox = hessian(jbox)
#    hess_jmax = hess_jbox[:,:,nstep,nstep,nstep]
#    w,v = np.linalg.eig(hess_jmax)
#    major = 0
#    for i in range(1,3):
#      if np.abs(w[i])>np.abs(w[i-1]):
#        major = i
#    hvec = v[:,major]
#    #print 'hvec = ',hvec
#    ##projection along hvec and then find the thickness 
#    rad = np.array([np.linalg.norm(np.cross(xlocs[i,:],hvec)) for i in range(0,nlocs)])
#    idx_dis = [i for i in range(0,nlocs) if (np.abs(rad[i]-rad[idx_maxj]) <= 1.0*dx)] #2.0*dx)]
#    ndis = len(idx_dis)
#    #d3 = np.array(sorted([xlocs[i,:].dot(i3) for i in idx_dis]))
#    d3 = np.array(sorted([xlocs[i,:].dot(hvec) for i in idx_dis]))
#    d_maxj = xlocs[idx_maxj,:].dot(hvec)
#    idx_maxj = (np.abs(d3-d_maxj)).argmin()
#    lower,upper = 0,ndis-1
#    for i in range(idx_maxj,1,-1):
#      if np.abs(d3[i-1]-d3[i]) >dx*1.414:
#        lower = i
#        break
#    for i in range(idx_maxj,ndis-2,1):
#      if np.abs(d3[i+1]-d3[i]) >dx*1.414:
#        upper = i
#        break
#    
#    #L3 = np.max([(max(d3)-min(d3)),dx])
#    L3 = np.abs(d3[upper]-d3[lower])
#    L3 = np.max([L3,dx])
#    i3 = hvec


    return [ncell,maxj,eps,(L1,L2,L3),(i1,i2,i3)]
  
def dimensions_jnull(num):
  js = jlist_null[num]
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
    idx_maxj = np.argmax(np.array(js.values()))
    #(3) compute cell averaged dissipation \epsilon
    eps  = np.average(eta*np.array(js.values()))
    #(4) compute (lambda,xi,l)
    # pca.fit(locs)
    pca = PCA(n_components=3).fit(xlocs)
    i1 = pca.components_[0]; i2= pca.components_[1]; i3 = pca.components_[2]
    nlocs = len(locs[:,0])
    d1 = np.array([xlocs[i,:].dot(i1) for i in range(0,nlocs)])
    d2 = np.array([xlocs[i,:].dot(i2) for i in range(0,nlocs)])
    #d3 = np.array([xlocs[i,:].dot(i3) for i in range(0,len(locs[:,0]))])
    # dimensions of the sheet
    L1 = np.max([(max(d1)-min(d1)),dx])
    L2 = np.max([(max(d2)-min(d2)),dx])

    ### method 5: using Minkowski functionals where 
    ###           thickness is ~ volume/surface_area
    L3 = 0.0
    volume = ncell*dv
    surface = 0.0
    locset = set(js.keys())
    refset0 = locset.copy() # xy,yz faces
    refset1 = locset.copy() # xz face
    for item in locset:
      (iz,iy,ix) = (item[0],item[1],item[2])
      for i in range(-1,2,2):
	refset0.add((iz,iy,ix+i))
	refset0.add((iz+i,iy,ix))
	refset1.add((iz,iy+i,ix))
    nface0 = len(refset0.difference(locset))
    nface1 = len(refset1.difference(locset))
    surface = nface0*dx*dy+nface1*dx*dz

    L3 = volume*3.0/surface
    return [ncell,maxj,eps,(L1,L2,L3),(i1,i2,i3)]

def get_xcenter(x1f,x2f,x3f):
  x = x1f+0.5*(x1f[1]-x1f[0])
  y = x2f+0.5*(x2f[1]-x2f[0])
  z = x3f+0.5*(x3f[1]-x3f[0])
  return x[:-1],y[:-1],z[:-1]
  
def loadData(fname):
  jlist_sorted = pickle.load( open(fname, "rb" ) )
  return jlist_sorted

def loadNull(fname):
  null = pickle.load( open(fname, "rb" ) )
  return null

def identifier(jlist,nlist):
  j_null=[]
  #j_nonu = []
  null_copy = nlist[:]
  for i in range(0,len(jlist)):
    points = set(null_copy).intersection(jlist[i].keys()) 
    if bool(points): # found identity 
      j_null.append(jlist[i])
    #else:  # no identity 
    #  j_nonu +=jlist[i]
    null_copy=[x for x in null_copy if x not in list(points)]

  return j_null #, j_nonu

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
    # load jlist
    fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.athdf.jlist.p'	    
    jlist_sorted = loadData(fname)
    nlist = len(jlist_sorted)
    pca = PCA(n_components=3)
    # load null pts 
    fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.athdf.null.p'	    
    null = loadNull(fname)

    #print 'len(jlist_sorted) = ',len(jlist_sorted)
    #print 'len(null) = ',len(null)
    # separate sheet with null-pts and w/o null-pts
    jlist_null = identifier(jlist_sorted,null)

    # only consider sheet with more than nlim grid points
    nlist = len(jlist_null)
    nlist = np.max([i for i in range(nlist) if len(jlist_null[i])>=nlim])
    print 'nlist=',nlist
    
    #create a pool and map the target function with multi-arguments
    tstart = time.time() 
    npr = 8
    p = mp.Pool(processes=npr)
    result = p.map(dimensions_jnull,range(nlist),chunksize=npr)
  
    p.close()
    p.join()
  
    tend = time.time()
    print 'time spent with ',str(npr), ' processors = ', tend-tstart, ' seconds'
  
    ncells0 = [result[i][0] for i in range(nlist)]
    jmax0   = [result[i][1] for i in range(nlist)]
    diss0   = [result[i][2] for i in range(nlist)]
    size0   = [result[i][3] for i in range(nlist)]
    orient0 = [result[i][4] for i in range(nlist)]
    pname = fdir+targ+'/'+'jprop_null.'+str(frame).zfill(5)+'.p'
    pickle.dump([ncells0,jmax0,diss0,size0,orient0], open(pname, "wb" ),2)

    ncells += ncells0 #[result[i][0] for i in range(nlist)]
    jmax   += jmax0   #[result[i][1] for i in range(nlist)]
    diss   += diss0   #[result[i][2] for i in range(nlist)]
    size   += size0   #[result[i][3] for i in range(nlist)]
    orient += orient0 #[result[i][4] for i in range(nlist)]

  pname = fdir+targ+'/'+'jprop_null_ts='+str(ts)+'_te='+str(te)+'.p'
  pickle.dump([ncells,jmax,diss,size,orient], open(pname, "wb" ),2)
 


# end of the script

