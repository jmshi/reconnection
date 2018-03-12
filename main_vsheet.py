"""
 driver of the whole pipe line
 for finding vorticity sheet
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import athena4_read as ath
import athena_read as ath
import scipy.ndimage.measurements as measurements
import scipy.ndimage.morphology as morphology
import cPickle as pickle
import time
import multiprocessing as mp
from itertools import product
import sys
#from scipy.interpolate import RegularGridInterpolator
from skimage.transform import resize
from collections import defaultdict
from sklearn.decomposition import PCA

def getVsheetProperties(x,y,z,vlist_sorted,nlim=27,eta=2.5e-4):
  """ using input sorted vlist and grid info to calculate
      properties of each current sheet
  """
  nlist = len(vlist_sorted)
  pca = PCA(n_components=3)
  
  # only consider sheet with more than 27 grid points
  nlist = np.max([i for i in range(nlist) if len(vlist_sorted[i])>(nlim-1)])
  
  # allocate list for storage 
  #jmax = []; diss = []; size = []; ncells = []
  dx = x[1]-x[0]
  dy = y[1]-y[0]
  dz = z[1]-z[0]
  dv = dx*dy*dz
  
  def dimensions_vsheet(num):
    js = vlist_sorted[num]
    locs = np.array(js.keys()) 
    xlocs = np.array(zip(z[locs[:,0]],y[locs[:,1]],x[locs[:,2]]))

    if (len(xlocs[:,0])<nlim):
      #print 'reach < ',nlim
      return [0,0,0,(0,0,0)]
    else:
      #(1) num of cells in given sheet
      ncell = len(js)
      #(2) find j_max
      maxj = np.max(np.array(js.values()))
      #(3) compute dissipation \epsilon
      eps  = np.sum(eta*np.array(js.values()))*dv
      #(4) compute (lambda,xi,l)
      # pca.fit(locs)
      pca = PCA(n_components=3).fit(xlocs)
      i1 = pca.components_[0]; i2= pca.components_[1]; i3 = pca.components_[2]
      d1 = np.array([xlocs[i,:].dot(i1) for i in range(0,len(locs[:,0]))])
      d2 = np.array([xlocs[i,:].dot(i2) for i in range(0,len(locs[:,0]))])
      d3 = np.array([xlocs[i,:].dot(i3) for i in range(0,len(locs[:,0]))])
      # dimensions of the sheet
      L1 = np.max([(max(d1)-min(d1)),dx])
      L2 = np.max([(max(d2)-min(d2)),dx])
      L3 = np.max([(max(d3)-min(d3)),dx])
      return [ncell,maxj,eps,(L1,L2,L3),(d1,d2,d3)]
  
  #create a pool and map the target function with multi-arguments
  tstart = time.time() 
  
  npr = 8
  p = mp.Pool(processes=npr)
  result = p.map(dimensions_vsheet,range(nlist),chunksize=npr)
  
  p.close()
  p.join()
  
  tend = time.time()
  print 'time spent with ',str(npr), ' processors = ', tend-tstart, ' seconds'
  
  ncells = [result[i][0] for i in range(nlist)]
  jmax   = [result[i][1] for i in range(nlist)]
  diss   = [result[i][2] for i in range(nlist)]
  size   = [result[i][3] for i in range(nlist)]
  orient = [result[i][4] for i in range(nlist)]

  return ncells,jmax,diss,size,orient
  

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
  return vx,vy,vz,w2


def refine(nx1,nx2,nx3,coarse,factor):
  """refine the input data array by a factor of the grid points"""
  return resize(coarse, (nx3*factor[0], nx2*factor[1], nx1*factor[2]))
 

def estimateJth(j2):
  """
  estimate the j_threshold => 25% of the total 
  viscous dissipation having j>j_threshold;
  return jth (j_threshold^2)
  """
  j2min=np.min(j2)
  j2max=np.max(j2)
  j2th = 0.5*(j2min+j2max)
  diss4 = np.sum(j2)*0.25
  diss  = np.sum(j2[j2>=j2th])
  epsilon = 0.001
  NMAX  = 100
  count = 0
  while ((np.fabs(1.0-diss/diss4) > epsilon) or (count >=NMAX)): 
    count += 1 
    if (diss < diss4): # need to lower jth
      j2max = j2th
      j2th = 0.5*(j2min+j2th)
    else:              # need to raise jth
      j2min = j2th
      j2th = 0.5*(j2th+j2max)
    diss = np.sum(j2[j2>=j2th])
  
  print 'after ',count,' iterations: '
  print 'find threshold vorticity sheet magnitude'
  print 'jth^2 = ',j2th,' 0.25*qdiss_tot = ',diss4,' actual qdiss = ',diss
  return j2th



def findVorticitySheet(w2,fname):
  """
  for given w^2(vorticity squared) find
  the vorticity sheet above w_threshold based 
  on connected-component labeling algorithm.
  """

  print '============ start to find current sheet ==============='
  tstart = time.time()
  jth = estimateJth(w2)
  tend = time.time()
  print 'time cost for the estimate jth: ',tend-tstart,' seconds'

  # save the jth, <vsheet>, and <jtot> in table
  fhandler=open(fname[:-11]+'wth.tab','a')
  frame = int(fname[-11:-6])
  tmp = np.array([frame,jth,np.average(w2[w2>=jth]),np.average(w2)]).reshape(1,4)
  np.savetxt(fhandler,tmp,fmt='%.10e')
  fhandler.close()

  tstart = time.time()
  s=morphology.generate_binary_structure(3,1) # 1: 6-pts connection ; 2: 26(?)-pts connection
  data = np.copy(w2)
  data[data < jth] = 0.
  labeled_array, num_features = measurements.label(data, structure=s)
  tend = time.time()
  print 'identified ',num_features,' vorticity sheets in total!'
  print 'time cost for the identify features: ',tend-tstart,' seconds'

  # construct current sheet list based on the label_array
  #   sort and store each current sheet as a dictionary which 
  #   consists cell location in the cube and current density 
  #   magnitude as key and value. The whole set of current 
  #   sheets then form a list.  
  tstart = time.time()
  whole_dict = defaultdict(list)
  whole_cond = np.argwhere(labeled_array > 0)
  tmid = time.time()
  for loc in whole_cond:
    whole_dict[labeled_array[loc[0],loc[1],loc[2]]].append([(loc[0],loc[1],loc[2]),w2[loc[0],loc[1],loc[2]]])
  #allkeys = [(loc[0],loc[1],loc[2]) for loc in np.argwhere(labeled_array == 1)]
  #alldata = [w2[loc[0],loc[1],loc[2]] for loc in np.argwhere(labeled_array == 1)]
  #ndict = dict(zip(allkeys,alldata))
  tend = time.time()
  print 'time cost: ',tend-tstart,' seconds'
  print 'time cost for argwhere = ', tmid-tstart
  tstart=time.time()
  vlist = []
  for i in range(1,num_features):
    vsheet={}
    for item in whole_dict[i]:
      vsheet[item[0]] = item[1]
    if vsheet:
      vlist.append(vsheet)
    
  tmid = time.time()
  # sort the list based on length
  vlist_sorted = sorted(vlist,key= lambda sheet: -len(sheet))  
  tend = time.time()
  print 'time cost for populating =: ',tmid-tstart,' seconds'
  print 'time cost for sorting    = ', tend-tmid,' seconds'
  print 'first top 5 largest sheets: ',[len(vlist_sorted[i]) for i in range(0,5)]
  print '# of sheets poccess more than 10 cells: ', [i for i in range(num_features-1) if (len(vlist_sorted[i]) >= 10)][-1]
  #tstart = time.time()
  #vlist = []
  #for i in range(1,num_features+1):
  #  vsheet = {}
  #  for loc in (np.argwhere(labeled_array == i)):
  #    vsheet[(loc[0],loc[1],loc[2])] = j2[loc[0],loc[1],loc[2]]
  #  if vsheet: 
  #    vlist.append(vsheet)
  #tend = time.time()
  #print 'time cost for construct vlist: ',tend-tstart,' seconds'
  ## sort the list based on length
  #tstart = time.time()
  #vlist_sorted = sorted(vlist,key = lambda sheet: -len(sheet))
  #tend = time.time()
  #print 'time cost for sort vlist: ',tend-tstart,' seconds'
  #print 'the top 5 sheets in size: ',[len(vlist_sorted[i]) for i in range(0,5)]

  # store dict list with pickle
  # note picle dump has bad performance and we might 
  # need to change the format to hdf5 later
  tstart = time.time()
  pickle.dump(vlist_sorted, open( fname+".vlist.p", "wb" ),2)
  tend = time.time()
  print 'time cost for pickle dump: ',tend-tstart,' seconds'

  



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
      vx,vy,vz,w2 = loadData(fname)
      findVorticitySheet(w2,fname)
 


# end of the script

