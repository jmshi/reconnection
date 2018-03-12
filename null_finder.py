"""
 driver of the whole pipe line
 for finding current sheet and null pts
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



## the following func is moved to get_jsheet_prop.py_
def getJsheetProperties(x,y,z,jlist_sorted,nlim=27,eta=2.5e-4):
  """ using input sorted jlist and grid info to calculate
      properties of each current sheet
  """
  nlist = len(jlist_sorted)
  pca = PCA(n_components=3)
  
  # only consider sheet with more than 27 grid points
  nlist = np.max([i for i in range(nlist) if len(jlist_sorted[i])>(nlim-1)])
  
  # allocate list for storage 
  #jmax = []; diss = []; size = []; ncells = []
  dx = x[1]-x[0]
  dy = y[1]-y[0]
  dz = z[1]-z[0]
  dv = dx*dy*dz
  
  def dimensions_jsheet(num):
    js = jlist_sorted[num]
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
  result = p.map(dimensions_jsheet,range(nlist),chunksize=npr)
  
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
  time,data=ath.athdf(fname,quantities=['Bcc1'])
  bx = data['Bcc1']
  time,data=ath.athdf(fname,quantities=['Bcc2'])
  by = data['Bcc2']
  time,data=ath.athdf(fname,quantities=['Bcc3'])
  bz = data['Bcc3']
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
  jx,jy,jz = curl(bx,by,bz,dx,dy,dz)
  j2 = jx**2+jy**2+jz**2
  return bx,by,bz,j2

def loadBonly(fname='Unstra.out2.00008.athdf'):
  """load 3d bfield only"""
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
  return bx,by,bz

def refine(nx1,nx2,nx3,coarse,factor):
  """refine the input data array by a factor of the grid points"""
  return resize(coarse, (nx3*factor[0], nx2*factor[1], nx1*factor[2]))
 

def estimateJth(j2):
  """
  estimate the j_threshold => 25% of the total 
  ohmic dissipation having j>j_threshold;
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
  print 'find threshold current sheet magnitude'
  print 'jth^2 = ',j2th,' 0.25*qdiss_tot = ',diss4,' actual qdiss = ',diss
  return j2th



def findCurrentSheet(j2,fname):
  """
  for given j^2(current density squared) find
  the current sheet above j_threshold based 
  on connected-component labeling algorithm.
  """

  print '============ start to find current sheet ==============='
  tstart = time.time()
  jth = estimateJth(j2)
  tend = time.time()
  print 'time cost for the estimate jth: ',tend-tstart,' seconds'

  # save the jth, <jsheet>, and <jtot> in table
  fhandler=open(fname[:-11]+'jth.tab','a')
  frame = int(fname[-11:-6])
  tmp = np.array([frame,jth,np.average(j2[j2>=jth]),np.average(j2)]).reshape(1,4)
  np.savetxt(fhandler,tmp,fmt='%.10e')
  fhandler.close()

  tstart = time.time()
  s=morphology.generate_binary_structure(3,1) # 1: 6-pts connection ; 2: 34(?)-pts connection
  data = np.copy(j2)
  data[data < jth] = 0.
  labeled_array, num_features = measurements.label(data, structure=s)
  tend = time.time()
  print 'identified ',num_features,' current sheets in total!'
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
    whole_dict[labeled_array[loc[0],loc[1],loc[2]]].append([(loc[0],loc[1],loc[2]),j2[loc[0],loc[1],loc[2]]])
  #allkeys = [(loc[0],loc[1],loc[2]) for loc in np.argwhere(labeled_array == 1)]
  #alldata = [j2[loc[0],loc[1],loc[2]] for loc in np.argwhere(labeled_array == 1)]
  #ndict = dict(zip(allkeys,alldata))
  tend = time.time()
  print 'time cost: ',tend-tstart,' seconds'
  print 'time cost for argwhere = ', tmid-tstart
  tstart=time.time()
  jlist = []
  for i in range(1,num_features):
    jsheet={}
    for item in whole_dict[i]:
      jsheet[item[0]] = item[1]
    if jsheet:
      jlist.append(jsheet)
    
  tmid = time.time()
  # sort the list based on length
  jlist_sorted = sorted(jlist,key= lambda sheet: -len(sheet))  
  tend = time.time()
  print 'time cost for populating =: ',tmid-tstart,' seconds'
  print 'time cost for sorting    = ', tend-tmid,' seconds'
  print 'first top 5 largest sheets: ',[len(jlist_sorted[i]) for i in range(0,5)]
  print '# of sheets poccess more than 10 cells: ', [i for i in range(num_features-1) if (len(jlist_sorted[i]) >= 10)][-1]
  #tstart = time.time()
  #jlist = []
  #for i in range(1,num_features+1):
  #  jsheet = {}
  #  for loc in (np.argwhere(labeled_array == i)):
  #    jsheet[(loc[0],loc[1],loc[2])] = j2[loc[0],loc[1],loc[2]]
  #  if jsheet: 
  #    jlist.append(jsheet)
  #tend = time.time()
  #print 'time cost for construct jlist: ',tend-tstart,' seconds'
  ## sort the list based on length
  #tstart = time.time()
  #jlist_sorted = sorted(jlist,key = lambda sheet: -len(sheet))
  #tend = time.time()
  #print 'time cost for sort jlist: ',tend-tstart,' seconds'
  #print 'the top 5 sheets in size: ',[len(jlist_sorted[i]) for i in range(0,5)]

  # store dict list with pickle
  # note picle dump has bad performance and we might 
  # need to change the format to hdf5 later
  tstart = time.time()
  pickle.dump(jlist_sorted, open( fname+".jlist.p", "wb" ),2)
  tend = time.time()
  print 'time cost for pickle dump: ',tend-tstart,' seconds'

  

def positive_negative(inputlist):
  """ return true if all items in the input list having same sign (+ or -)"""
  if ((all(item>0 for item in list(inputlist))) or (all(item<0 for item in list(inputlist)))):
    return True  # elements with same sign (+ or -)
  else:
    return False # elements with diff sign or all zeros;

def f(location):
  """check the corner values of b-field
     return location if not sharing same signs"""
  k = location[0];j=location[1];i=location[2]
  corner=bx[k:k+2,j:j+2,i:i+2]
  is_bx = positive_negative(np.ravel(corner))
  corner=by[k:k+2,j:j+2,i:i+2]
  is_by = positive_negative(np.ravel(corner))
  corner=bz[k:k+2,j:j+2,i:i+2]
  is_bz = positive_negative(np.ravel(corner))
  if ((is_bx or is_by or is_bz) == False):
    return (k,j,i) # find candidate
  else:
    return None    # not candidate

def bilinear(x,y,f00,f01,f10,f11):
  """ function which takes the cell location and solve for 
      intersection points on six surfaces """
  a = f00; b = f10-f00; c = f01-f00; d = f00+f11-f10-f01
  fxy = a + b*x + c*y + d*x*y
  return fxy

def issym(b3):
  """test if a list has equal number of positive
     and negative values; zeros belong to both. """
  npos = 0; nneg = 0
  for item in b3:
    if (item >= 0): 
      npos +=1
    if (item <= 0):
      nneg +=1
  if (npos==nneg):
    return True
  else:
    return False

###############################################
# intersection curves Bx=By=0
###############################################
def intersect_bxby(location):
  npts = 0 # # of end points
  b3 = [] # third component values
  sk = location[0];sj=location[1];si=location[2]
  f = bx; g = by; h = bz
    
  for nface in range(6):
    if (nface == 0):  # x = 0 face (y,z) 
      f00,f01,f10,f11 = f[sk,sj,si],f[sk+1,sj,si],f[sk,sj+1,si],f[sk+1,sj+1,si]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk+1,sj,si],g[sk,sj+1,si],g[sk+1,sj+1,si]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk+1,sj,si],h[sk,sj+1,si],h[sk+1,sj+1,si]
    if (nface == 1): # x=1 face 
      f00,f01,f10,f11 = f[sk,sj,si+1],f[sk+1,sj,si+1],f[sk,sj+1,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj,si+1],g[sk+1,sj,si+1],g[sk,sj+1,si+1],g[sk+1,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj,si+1],h[sk+1,sj,si+1],h[sk,sj+1,si+1],h[sk+1,sj+1,si+1]
    if (nface == 2):  # y = 0 face (x,z) 
      f00,f01,f10,f11 = f[sk,sj,si],f[sk+1,sj,si],f[sk,sj,si+1],f[sk+1,sj,si+1]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk+1,sj,si],g[sk,sj,si+1],g[sk+1,sj,si+1]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk+1,sj,si],h[sk,sj,si+1],h[sk+1,sj,si+1]
    if (nface == 3): # y =1 face
      f00,f01,f10,f11 = f[sk,sj+1,si],f[sk+1,sj+1,si],f[sk,sj+1,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj+1,si],g[sk+1,sj+1,si],g[sk,sj+1,si+1],g[sk+1,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj+1,si],h[sk+1,sj+1,si],h[sk,sj+1,si+1],h[sk+1,sj+1,si+1]
    if (nface == 4): # z = 0 face (x,y)
      f00,f01,f10,f11 = f[sk,sj,si],f[sk,sj+1,si],f[sk,sj,si+1],f[sk,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk,sj+1,si],g[sk,sj,si+1],g[sk,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk,sj+1,si],h[sk,sj,si+1],h[sk,sj+1,si+1]
    if (nface == 5): # z=1 face
      f00,f01,f10,f11 = f[sk+1,sj,si],f[sk+1,sj+1,si],f[sk+1,sj,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk+1,sj,si],g[sk+1,sj+1,si],g[sk+1,sj,si+1],g[sk+1,sj+1,si+1] 
      h00,h01,h10,h11 = h[sk+1,sj,si],h[sk+1,sj+1,si],h[sk+1,sj,si+1],h[sk+1,sj+1,si+1]
        
    a1,b1,c1,d1 = f00,f10-f00,f01-f00,f00+f11-f10-f01
    a2,b2,c2,d2 = g00,g10-g00,g01-g00,g00+g11-g10-g01
    coeff = [(b1*d2-b2*d1), (a1*d2-a2*d1+b1*c2-b2*c1), (a1*c2-a2*c1)]
    sol1 = np.roots(coeff)
    sol2 = -(a2+b2*sol1)/(c2+d2*sol1)
    #  determine if intersection points
    for s,t in zip(sol1,sol2):
       if ((not isinstance(s,complex)) and (not isinstance(t,complex))): 
         if (s>=0 and s<=1 and t>=0 and t<=1):
           npts += 1
           # if yes, calculate the third components on those points
           b3.append(bilinear(s,t,h00,h01,h10,h11))

  #  determine if null pt: (1) intersection points in pair (2) third component flip signs
  #  if yes then return True; else return False
  if (npts>0 and npts%2==0): 
    return issym(b3) 
  else:
    return False

###############################################
# intersection curves By=Bz=0
###############################################
def intersect_bybz(location):
  npts = 0 # # of end points
  b3 = [] # third component values
  sk = location[0];sj=location[1];si=location[2]
  f = by; g = bz; h = bx
    
  for nface in range(6):
    if (nface == 0):  # x = 0 face (y,z) 
      f00,f01,f10,f11 = f[sk,sj,si],f[sk+1,sj,si],f[sk,sj+1,si],f[sk+1,sj+1,si]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk+1,sj,si],g[sk,sj+1,si],g[sk+1,sj+1,si]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk+1,sj,si],h[sk,sj+1,si],h[sk+1,sj+1,si]
    if (nface == 1): # x=1 face 
      f00,f01,f10,f11 = f[sk,sj,si+1],f[sk+1,sj,si+1],f[sk,sj+1,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj,si+1],g[sk+1,sj,si+1],g[sk,sj+1,si+1],g[sk+1,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj,si+1],h[sk+1,sj,si+1],h[sk,sj+1,si+1],h[sk+1,sj+1,si+1]
    if (nface == 2):  # y = 0 face (x,z) 
      f00,f01,f10,f11 = f[sk,sj,si],f[sk+1,sj,si],f[sk,sj,si+1],f[sk+1,sj,si+1]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk+1,sj,si],g[sk,sj,si+1],g[sk+1,sj,si+1]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk+1,sj,si],h[sk,sj,si+1],h[sk+1,sj,si+1]
    if (nface == 3): # y =1 face
      f00,f01,f10,f11 = f[sk,sj+1,si],f[sk+1,sj+1,si],f[sk,sj+1,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj+1,si],g[sk+1,sj+1,si],g[sk,sj+1,si+1],g[sk+1,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj+1,si],h[sk+1,sj+1,si],h[sk,sj+1,si+1],h[sk+1,sj+1,si+1]
    if (nface == 4): # z = 0 face (x,y)
      f00,f01,f10,f11 = f[sk,sj,si],f[sk,sj+1,si],f[sk,sj,si+1],f[sk,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk,sj+1,si],g[sk,sj,si+1],g[sk,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk,sj+1,si],h[sk,sj,si+1],h[sk,sj+1,si+1]
    if (nface == 5): # z=1 face
      f00,f01,f10,f11 = f[sk+1,sj,si],f[sk+1,sj+1,si],f[sk+1,sj,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk+1,sj,si],g[sk+1,sj+1,si],g[sk+1,sj,si+1],g[sk+1,sj+1,si+1] 
      h00,h01,h10,h11 = h[sk+1,sj,si],h[sk+1,sj+1,si],h[sk+1,sj,si+1],h[sk+1,sj+1,si+1]
        
    a1,b1,c1,d1 = f00,f10-f00,f01-f00,f00+f11-f10-f01
    a2,b2,c2,d2 = g00,g10-g00,g01-g00,g00+g11-g10-g01
    coeff = [(b1*d2-b2*d1), (a1*d2-a2*d1+b1*c2-b2*c1), (a1*c2-a2*c1)]
    sol1 = np.roots(coeff)
    sol2 = -(a2+b2*sol1)/(c2+d2*sol1)
    #  determine if intersection points
    for s,t in zip(sol1,sol2):
       if ((not isinstance(s,complex)) and (not isinstance(t,complex))): 
         if (s>=0 and s<=1 and t>=0 and t<=1):
           npts += 1
           # if yes, calculate the third components on those points
           b3.append(bilinear(s,t,h00,h01,h10,h11))

  #  determine if null pt: (1) intersection points in pair (2) third component flip signs
  #  if yes then return True; else return False
  if (npts>0 and npts%2==0): 
    return issym(b3) 
  else:
    return False

###############################################
# intersection curves Bx=Bz=0
###############################################
def intersect_bxbz(location):
  npts = 0 # # of end points
  b3 = [] # third component values
  sk = location[0];sj=location[1];si=location[2]
  f = bx; g = bz; h = by
    
  for nface in range(6):
    if (nface == 0):  # x = 0 face (y,z) 
      f00,f01,f10,f11 = f[sk,sj,si],f[sk+1,sj,si],f[sk,sj+1,si],f[sk+1,sj+1,si]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk+1,sj,si],g[sk,sj+1,si],g[sk+1,sj+1,si]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk+1,sj,si],h[sk,sj+1,si],h[sk+1,sj+1,si]
    if (nface == 1): # x=1 face 
      f00,f01,f10,f11 = f[sk,sj,si+1],f[sk+1,sj,si+1],f[sk,sj+1,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj,si+1],g[sk+1,sj,si+1],g[sk,sj+1,si+1],g[sk+1,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj,si+1],h[sk+1,sj,si+1],h[sk,sj+1,si+1],h[sk+1,sj+1,si+1]
    if (nface == 2):  # y = 0 face (x,z) 
      f00,f01,f10,f11 = f[sk,sj,si],f[sk+1,sj,si],f[sk,sj,si+1],f[sk+1,sj,si+1]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk+1,sj,si],g[sk,sj,si+1],g[sk+1,sj,si+1]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk+1,sj,si],h[sk,sj,si+1],h[sk+1,sj,si+1]
    if (nface == 3): # y =1 face
      f00,f01,f10,f11 = f[sk,sj+1,si],f[sk+1,sj+1,si],f[sk,sj+1,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj+1,si],g[sk+1,sj+1,si],g[sk,sj+1,si+1],g[sk+1,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj+1,si],h[sk+1,sj+1,si],h[sk,sj+1,si+1],h[sk+1,sj+1,si+1]
    if (nface == 4): # z = 0 face (x,y)
      f00,f01,f10,f11 = f[sk,sj,si],f[sk,sj+1,si],f[sk,sj,si+1],f[sk,sj+1,si+1]
      g00,g01,g10,g11 = g[sk,sj,si],g[sk,sj+1,si],g[sk,sj,si+1],g[sk,sj+1,si+1]
      h00,h01,h10,h11 = h[sk,sj,si],h[sk,sj+1,si],h[sk,sj,si+1],h[sk,sj+1,si+1]
    if (nface == 5): # z=1 face
      f00,f01,f10,f11 = f[sk+1,sj,si],f[sk+1,sj+1,si],f[sk+1,sj,si+1],f[sk+1,sj+1,si+1]
      g00,g01,g10,g11 = g[sk+1,sj,si],g[sk+1,sj+1,si],g[sk+1,sj,si+1],g[sk+1,sj+1,si+1] 
      h00,h01,h10,h11 = h[sk+1,sj,si],h[sk+1,sj+1,si],h[sk+1,sj,si+1],h[sk+1,sj+1,si+1]
        
    a1,b1,c1,d1 = f00,f10-f00,f01-f00,f00+f11-f10-f01
    a2,b2,c2,d2 = g00,g10-g00,g01-g00,g00+g11-g10-g01
    coeff = [(b1*d2-b2*d1), (a1*d2-a2*d1+b1*c2-b2*c1), (a1*c2-a2*c1)]
    sol1 = np.roots(coeff)
    sol2 = -(a2+b2*sol1)/(c2+d2*sol1)
    #  determine if intersection points
    for s,t in zip(sol1,sol2):
       if ((not isinstance(s,complex)) and (not isinstance(t,complex))): 
         if (s>=0 and s<=1 and t>=0 and t<=1):
           npts += 1
           # if yes, calculate the third components on those points
           b3.append(bilinear(s,t,h00,h01,h10,h11))

  #  determine if null pt: (1) intersection points in pair (2) third component flip signs
  #  if yes then return True; else return False
  if (npts>0 and npts%2==0): 
    return issym(b3) 
  else:
    return False



def findNullPoint(bx,by,bz,fname, nproc=0):
  """
  for given magnetic field, find all the neutral
  points (null points) based on the trilinear 
  method (Haynes & Parnell 2007, Phys of Plasmas
  14, 082107)
  """
  print '============ start to find null points ==============='
  tstart = time.time()  
  nx = bx.shape[2]
  ny = bx.shape[1]
  nz = bx.shape[0] 
  #contruct list of argument for parallelization
  klist,jlist,ilist = zip(*product(range(nz-1),range(ny-1),range(nx-1)))
  tend = time.time()
  print 'time cost for constructing location: ',tend-tstart,' seconds'
  

  # (1) reduction process
  tstart = time.time()  
  result_list=[]
  candidate_list=[]
  #create a pool and map the target function with multi-arguments
  p = mp.Pool(processes=nproc)
  result_list = p.map(f,zip(klist,jlist,ilist))
  tend = time.time()
  print 'time cost for reduction [step1]: ',tend-tstart,' seconds'
 
  tstart = time.time()  
  candidate_list = [item for item in result_list if item is not None]
  tend = time.time()
  print 'time cost for reduction [step2]: ',tend-tstart,' seconds'
  del result_list
  #p.terminate()
  #p.close()

  print 'identified ',len(candidate_list), ' candidate cells after reduction'
 

  # (2) selection process (by counting intersection points and 
  #     parity of the third component)
  tstart = time.time()  
  # first pass: finding intersection of Bx=By=0

  #create a pool and map the target function with multi-arguments
  #p = mp.Pool(processes=nproc)
  mask_list = p.map(intersect_bxby, candidate_list)
  cell = [candidate_list[i] for i in range(len(candidate_list)) if mask_list[i]]
  tend = time.time()
  print 'time cost for Bx=By=0: ',tend-tstart,' seconds'
  print len(cell),' cells left'
  
  tstart = time.time() 
  # second pass: finding intersection of By=Bz=0
  candidate_list = cell
  mask_list = p.map(intersect_bxby, candidate_list)
  cell = [candidate_list[i] for i in range(len(candidate_list)) if mask_list[i]]
  tend = time.time()
  #for cell_one in cell:
  #  if (not intersect_bybz(cell_one)):
  #    cell.remove(cell_one)
  #tend = time.time()
  print 'time cost for By=Bz=0: ',tend-tstart,' seconds'
  print len(cell),' cells left'
  
  tstart = time.time() 
  # third pass: finding intersection of Bx=Bz=0
  candidate_list = cell
  mask_list = p.map(intersect_bxby, candidate_list)
  cell = [candidate_list[i] for i in range(len(candidate_list)) if mask_list[i]]
  tend = time.time()
  #for cell_one in cell:
  #  if (not intersect_bxbz(cell_one)):
  #    cell.remove(cell_one) 
  #tend = time.time()
  print 'time cost for the Bx=Bz=0: ',tend-tstart,' seconds'
  print len(cell),' cells left'

  print 'identified ',len(cell),' magnetic null points'

  p.close()


  #for item in cell:
  #    a,b,c = item[2],item[1],item[0]
  #    print x[a],y[b],z[c]
      
  # dump the magnetic null
  #pickle.dump(cell, open( "NULL.p", "wb" ) )

  pickle.dump(cell, open( fname+".null.p", "wb" ),2)

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
      bx,by,bz = loadBonly(fname)
      #nx,ny,nz = len(x),len(y),len(z)
      #factor = (4,4,4)
      #tstart = time.time()
      #bx = refine(nx,ny,nz,bx,factor)
      #by = refine(nx,ny,nz,by,factor)
      #bz = refine(nx,ny,nz,bz,factor)
      #j2 = refine(nx,ny,nz,j2,factor)
      #tend   = time.time()
      #print 'time cost for ', factor,' times refinement: ',tend-tstart,' seconds'
      #print 'now the input data shape: ',bx.shape,by.shape,bz.shape,j2.shape
      #del bx,by,bz,j2
      #return
      #findCurrentSheet(j2,fname)
      findNullPoint(bx,by,bz,fname,16)
 


# end of the script

