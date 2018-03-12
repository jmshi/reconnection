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



def loadData(fname):
  jlist_sorted = pickle.load( open(fname, "rb" ) )
  return jlist_sorted
  
def track(varlist):
  rank = varlist[0] # rank of jsheet under investigation
  js0  = varlist[1] # content of the jsheet at the rank 
  ts,te = varlist[2],varlist[3] # time frames to lookup
  dt = np.pi*0.001
  fdir = '/tigress/jiming/reconnect/athena/bin/'
  basename = 'box.'
  tstride = 1
  size = len(js0)
  js_track = []
  js_track += js0  #set initial jsheet dict
  js_time = []
  js_time += [ts*dt] #set initial jsheet time
  size0 = len(js0)
 
  orig_stdout = sys.stdout
  fname = fdir+targ+'/'+basename+str(ts).zfill(3)+'.r'+str(rank).zfill(1)+'.flow'	    
  fhandler=open(fname,'w') # for print outs
  sys.stdout = fhandler


  nmerg,nself=0,0
  jtrack = {} # dictionary with time as key and individual jsheet as values
  jtrack[ts] = [js0] # initialize tracklist
  flag = -1 # 0: lost 1: merger 2:self-destructor -1: overtime
  for frame in np.arange(ts+1,te+1,tstride):
    fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.vtk.jlist.p'	    
    jlist = loadData(fname)
    js_candid,js_merg,js_self,js_tag =[],[],[],[] #set seg. of jsheet to zero
    ncandid = 0
    js0_size = len(js0)
    jcandid_size = 0
    upper = max([num for num,item in enumerate(jlist) if len(item)>=200]) #use 1000 for stronger cases
    for js in jlist[0:upper]:
      js0_in_js = set(js0).intersection(set(js))
      inside = len(js0_in_js)
      outside = len(js)-inside
      if inside > 0: # put into candidate list
        js_candid += [js]
        ncandid += 1
	jcandid_size +=len(js)
      else:
        pass
    
    print frame, ' ',ncandid,' ',jcandid_size
    #print 'at t= ',frame*dt, 'ncandid= ',ncandid,' candid_size= ',jcandid_size
    if ncandid ==0: # no candidates in the tracking volume print out warning
      #print '[warning] j(rank,size0,size1)= (',rank,',',size0,',',js0_size,') after t= ',(frame-ts)*dt, ' gets lost'
      flag = 0
      break
    else: 
      sorted_candid = sorted(js_candid,key=len)
      lenmax = len(sorted_candid[-1])
      if ncandid >1:
	lensec = len(sorted_candid[-2])
      else:
	lensec = 0
      tmp = np.array([ts,frame,ncandid,lenmax,lensec],dtype='i')
      # 1) major merger
      if jcandid_size > js0_size*1.5:
        #print '[merging] j(',rank,',',size0,',',jcandid_size,') ',tmp 
	flag = 1
        #nmerg += 1
        #jtrack[frame] = js_candid
        js0 = sorted_candid[-1]
	#break
      # 2) self-destructor
      elif (ncandid > 1) and ((jcandid_size < 0.45*size0)): #or (lenmax < 0.85*jcandid_size)):
        #print '[selfdes] j(',rank,',',size0,',',jcandid_size,') ',tmp 
	flag = 2
	js0 = sorted_candid[-1]
      # 3) successor rest
      else: 
        js0 = sorted_candid[-1]	

  #if flag == -1:
    #print '[overtime] j(',rank,',',size0,',',jcandid_size,') ',tmp 


  sys.stdout = orig_stdout
  fhandler.close()
  return
  #print  'rank,size,nself,nmerg,time = ',rank,size,nself,nmerg,frame*dt 
  

if __name__=='__main__':
    if len( sys.argv ) < 4:
        print "Please specify input targ,ts,te,tstride" 
        exit(  )
    targ = sys.argv[1]
    #targ = x2y4z1r64pm1re4000
    ts,te,tstride = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
    dt = np.pi*0.001
    rank = 0
    fdir = '/tigress/jiming/reconnect/athena/bin/'
    basename = 'box.'
    frame=ts
    fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.vtk.jlist.p'	    
    jlist_sorted = loadData(fname)


    print 'ts= ',ts
    p = mp.Pool(processes=5)
    varlist = [[rank,jlist_sorted[rank],ts,te] for rank in range(0,5)]
    result = p.map(track,tuple(varlist))

 #   for js0 in jlist_sorted[4:5]: # only analyze top 10
 #     varlist = [rank,js0,ts,te]
 #     track(varlist)
 #     rank +=1
      #size = len(js0)
      #js_track = []
      #js_track += js0  #set initial jsheet dict
      #js_time = []
      #js_time += [ts*dt] #set initial jsheet time
      #for frame in np.arange(ts+1,te+1,tstride):
      #  fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.vtk.jlist.p'	    
      #  jlist = loadData(fname)
      #  js_self,js_merg,js_tag =[],[],[] #set seg. of jsheet to zero
      #  nself,nmerg = 0,0
      #  for js in jlist[0:100]:
      #    js0_in_js = set(js0).intersection(set(js))
      #    inside = len(js0_in_js)
      #    outside = len(js)-inside
      #    if inside >= outside: # successor found; append list; count+1
      #      js_self += js
      #      js_tag  += [frame*dt]
      #      nself += 1
      #    elif inside >0: #potential merger; start counting
      #      js_merg += js
      #      js_tag += [frame*dt]
      #      nmerg += 1
      #    else:
      #      pass
      #    
      #  if nself == 1 and nmerg == 0: # single successor keep tracking
      #    js_track += js_self
      #    js_time += [dt*frame]
      #    js0 = js_self  # set current jsheet as initial for next step
      #  elif nself == 1 and nmerg > 0: # incoming sheet to merge
      #    js_track += js_self
      #    js_time += js_tag
      #    js_track += js_merg
      #    print 'jsheet = ',rank-1, size,' merged @Dt = ',frame*dt,'nself,nmerg = ',nself,nmerg
      #    break  # break out the ts+1:te+1 loop, go to next init jsheet
      #  elif nself >1: # self-disruption
      #    js_track += js_self
      #    js_track += js_merg
      #    js_time += js_tag
      #    print 'jsheet = ',rank-1, size, 'self-des @ Dt = ',frame*dt, 'nself,nmerg = ',nself,nmerg
      #    break
      #  else:
      #    print '[warning] rank,size,nself,nmerg,time = ',rank-1,size,nself,nmerg,frame*dt 
      #    break

      #print  'rank,size,nself,nmerg,time = ',rank-1,size,nself,nmerg,frame*dt 




# end of the script

