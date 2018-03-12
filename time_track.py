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
  tstamp = varlist[0]
  rank = varlist[1] # rank of jsheet under investigation
  js0  = varlist[2] # content of the jsheet at the rank 
  ts,te = varlist[3],varlist[4] # time frames to lookup
  dt = np.pi*0.001
  fdir = '/tigress/jiming/reconnect/athena/bin/'
  basename = 'box.'
  rank +=1
  tstride = 1
  size = len(js0)
  js_track = []
  js_track += js0  #set initial jsheet dict
  js_time = []
  js_time += [ts*dt] #set initial jsheet time
 
  #fname = fdir+targ+'/'+basename+str(ts).zfill(5)+'.tab'	    
  #fhandler1=open(fname,'a') # for print outs
  for frame in np.arange(ts+1,te+1,tstride):
    fname = fdir+targ+'/'+basename+str(frame).zfill(5)+'.vtk.jlist.p'	    
    jlist = loadData(fname)
    js_self,js_merg,js_tag =[],[],[] #set seg. of jsheet to zero
    nself,nmerg = 0,0
    upper = max([num for num,item in enumerate(jlist) if len(item)>=1000])
    for js in jlist[0:upper]:
      js0_in_js = set(js0).intersection(set(js))
      inside = len(js0_in_js)
      outside = len(js)-inside
      if inside > outside: # successor found; append list; count+1
        js_self += [js]
        js_tag  += [frame*dt]
        nself += 1
      elif inside >0: #potential merger; start counting
        js_merg += [js]
        js_tag += [frame*dt]
        nmerg += 1
      else:
        pass
      
    if js_merg == []:
      lenmergmax = 0 
    else:
      lenmergmax = max([len(item) for item in js_merg])
    if js_self == []:
      lenselfmax = 0 
      lenselfsec = 0
    else:
      sorted_self = sorted(js_self,key=len)
      lenselfmax = len(sorted_self[-1]) #max([len(item) for item in js_self])
      if len(js_self)>=2:
	lenselfsec = len(sorted_self[-2])
      else:
	lenselfsec = 0

    if nself == 1 and nmerg == 0: # single successor keep tracking
      #js_track += js_self
      #js_time += [dt*frame]
      js0 = js_self[0]  # set current jsheet as initial for next step
    elif nself == 1 and nmerg > 0: # incoming sheet to merge
      #js_track += js_self
      #js_time += js_tag
      #js_track += js_merg
      flag = 0
      tmp = np.array([tstamp,rank-1,flag,size,(frame-ts),nself,lenselfmax,lenselfsec,nmerg,lenmergmax],dtype='i')
      #np.savetxt(fhandler1,tmp,fmt='%i %i %i %f %i %i %i %i')
      #print 'jsheet = ',rank-1, size,' merged @Dt = ',frame*dt,'nself,nmerg = ',nself,nmerg
      print tmp
      break  # break out the ts+1:te+1 loop, go to next init jsheet
    elif nself >1: # self-disruption
      #js_track += js_self
      #js_track += js_merg
      #js_time += js_tag
      if lenselfsec >=800: 
        flag = 1
        tmp = np.array([tstamp,rank-1,flag,size,(frame-ts),nself,lenselfmax,lenselfsec,nmerg,lenmergmax],dtype='i')
        #print 'jsheet = ',rank-1, size, 'self-des @ Dt = ',frame*dt, 'nself,nmerg = ',nself,nmerg
        print tmp
        break
      else: 
	js0 = sorted_self[-1] 
    elif nself==0 & nmerg==1: #somehow large displacement occurs
      js0 = js_merg[0]
    else:
      flag = 2
      tmp = np.array([tstamp,rank-1,flag,size,(frame-ts),nself,lenselfmax,lenselfsec,nmerg,lenmergmax],dtype='i')
      #print '[warning] rank,size,nself,nmerg,time = ',rank-1,size,nself,nmerg,frame*dt 
      print tmp
      break

  return
  #print  'rank,size,nself,nmerg,time = ',rank-1,size,nself,nmerg,frame*dt 
  

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
    p = mp.Pool(processes=6)
    varlist = [[ts,rank,jlist_sorted[rank],ts,te] for rank in range(0,6)]
    result = p.map(track,tuple(varlist))

#    for js0 in jlist_sorted[0:2]: # only analyze top 10
#      varlist = [rank,js0,ts,te]
#      track(varlist)
#      rank +=1
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

