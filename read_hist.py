import numpy as np


##
# history dump reader
#####

def read_hist1d(fname,adb=False):
    def time_trim(ahist):
      count = 0; tmp = np.copy(ahist)
      for i in np.arange(1,len(ahist['time'])):
        if ahist['time'][i] > tmp['time'][count]:
          tmp['time'][count+1] = ahist['time'][i]
          tmp['em1'][count+1] = ahist['em1'][i]
          tmp['em2'][count+1] = ahist['em2'][i]
          tmp['em3'][count+1] = ahist['em3'][i]
          tmp['maxwell'][count+1] = ahist['maxwell'][i]
          tmp['reynolds'][count+1] = ahist['reynolds'][i]
          count +=1
      #print "finishing counting count= ",count
      tt  = tmp['time'][0:count]
      em1 = tmp['em1'][0:count]
      em2 = tmp['em2'][0:count]
      em3 = tmp['em3'][0:count]
      maxwell  = tmp['maxwell'][0:count]
      reynolds = tmp['reynolds'][0:count]
      return tt,em1,em2,em3,maxwell,reynolds

    dtype = np.dtype([('time', 'f8'), ('em1', 'f8'),('em2','f8'),('em3','f8'),('maxwell','f8'),('reynolds','f8')])
    if adb:
      ahist = np.loadtxt(fname, dtype=dtype, skiprows=2, usecols=(0,10,11,12,13,14))
    else:
      ahist = np.loadtxt(fname, dtype=dtype, skiprows=2, usecols=(0,9,10,11,12,13))
    tt,em1,em2,em3,maxwell,reynolds = time_trim(ahist)
    return tt,em1,em2,em3,maxwell,reynolds

def read_hist1d_ek(fname,adb=False):
    def time_trim(ahist):
      count = 0; tmp = np.copy(ahist)
      for i in np.arange(1,len(ahist['time'])):
        if ahist['time'][i] > tmp['time'][count]:
          tmp['time'][count+1] = ahist['time'][i]
          tmp['ek1'][count+1] = ahist['ek1'][i]
          tmp['ek2'][count+1] = ahist['ek2'][i]
          tmp['ek3'][count+1] = ahist['ek3'][i]
          tmp['maxwell'][count+1] = ahist['maxwell'][i]
          tmp['reynolds'][count+1] = ahist['reynolds'][i]
          count +=1
      #print "finishing counting count= ",count
      tt  = tmp['time'][0:count]
      ek1 = tmp['ek1'][0:count]
      ek2 = tmp['ek2'][0:count]
      ek3 = tmp['ek3'][0:count]
      maxwell  = tmp['maxwell'][0:count]
      reynolds = tmp['reynolds'][0:count]
      return tt,ek1,ek2,ek3,maxwell,reynolds

    dtype = np.dtype([('time', 'f8'), ('ek1', 'f8'),('ek2','f8'),('ek3','f8'),('maxwell','f8'),('reynolds','f8')])
    if adb:
      ahist = np.loadtxt(fname, dtype=dtype, skiprows=2, usecols=(0,6,7,8,13,14))
    else:
      ahist = np.loadtxt(fname, dtype=dtype, skiprows=2, usecols=(0,6,7,8,12,13))
    tt,ek1,ek2,ek3,maxwell,reynolds = time_trim(ahist)
    return tt,ek1,ek2,ek3,maxwell,reynolds

def read_hist1d_energy(fname):
    def time_trim(ahist):
      count = 0; tmp = np.copy(ahist)
      for i in np.arange(1,len(ahist['time'])):
        if ahist['time'][i] > tmp['time'][count]:
          tmp['time'][count+1] = ahist['time'][i]
          tmp['ek1'][count+1] = ahist['ek1'][i]
          tmp['ek2'][count+1] = ahist['ek2'][i]
          tmp['ek3'][count+1] = ahist['ek3'][i]
          tmp['em1'][count+1] = ahist['em1'][i]
          tmp['em2'][count+1] = ahist['em2'][i]
          tmp['em3'][count+1] = ahist['em3'][i]
          tmp['etot'][count+1] = ahist['etot'][i]
          count +=1
      #print "finishing counting count= ",count
      tt  = tmp['time'][0:count]
      ek1 = tmp['ek1'][0:count]
      ek2 = tmp['ek2'][0:count]
      ek3 = tmp['ek3'][0:count]
      em1 = tmp['em1'][0:count]
      em2 = tmp['em2'][0:count]
      em3 = tmp['em3'][0:count]
      etot = tmp['etot'][0:count]
      return tt,ek1,ek2,ek3,em1,em2,em3,etot

    dtype = np.dtype([('time', 'f8'), ('ek1', 'f8'),('ek2','f8'),('ek3','f8'),('etot','f8'),('em1', 'f8'),('em2','f8'),('em3','f8')])
    ahist = np.loadtxt(fname, dtype=dtype, skiprows=2, usecols=(0,6,7,8,9,10,11,12))
    tt,ek1,ek2,ek3,em1,em2,em3,etot = time_trim(ahist)
    return tt,ek1,ek2,ek3,em1,em2,em3,etot

def read_hist1d_mdot(fname):
    def time_trim(ahist):
      count = 0; tmp = np.copy(ahist)
      for i in np.arange(1,len(ahist['time'])):
        if ahist['time'][i] > tmp['time'][count]:
          tmp['time'][count+1] = ahist['time'][i]
          tmp['mdot1'][count+1] = ahist['mdot1'][i]
          tmp['mdot2'][count+1] = ahist['mdot2'][i]
          count +=1
      #print "finishing counting count= ",count
      tt  = tmp['time'][0:count]
      mdot1 = tmp['mdot1'][0:count]
      mdot2 = tmp['mdot2'][0:count]
      return tt,mdot1,mdot2

    dtype = np.dtype([('time', 'f8'), ('mdot1', 'f8'),('mdot2','f8')])
    ahist = np.loadtxt(fname, dtype=dtype, skiprows=2, usecols=(0,10,11))
    tt,mdot1,mdot2 = time_trim(ahist)
    #tt,mdot1,mdot2 = ahist['time'],ahist['mdot1'],ahist['mdot2']
    return tt,mdot1,mdot2

def time_trim(ahist):
    count = 0; tmp = np.copy(ahist)
    for i in np.arange(1,len(ahist['time'])):
      if ahist['time'][i] > tmp['time'][count]:
        tmp['time'][count+1] = ahist['time'][i]
        tmp['em1'][count+1] = ahist['em1'][i]
        tmp['em2'][count+1] = ahist['em2'][i]
        tmp['em3'][count+1] = ahist['em3'][i]
        tmp['maxwell'][count+1] = ahist['maxwell'][i]
        tmp['reynolds'][count+1] = ahist['reynolds'][i]
        count +=1
    #print "finishing counting count= ",count
    tt  = tmp['time'][0:count]
    em1 = tmp['em1'][0:count]
    em2 = tmp['em2'][0:count]
    em3 = tmp['em3'][0:count]
    maxwell  = tmp['maxwell'][0:count]
    reynolds = tmp['reynolds'][0:count]
    return tt,em1,em2,em3,maxwell,reynolds





