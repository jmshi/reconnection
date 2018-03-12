#!/usr/bin/env python
import numpy as np
 
fdir = '/tigress/jiming/reconnect/athena/bin/'
targlist = ['x2y4z1r64pm1re4000','x2y4z1r128pm1re4000','x2y4z1r256pm1re4000'
ts,te = 30,50
for targ in targlist:
  python main.py targ ts te 
