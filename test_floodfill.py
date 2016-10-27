
import numpy as np
import athena4_read as ath
 
def flood_fill_tryout(image,x,y,nx,ny,value):
     edge = [(x, y)]
     image[x,y] = value
     while edge:
         newedge = []
         for (x, y) in edge:
             print "flood fill at (ix,iy) = (",x,",",y,")"
             for (s, t) in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
                 if ((s >=0 and s<nx and t>=0 and t<ny) and (image[s,t] >= value)):
                     image[s,t]=value
                     newedge.append((s, t))
                 else:
		     image[s,t]=0.0
                     continue
         edge = newedge
     return
#
# use a queue which consist the nearby 4 neightbor
#
def flood_fill_stack(image,x,y,color):
     if image[x,y] < color: 
       print "initial cell is below color = ",color
       return
     	  
     nx = image.shape[0]
     ny = image.shape[1]
     work = np.zeros([nx,ny]) # 0: unprocessed 1: processed
     queue  = [(x, y)] # flood fill stack
     jsheet = {} #{(x,y):image[x,y]} # current sheet dict loc:j^2 value
     count = 0
     while queue:
       (x,y) = queue[0] 
       queue.remove((x,y))
       #print 'process (x,y) = ',x,',',y
       if (image[x,y] >= color):
         if (work[x,y] != 1):
           work[x,y] = 1
           jsheet[(x,y)]=image[x,y]
         for (s, t) in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
           if (s>=0 and s<nx and t>=0 and t<ny and (work[s,t] !=1)):
             queue.append((s, t))
       else:
         work[x,y] = 1

       count += 1
       if count > 300:
	   print "exceed the hard limit count = ",count
           break
     return jsheet

#
# use a queue which search west/east neighbors continuously 
#
def flood_fill_westack(image,x,y,color):
     if image[x,y] < color: 
       print "initial cell is below color = ",color
       return
     	  
     nx = image.shape[0]
     ny = image.shape[1]
     work = np.zeros([nx,ny]) # 0: unprocessed 1: processed
     queue  = [(x, y)] # flood fill stack
     jsheet = {} # should be a list of dictionaries
                 # and each dict represent one current sheet
		 # {(x_0,y_0):image[x_0,y_0], (x_1,y_1):image[x_1,y_1],...} 
		 # which stores loc:j^2 pair
		 # for the moment(to test the stack), we use single dict
     count = 0
     while queue:
       for (x,y) in queue: 
	 queue.remove((x,y))
	 west = x
	 east = x
	 while(west-1>=0):
	   west -=1
	   if (work[west,y] >0 or image[west,y] < color): ##processed or below the color value
	     west +=1
	     break
         while(east+1<nx):
	   east +=1
	   if (work[east,y] >0 or image[east,y] < color): ##processed or below the color value
	     east -=1
	     break
         for i in range(west,east+1):
	   work[i,y] = 1
	   jsheet[(i,y)] = image[i,y]
	   for (s,t) in ((i,y-1),(i,y+1)):
             if (t>=0 and t<ny and (work[s,t] ==0) and image[s,t] >=color):
               queue.append((s,t))

       count += 1
       if count > 300:
	   print "exceed the hard limit count = ",count
           break
     return jsheet

#
# use a queue which search west/east neighbors continuously 
# now eliminate the use of work array by set image cell = 0 < color
# Called by a driver which loops over all local maxima in the image
# still 2d construction here.
#
def flood_fill_2d(image,x,y,color,jlist):
     if image[x,y] < color: 
       print "initial cell is below color = ",color
       return
     	  
     nx = image.shape[0]
     ny = image.shape[1]
     queue  = [(x, y)] # flood fill stack
     jsheet = {} # should be add to jlist to form list of current sheet(dictionaries)
                 # and each dict represent one current sheet
		 # {(x_0,y_0):image[x_0,y_0], (x_1,y_1):image[x_1,y_1],...} 
		 # which stores loc:j^2 pair
		 # for the moment(to test the stack), we use single dict
     #count = 0
     while queue:
       for (x,y) in queue: 
	 queue.remove((x,y))
	 west = x
	 east = x
	 while(west-1>=0):
	   west -=1
	   if (image[west,y] < color): ##processed or below the color value
	     west +=1
	     break
         while(east+1<nx):
	   east +=1
	   if (image[east,y] < color): ##processed or below the color value
	     east -=1
	     break
         for i in range(west,east+1):
	   if (image[i,y] >= color):
	     jsheet[(i,y)] = image[i,y]
	     image[i,y] = 0.0 # remove it from image
	     for (s,t) in ((i,y-1),(i,y+1)):
               if (t>=0 and t<ny and image[s,t] >=color):
                 queue.append((s,t))

       #count += 1
       #if count > 300:
       #    print "exceed the hard limit count = ",count
       #    break

     if jsheet: # found current sheet
       jlist.append(jsheet)

     return

#
# Now for 3d construction here.
# use a queue which search west/east neighbors continuously 
# now eliminate the use of work array by set image cell = 0 < color
# Called by a driver which loops over all local maxima in the image
#
def flood_fill_3d(image,x,y,z,color,jlist):
     if image[x,y,z] < color: 
       print "initial cell is below color = ",color
       return
     	  
     nx = image.shape[0]
     ny = image.shape[1]
     nz = image.shape[2]
     queue  = [(x, y, z)] # flood fill stack
     jsheet = {} # should be add to jlist to form list of current sheet(dictionaries)
                 # and each dict represent one current sheet
		 # {(x_0,y_0,z_0):image[x_0,y_0,z_0], (x_1,y_1,z_1):image[x_1,y_1,z_1],...} 
		 # which stores loc:j^2 pair
		 # for the moment(to test the stack), we use single dict
     #count = 0
     while queue:
       for (x,y,z) in queue: 
	 queue.remove((x,y,z))
	 west = x
	 east = x
	 north = y
	 south = y
	 while(west-1>=0):
	   west -=1
	   if (image[west,y] < color): ##processed or below the color value
	     west +=1
	     break
         while(east+1<nx):
	   east +=1
	   if (image[east,y] < color): ##processed or below the color value
	     east -=1
	     break
         for i in range(west,east+1):
	   if (image[i,y] >= color):
	     jsheet[(i,y)] = image[i,y]
	     image[i,y] = 0.0 # remove it from image
	     for (s,t) in ((i,y-1),(i,y+1)):
               if (t>=0 and t<ny and image[s,t] >=color):
                 queue.append((s,t))

       #count += 1
       #if count > 300:
       #    print "exceed the hard limit count = ",count
       #    break

     if jsheet: # found current sheet
       jlist.append(jsheet)

     return
if __name__=='__main__':
  
    t,x,y,z,data=ath.vtk('Unstra.0080.vtk')
    d2d = data['density'][0,:,:]
    flood_fill(d2d,0,0,256,64,1.0)
    print d2d
