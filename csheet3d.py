import athena_read as ath
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import cPickle as pickle
from my_colors import get_jhcolors
############################################
import matplotlib.gridspec as gridspec


import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def plot_jhist(dv,ncells,jmax,jmaxnorm,diss,dissnorm,size,orient,theta,ntime,vsheet=False,size2=None):
    nframe = 2
    nplot  = 3
    matplotlib.rcParams['figure.figsize'] = (12, 12.0*nframe/nplot)
    fig = plt.figure()
    fraction=0.046;pad=0.04
    gs = gridspec.GridSpec(2, 3)
    # calc dissipation rate per sheet for normalization or weighted average
    ncelltot= float(np.sum(np.array(ncells)))
    diss_sheet = np.array(ncells)*np.array(diss)
    disstot= float(np.sum(diss_sheet))
    for i in np.arange(6):
      if (i==0):
        hdata = np.array(size)[:,0]
        vmin=-2; vmax=0.5  #in log10
        nbins = 100
        xlab = r'$l/H$'
	# modify to plot the aspect ratio: xi/lambda
	hdata = np.array(size)[:,1]/np.array(size)[:,2]
        vmin=-0.5; vmax=3.0  #in log10
        nbins = 100
        xlab = r'$\xi/\lambda$'
	#
        print 'non-weighted averaged (l,xi,lambda,aspect) = (',np.average(np.array(size)[:,0]),' ,',\
                                           np.average(np.array(size)[:,1]),' ,',\
                                           np.average(np.array(size)[:,2]),',',\
					   np.average(np.array(size)[:,1]/np.array(size)[:,2]),')'
        print 'vol-weighted averaged (l,xi,lambda,aspect) = (',np.sum(np.array(size)[:,0]*np.array(ncells))/ncelltot,' ,',\
                                           np.sum(np.array(size)[:,1]*np.array(ncells))/ncelltot,' ,',\
                                           np.sum(np.array(size)[:,2]*np.array(ncells))/ncelltot,', ',\
					   np.sum(np.array(size)[:,1]/np.array(size)[:,2]*np.array(ncells))/ncelltot,')'      
        print 'dis-weighted averaged (l,xi,lambda,aspect) = (',np.sum(np.array(size)[:,0]*diss_sheet)/disstot,' ,',\
                                           np.sum(np.array(size)[:,1]*diss_sheet)/disstot,' ,',\
                                           np.sum(np.array(size)[:,2]*diss_sheet)/disstot,', ',\
					   np.sum(np.array(size)[:,1]/np.array(size)[:,2]*diss_sheet)/disstot,')'
        index,coef=-3.2,1e-2
      if (i==1):
        hdata = np.array(size)[:,0]
        vmin=-3;vmax=1
        nbins = 100
        xlab1 = r'$\xi$'
        xlab =r'$l$'
        xlab2 = r'$\lambda$'
        hdata1 = np.array(size)[:,1]
        hdata2 = np.array(size)[:,2]
        index,coef=-3.0,1e-2
	if size2:
	  print "reach size2 loop"
	  hdata3 = np.array(size2)[:,2]
      if (i==2):
        hdata = (theta)
        vmin = np.log10(np.min(hdata))
        vmax = np.log10(np.max(hdata))
        nbins = 100
	if vsheet:
          xlab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
          print 'average theta/E_nu = ', np.average(hdata)
	else:
          xlab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
          print 'average theta/E_eta = ', np.average(hdata)
        index,coef=-2.0,1e-5
      if (i==3):
        hdata = (np.array(jmaxnorm))
        vmin = np.log10(2) #np.log10(np.min(hdata))
        vmax = np.log10(20) #np.log10(np.max(hdata))
        nbins=100
	if vsheet:
          xlab = r'$\omega_{max}/\omega_{rms}$'
          print 'average w_max/w_rms = ',np.average(hdata)
	else:
          xlab = r'$j_{max}/j_{rms}$'
          print 'average j_max/j_rms = ',np.average(hdata)
      if (i==4):
        hdata = np.array(dissnorm)
        vmin = np.log10(5) #np.log10(np.min(hdata))
        vmax = np.log10(50) #np.log10(np.max(hdata))
        nbins = 100
	if vsheet:
          xlab = r'$\langle\varepsilon_{\nu,sh}\rangle_i /\langle\varepsilon_{\nu}\rangle$'
          print 'average epsilon = ', np.average(hdata)
	else:
          xlab = r'$\langle\varepsilon_{\eta,sh}\rangle_i /\langle\varepsilon_{\eta}\rangle$'
          print 'average epsilon = ', np.average(hdata)
      if (i==5):
        hdata = np.array(ncells)*dv #/((x[1]-x[0])**3)
        #vmin = -3; vmax = 0
        vmin = np.log10(np.min(hdata))
        vmax = np.log10(np.max(hdata))
        nbins = 100
        xlab = r'$v_{sh}/V$'
        print 'average volume fraction per sheet = ',np.average(hdata)/8.0
        print 'dis-weighted average volume fraction per sheet = ',np.average(hdata*diss_sheet)\
                                                             /8.0/np.average(diss_sheet)
        print 'average volume fraction for all sheets = ',np.sum(hdata)/float(ntime)/8.0
        print 'time averaged number of current sheets = ', len(ncells)/float(ntime)
        #print 'dis-weighted average volume fraction for all sheets = ',np.sum(hdata*diss_sheet)\
        #                                                        /float(ntime)/8.0/(disstot)

      hist, bins = np.histogram(hdata, bins=np.logspace(vmin,vmax,nbins),density=1)
      center = 10**((np.log10(bins[:-1])+np.log10(bins[1:]))/2)
      if (i==1):
         hist1, bins1 = np.histogram(hdata1, bins=np.logspace(vmin,vmax,nbins),density=1)
         hist2, bins2 = np.histogram(hdata2, bins=np.logspace(vmin,vmax,nbins),density=1)
	 if size2:
           hist3, bins3 = np.histogram(hdata3, bins=np.logspace(vmin,vmax,nbins),density=1)
      
      fig.add_subplot(gs[i])
      if (i< 3): 
        plt.plot(center,hist,'k-',label=xlab)
        #n, bins, patches = plt.hist(hdata, bins=np.logspace(vmin,vmax,nbins), \
        #                            normed=1, facecolor='green', alpha=0.5,label=xlab)
        if i==1:
          plt.plot(center,hist1,'b-',label=xlab1)
          plt.plot(center,hist2,'g-',label=xlab2)
          #n, bins, patches = plt.hist(hdata1, bins=np.logspace(vmin,vmax,nbins), \
          #                          normed=1, facecolor='blue', alpha=0.35,label=xlab1)
	  if size2:
            plt.plot(center,hist3,'g--')

          plt.legend(loc=1,fontsize=20)
          xlab = r'$(l,\xi,\lambda)/H$'
        #plt.text(vmin,7000,'total counts: '+str(len(hdata)),size=20,color='g')
        plt.plot(center,coef*center**index,'k:')
      else:
        plt.plot(center,hist,'k-')
        #n, bins, patches = plt.hist(hdata, bins=np.logspace(vmin,vmax,nbins), normed=1, facecolor='blue', alpha=0.5)

      plt.xlabel(xlab,size=20)
      plt.ylabel('probability density',size=20)
      plt.xlim([10**vmin,10**vmax])
      plt.gca().set_xscale("log")
      plt.gca().set_yscale("log",nonposy='clip')
      #if(i==5):
      #  plt.xlim([1e4,10**4.3])

    plt.tight_layout()

## only use the following for the paper
def plot_jhist_2panel(dv,ncells,jmax,jmaxnorm,diss,dissnorm,size,orient,theta,ntime,vsheet=False,size2=None):
    nframe = 2
    nplot  = 1
    matplotlib.rcParams['figure.figsize'] = (8.0,4.0)
    fig = plt.figure()
    fraction=0.046;pad=0.04
    gs = gridspec.GridSpec(nplot, nframe)
    # calc dissipation rate per sheet for normalization or weighted average
    ncelltot= float(np.sum(np.array(ncells)))
    diss_sheet = np.array(ncells)*np.array(diss)
    disstot= float(np.sum(diss_sheet))
    for i in np.arange(2):
      if (i==0):
        hdata = np.array(size)[:,0]
        vmin=-3;vmax=1
        nbins = 100
        xlab1 = r'$\xi$'
        xlab =r'$l$'
        xlab2 = r'$\lambda_2$'
        hdata1 = np.array(size)[:,1]
        hdata2 = np.array(size)[:,2]
        index,coef=-3.0,1e-2
	if size2:
	  print "reach size2 loop"
	  hdata3 = np.array(size2)[:,2]
	  xlab3 = r'$\lambda_1$'
      if (i==1):
        hdata = (theta)
        vmin = np.log10(np.min(hdata))
        vmax = np.log10(np.max(hdata))
        nbins = 100
	if vsheet:
          xlab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
          print 'average theta/E_nu = ', np.average(hdata)
	else:
          xlab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
          print 'average theta/E_eta = ', np.average(hdata)
        index,coef=-2.0,1e-5

      hist, bins = np.histogram(hdata, bins=np.logspace(vmin,vmax,nbins),density=1)
      center = 10**((np.log10(bins[:-1])+np.log10(bins[1:]))/2)
      if (i==0):
         hist1, bins1 = np.histogram(hdata1, bins=np.logspace(vmin,vmax,nbins),density=1)
         hist2, bins2 = np.histogram(hdata2, bins=np.logspace(vmin,vmax,nbins),density=1)
	 if size2:
           hist3, bins3 = np.histogram(hdata3, bins=np.logspace(vmin,vmax,nbins),density=1)
      
      fig.add_subplot(gs[i])
      if (i< 3): 
        plt.plot(center,hist,'k-',label=xlab)
        #n, bins, patches = plt.hist(hdata, bins=np.logspace(vmin,vmax,nbins), \
        #                            normed=1, facecolor='green', alpha=0.5,label=xlab)
        if i==0:
          plt.plot(center,hist1,'b-',label=xlab1)
	  if size2:
            plt.plot(center,hist3,'r-',label=xlab3)
          plt.plot(center,hist2,'g-',label=xlab2)
          #n, bins, patches = plt.hist(hdata1, bins=np.logspace(vmin,vmax,nbins), \
          #                          normed=1, facecolor='blue', alpha=0.35,label=xlab1)

          plt.legend(loc=1,fontsize=20,handlelength=1.5,handletextpad=0.05)
          xlab = r'$(l,\xi,\lambda)/H$'
        #plt.text(vmin,7000,'total counts: '+str(len(hdata)),size=20,color='g')
        plt.plot(center,coef*center**index,'k:')
      else:
        plt.plot(center,hist,'k-')
        #n, bins, patches = plt.hist(hdata, bins=np.logspace(vmin,vmax,nbins), normed=1, facecolor='blue', alpha=0.5)

      plt.xlabel(xlab,fontsize=20)
      plt.ylabel('probability density',fontsize=20)
      plt.xlim([10**vmin,10**vmax])
      plt.gca().set_xscale("log")
      plt.gca().set_yscale("log",nonposy='clip')
      #if(i==5):
      #  plt.xlim([1e4,10**4.3])

    plt.tight_layout()

def plot_jcorr(ncells,jmax,jmaxnorm,diss,dissnorm,size,orient,theta,disscut=4e-5,ntime=21,project=np.array([1,0,0]),cmap='jet',vsheet=False):
    """
    plot the correlation between quantities within current sheet
    """
    nframe = 3
    nplot  = 3
    enlarge = 2
    matplotlib.rcParams['figure.figsize'] = (10*enlarge, enlarge*8.0*nframe/nplot)
    fig = plt.figure()
    fraction=0.046;pad=0.04
    gs = gridspec.GridSpec(nframe, nplot)
    mask = (theta > disscut)
    # calc dissipation rate per sheet for normalization or weighted average
    ncelltot= float(np.sum(np.array(ncells)))
    diss_sheet = np.array(ncells)*np.array(diss)
    disstot= float(np.sum(diss_sheet))
    
    for i in np.arange(9):
    #######################################
    # (1) orientation
    #######################################
      if (i==0):
        vec2d = np.array(orient)[mask,0]#np.cross(np.array(orient)[mask,0],project)
        hdata_x = vec2d[:,2]
        hdata_y = vec2d[:,1]
        xlab = r'$\hat{l}_x$'
        ylab = r'$\hat{l}_y$'
        xmin,xmax = -1,1
        ymin,ymax = -1,1
	vec2d = np.array(orient)[:,0]
	vecnorm = np.average(vec2d,axis=0)
	vecnorm /= np.linalg.norm(vecnorm)
        print 'non-weighted averaged (\hat{l}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(np.array(ncells)*vec2d.T,axis=1)/ncelltot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'vol-weighted averaged (\hat{l}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(diss_sheet*vec2d.T,axis=1)/disstot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'dis-weighted averaged (\hat{l}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
      if (i==1):
        vec2d = np.array(orient)[mask,1]
        hdata_x = vec2d[:,2]
        hdata_y = vec2d[:,1]
        xlab = r'$\hat{\xi}_x$'
        ylab = r'$\hat{\xi}_y$'
        xmin,xmax = -1,1
        ymin,ymax = -1,1
	vec2d = np.array(orient)[:,1]
	vecnorm = np.average(vec2d,axis=0)
	vecnorm /= np.linalg.norm(vecnorm)
        print 'non-weighted averaged (\hat{xi}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(np.array(ncells)*vec2d.T,axis=1)/ncelltot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'vol-weighted averaged (\hat{xi}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(diss_sheet*vec2d.T,axis=1)/disstot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'dis-weighted averaged (\hat{xi}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
      if (i==2):
        vec2d = np.array(orient)[mask,2]
        hdata_x = vec2d[:,2]
        hdata_y = vec2d[:,1]
        xlab = r'$\hat{\lambda}_x$'
        ylab = r'$\hat{\lambda}_y$'
        xmin,xmax = -1,1
        ymin,ymax = -1,1
	vec2d = np.array(orient)[:,2]
	vecnorm = np.average(vec2d,axis=0)
	vecnorm /= np.linalg.norm(vecnorm)
        print 'non-weighted averaged (\hat{\lambda}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(np.array(ncells)*vec2d.T,axis=1)/ncelltot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'vol-weighted averaged (\hat{\lambda}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(diss_sheet*vec2d.T,axis=1)/disstot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'dis-weighted averaged (\hat{\lambda}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
    #######################################
    # (2) 2d histogram of sheet dimensions
    #######################################
      if (i==3):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(size)[:,1]
        xlab = r'$l$'
        ylab = r'$\xi$'
        xmin,xmax = -2.7,1
        ymin,ymax = -2.7,0
        nbin = 100
      if (i==4):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(size)[:,2]
        xlab = r'$l$'
        ylab = r'$\lambda$'
        xmin,xmax = -2.7,1
        ymin,ymax = -2.7,0
        nbin = 100
      if (i==5):
        hdata_x = np.array(size)[:,1]
        hdata_y = np.array(size)[:,2]
        xlab = r'$\xi$'
        ylab = r'$\lambda$'
        xmin,xmax = -2.7,1
        ymin,ymax = -2.7,0
        nbin = 100
	# modify to plot the corr btw xi/lambda and lambda
        hdata_x = np.array(size)[:,1]
	hdata_y = np.array(size)[:,1]/np.array(size)[:,2]
        xlab = r'$\xi$'
        ylab = r'$\xi/\lambda$'
        xmin,xmax = -3,1.
        ymin,ymax = 0.0, 3. #-0.5,3.0
        nbin = 100
    #######################################
    # (3) dissipation below
    #######################################
      if (i==6):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(theta)
        xlab = r'$l$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -2.7,1
        nbin = 100
      if (i==7):
        hdata_x = np.array(size)[:,1]
        hdata_y = np.array(theta)
        xlab = r'$\xi$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -2.7,1
        nbin = 100
	# modify to plot corr btwn theta and aspect ratio
	hdata_x = np.array(size)[:,1]/np.array(size)[:,2]
        hdata_y = np.array(theta)
        xlab = r'$\xi/\lambda$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -0.5,2.5
        nbin = 100
      if (i==8):
        hdata_x = np.array(size)[:,2]
        hdata_y = np.array(theta)
        xlab = r'$\lambda$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -2.7,1
        nbin = 100
        
      if i !=-1:
        fig.add_subplot(gs[i],aspect='equal')  
      else:
        fig.add_subplot(gs[i])  
     
      if i<3:
        plt.scatter(hdata_x,hdata_y,s=0.02, marker = '.' )
        plt.ylim([ymin,ymax])
        plt.xlim([xmin,xmax])
        plt.xlabel(xlab,size=20)
        plt.ylabel(ylab,size=20)
        if i==0: # add tilt angle measured in Zhdankin paper -17.5degree
            tilt = -17.5/180.*np.pi
            rx,ry = 2.*np.sin(tilt),2.*np.cos(tilt)
            plt.plot([rx,-rx],[ry,-ry],'k:')
      elif i<6:
	#plt.scatter(hdata_x,hdata_y,s=0.01,marker='.')
	#hst are binned in the dissipation values (normalized with total dissipation)
        hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=(np.logspace(ymin,ymax,nbin),\
                                                         np.logspace(xmin,xmax,nbin)),\
                                                         weights=theta) #diss_sheet/disstot)
        xx,yy = np.meshgrid(xe,ye)
        loghst = np.log10(hst+1e-12)
        plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-5,vmax=-1,cmap=cmap)
			#vmax=np.max(loghst),cmap=cmap)

        #print 'min/max of hst = ',np.min(loghst),np.max(loghst)
	#print 'sum(hst) = ',np.sum(hst)
        #print 'min/max of xe = ',np.min(xe),np.max(xe)
        #print 'min/max of ye = ',np.min(ye),np.max(ye)
	
        plt.colorbar(pad=pad,fraction=fraction*0.8)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim([10**xmin,10**xmax])
        plt.ylim([10**ymin,10**ymax])
        plt.xlabel(xlab,size=20)
        plt.ylabel(ylab,size=20)
	if vsheet:
	  plt.title(r'$\theta_{\nu}/\mathcal{E}_{\nu}$',fontsize=15)
	else: 
	  plt.title(r'$\theta_{\eta}/\mathcal{E}_{\eta}$',fontsize=15)
        if i==3:
          index,coef=1.0,0.4
          plt.plot(xe,coef*xe**index,'k:')
      #plt.title(r'$\mathrm{Histogram\ of\ d_i:}$')
      else:
        hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=[np.logspace(ymin,ymax,nbin),\
                                                         np.logspace(xmin,xmax,nbin)],\
                                                         normed=True)
        xx,yy = np.meshgrid(xe,ye)
        loghst = np.log10(hst+1e-12)
        plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-3,vmax=9)
        #print 'min/max of hst = ',np.min(loghst),np.max(loghst)
        plt.colorbar(pad=pad,fraction=fraction)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim([10**xmin,10**xmax])
        plt.ylim([10**ymin,10**ymax])
        plt.xlabel(xlab,size=20)
        plt.ylabel(ylab,size=20)
	plt.title(r'$\log$'+'(PDF)',fontsize=15)
        if i==6:
	    index,coef=2.0,1.5e-3
            plt.plot(ye,coef*ye**index,'k:')
        
    plt.tight_layout()
    

       # linear regression
      #n_abovezero = np.min([len(hdata_x[hdata_x > 0]),len(hdata_y[hdata_y > 0])])
      #clean_data = np.log10(np.array(zip(hdata_x[:n_abovezero],hdata_y[:n_abovezero])))
      #print clean_data
      #coeff = np.polyfit(clean_data[:,0], clean_data[:,1], 1)
      #print coeff
      #break
      #yfit = 10**(coeff[0]*clean_data[:,0]+coeff[1])
      #ax.plot(clean_data[:,0],yfit,'r-')

# use the following to generate the figure used in the paper
def plot_jcorr_4panel(ncells,jmax,jmaxnorm,diss,dissnorm,size,orient,theta,disscut=4e-5,ntime=21,project=np.array([1,0,0]),cmap='jet',vsheet=False):
    """
    plot the correlation between quantities within current sheet
    """
    nframe = 2
    nplot  = 2
    enlarge = 1
    matplotlib.rcParams['figure.figsize'] = (8*enlarge, enlarge*7.0*nframe/nplot)
    fig = plt.figure()
    fraction=0.046;pad=0.04
    gs = gridspec.GridSpec(nframe, nplot)
    mask = (theta > disscut)
    # calc dissipation rate per sheet for normalization or weighted average
    ncelltot= float(np.sum(np.array(ncells)))
    diss_sheet = np.array(ncells)*np.array(diss)
    disstot= float(np.sum(diss_sheet))
    
    for i in np.arange(4):
    #######################################
    # (1) orientation
    # (1) replace it with aspect ratio vs. xi
    #######################################
      if (i==2):
        #vec2d = np.array(orient)[mask,0]#np.cross(np.array(orient)[mask,0],project)
        #hdata_x = vec2d[:,2]
        #hdata_y = vec2d[:,1]
        #xlab = r'$\hat{l}_x$'
        #ylab = r'$\hat{l}_y$'
        #xmin,xmax = -1,1
        #ymin,ymax = -1,1
	# modify to plot the corr btw xi/lambda and lambda
        hdata_x = np.array(size)[:,1]
	hdata_y = np.array(size)[:,1]/np.array(size)[:,2]
        xlab = r'$\xi$'
        ylab = r'$\xi/\lambda$'
        xmin,xmax = -3,1.
        ymin,ymax = 0.0, 3. #-0.5,3.0
        nbin = 100
    #######################################
    # (2) 2d histogram of sheet dimensions
    #######################################
      if (i==0):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(size)[:,1]
        xlab = r'$l$'
        ylab = r'$\xi$'
        xmin,xmax = -3.0,1.0 #2.7,1
        ymin,ymax = -3.0,1.0 #-2.7,0
        nbin = 100
      if (i==1):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(size)[:,2]
        xlab = r'$l$'
        ylab = r'$\lambda_2$'
        xmin,xmax = -3.0,1.0 #-2.7,1
        ymin,ymax = -3.0,-1.8 #-2.7,0
        nbin = 100
    #######################################
    # (3) dissipation below
    #######################################
      if (i==3):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(theta)
        xlab = r'$l$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -2.7,1
        nbin = 100
        
      #if i==2: 
      #	fig.add_subplot(gs[i] ,aspect='equal')  
      #else:
      fig.add_subplot(gs[i])  

     
      #if i==2:
      #  plt.scatter(hdata_x,hdata_y,s=0.03, marker = '.' )
      #  plt.ylim([ymin,ymax])
      #  plt.xlim([xmin,xmax])
      #  plt.xlabel(xlab,size=20)
      #  plt.ylabel(ylab,size=20)
      #  # add tilt angle measured in Zhdankin paper -17.5degree
      #  tilt = -17.5/180.*np.pi
      #  rx,ry = 2.*np.sin(tilt),2.*np.cos(tilt)
      #  plt.plot([rx,-rx],[ry,-ry],'k:')
      #if i<2:
      if i<3:
	#plt.scatter(hdata_x,hdata_y,s=0.01,marker='.')
	#hst are binned in the dissipation values (normalized with total dissipation)
        hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=(np.logspace(ymin,ymax,nbin),\
                                                         np.logspace(xmin,xmax,nbin)),\
                                                         weights=theta) #diss_sheet/disstot)
        xx,yy = np.meshgrid(xe,ye)
        loghst = np.log10(hst+1e-12)
        plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-5,vmax=-1,cmap=cmap)
			#vmax=np.max(loghst),cmap=cmap)

        #print 'min/max of hst = ',np.min(loghst),np.max(loghst)
	#print 'sum(hst) = ',np.sum(hst)
        #print 'min/max of xe = ',np.min(xe),np.max(xe)
        #print 'min/max of ye = ',np.min(ye),np.max(ye)
	
        plt.colorbar(pad=pad,fraction=fraction)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim([10**xmin,10**xmax])
        plt.ylim([10**ymin,10**ymax])
        plt.xlabel(xlab,size=20)
        plt.ylabel(ylab,size=20)
	if vsheet:
	  plt.title(r'$\theta_{\nu}/\mathcal{E}_{\nu}$',fontsize=15)
	else: 
	  plt.title(r'$\theta_{\eta}/\mathcal{E}_{\eta}$',fontsize=15)
        if i==0:
          index,coef=1.0,0.4
          plt.plot(xe,coef*xe**index,'k:')
        if i==1:
          index,coef=1./3.,0.008
	  xe = 10**(np.linspace(-3.0,0.4,100))
          plt.plot(xe,coef*xe**index,'k:')
          index,coef=0.,0.0065
	  xe = 10**(np.linspace(-1.8,1.0,50))
          plt.plot(xe,coef*xe**index,'k:')
        if i==2:
          index,coef=2./3.,85.
	  xe = 10**(np.linspace(-3.0,0.4,100))
          plt.plot(xe,coef*xe**index,'k:')
          index,coef=1.,160.
	  xe = 10**(np.linspace(-1.8,1.0,50))
          plt.plot(xe,coef*xe**index,'k:')
      #plt.title(r'$\mathrm{Histogram\ of\ d_i:}$')
      if i==3:
        hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=[np.logspace(ymin,ymax,nbin),\
                                                         np.logspace(xmin,xmax,nbin)],\
                                                         normed=True)
        xx,yy = np.meshgrid(xe,ye)
        loghst = np.log10(hst+1e-12)
        plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-3,vmax=9,cmap=cmap)
        #print 'min/max of hst = ',np.min(loghst),np.max(loghst)
        plt.colorbar(pad=pad,fraction=fraction)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim([10**xmin,10**xmax])
        plt.ylim([10**ymin,10**ymax])
        plt.xlabel(xlab,size=20)
        plt.ylabel(ylab,size=20)
	plt.title(r'$\log$'+'(PDF)',fontsize=15)
	index,coef=2.25,2.2e-3
        plt.plot(ye,coef*ye**index,'k:')
        
    plt.tight_layout()
    

# use the following to generate the figure used in the paper
def plot_jcorr_null(size,theta,size1,theta1,cmap1='hot_r',cmap='bone_r'):
    """
    compare sheets with null and sheets w/o null
    draw nonull sheets first then overlay with null sheets
    then draw colorbar of these two..
    """
    matplotlib.rcParams['figure.figsize'] = (5.5,5)
    fig = plt.figure()
    fraction=0.046;pad=0.04
    
    hdata_x = np.array(size)[:,0]
    hdata_y = np.array(size)[:,1]
    xlab = r'$l$'
    ylab = r'$\xi$'
    xmin,xmax = -3.0,1.0 #2.7,1
    ymin,ymax = -3.0,1.0 #-2.7,0
    nbin = 100
    # 1) plot all sheets
    cs0 = plt.scatter(hdata_x,hdata_y,s=3,marker='s',color='b')
#hst are binned in the dissipation values (normalized with total dissipation)
    #hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=(np.logspace(ymin,ymax,nbin),\
    #                                                     np.logspace(xmin,xmax,nbin)),\
    #                                                     weights=theta) #diss_sheet/disstot)
    #xx,yy = np.meshgrid(xe,ye)
    #loghst = np.log10(hst+1e-12)
    #cs1 = plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-8,vmax=-4,cmap=cmap)
	
    #cbar1 = fig.colorbar(cs1,pad=pad,fraction=fraction,orientation='horizontal',ticks=[-6,-5,-4])
    ##plt.colorbar(cs1,pad=pad,fraction=fraction,orientation='horizontal',location='top')
    
    #2) now plot the null sheets
    hdata_x = np.array(size1)[:,0]
    hdata_y = np.array(size1)[:,1]
    hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=(np.logspace(ymin,ymax,nbin),\
                                                         np.logspace(xmin,xmax,nbin)),\
                                                         weights=theta1) #diss_sheet/disstot)
    xx,yy = np.meshgrid(xe,ye)
    loghst = np.log10(hst+1e-12)
    cs2 = plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-5,vmax=-1,alpha=0.8,cmap=cmap1)
	
    #plt.colorbar(cs2,pad=pad,fraction=fraction,orientation='vertical')
    cbar2 = fig.colorbar(cs2,pad=pad,fraction=fraction,orientation='vertical',ticks=[-5,-4,-3,-2,-1])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([10**xmin,10**xmax])
    plt.ylim([10**ymin,10**ymax])
    plt.xlabel(xlab,size=20)
    plt.ylabel(ylab,size=20)
    plt.title(r'$\theta_{\eta}/\mathcal{E}_{\eta}$',fontsize=20)
    index,coef=1.0,0.4
    plt.plot(xe,coef*xe**index,'k:')


    plt.axes().set_aspect('equal')

# script to get big Epsilon normalized with total dissipation
def get_path(targname='x2y4z1r64pmre4000'):
  direname='/tigress/jiming/reconnect/athena/bin/'
  targname=targname+'/'
  return direname+targname

def get_normtheta(targ,ts,te,tstride,eta=2.5e-4,vsheet=False,thick=False,null=False,wth=1):
  direname = get_path(targ)
  # 1. load the volume averaged data
  if vsheet:
    fname = direname+'Unstra.out2.wth.tab'
  else:
    fname = direname+'Unstra.out2.jth.tab'

  dtype = np.dtype([('time', 'f8'), ('jth2', 'f8'),('jsh2','f8'),('jrms2','f8')])
  ahist = np.loadtxt(fname, dtype=dtype, usecols=(0,1,2,3))
  tt = ahist['time']
  jrms2 = ahist['jrms2']
  jsh2  = ahist['jsh2']
  # pick out ts to te 
  istart,iend = 0,0
  if wth==1:
    for i in np.arange(len(tt)):
      if ts == int(tt[i]):
        istart = i
	break
    for i in np.arange(len(tt)):
      if (te-1) == int(tt[i]):
        iend = i
	break
  else: #wth=2 where rho\omega^2 is calc instead omega alone
    for i in np.arange(len(tt)):
      if ts == int(tt[i]):
        istart = i
    for i in np.arange(len(tt)):
      if (te-1) == int(tt[i]):
        iend = i
  iend +=1
  if iend==istart:
    iend+=1
  tt = tt[istart:iend]
  jrms2 = jrms2[istart:iend]
  jsh2  = jsh2[istart:iend]
  jrms = np.sqrt(jrms2)
  if vsheet:
    print '<epsilon_{nu,sh}>/<epsilon_{nu}> = ',np.average(jsh2/jrms2)
  else:
    print '<epsilon_{eta,sh}>/<epsilon_{eta}> = ',np.average(jsh2/jrms2)
  # 2. now load current sheets for each frame
  # for each sheet we calc total dissipation
  # and normalize it with  jrms^2
  if vsheet:
    if wth == 1:
      basename='vprop.'
      appdname='.p'
    else:
      basename='vprop2.'
      appdname='.p'
  else:
    basename='jprop.'
    appdname='.p'
    if thick:
     basename='jprop_thick.'
    if null:
     basename='jprop_null.'
  first_time = True
  for i in np.arange(ts,te,tstride):
    fname = direname+basename+str(i).zfill(5)+appdname
    ncells,jmax,diss,size,orient = pickle.load( open(fname, "rb" ) )
    if i != int(tt[i-ts]):
        print 'at t = ',tt[i-ts]
        return 0
    else:
        
        # sum(eta jsh^2)/ sum(eta jrms2) dissipation rate per sheet normalized
        # with rms dissipation.equivalently it's theta normalized with 
        # instantaneous total dissipation (although the number of grid cells not included here)
        diss =np.array(diss)/eta
        ncells = np.array(ncells)
        if first_time:
          dissnorm = diss/jrms2[i-ts]
          theta = diss*ncells/jrms2[i-ts]
          jmaxnorm = np.array(jmax)/jrms[i-ts]
          first_time = False
        else:
          dissnorm = np.hstack((dissnorm,diss/jrms2[i-ts]))
          theta = np.hstack((theta,diss*ncells/jrms2[i-ts]))
          jmaxnorm = np.hstack((jmaxnorm,np.array(jmax)/jrms[i-ts]))
#         print 'at t = ',tt[i-ts]
#         print 'jsh2 from pickled file: ',np.sum(diss*ncells)/float(np.sum(ncells))
#         print 'jsh2 from tabulat file:',jsh2[i-ts]
#         print 'jrms2 from tabulat file:',jrms2[i-ts]
  return dissnorm,theta,jmaxnorm     
#   nbin = 100; nt = te-ts+1
#   jth2 = np.zeros(nbin)
#   cdfj2 = np.zeros(nbin)
#   cdfvol = np.zeros(nbin)
# theta1 = get_normtheta('x2y4z1r64pm1re4000',50,51,1)
# theta2 = get_normtheta('x2y4z1r64pm1re4000',50,101,1)


def get_curl(fname='Unstra.out2.00008.athdf'):
  """load 3d bfield and calc the current density"""
  # ---
  def curl(vx,vy,vz,dx,dy,dz):
    [dzvx,dyvx,dxvx] = np.gradient(vx)
    [dzvy,dyvy,dxvy] = np.gradient(vy)
    [dzvz,dyvz,dxvz] = np.gradient(vz)
    j2 = (dyvz/dy-dzvy/dz)**2
    j2 += (dzvx/dz-dxvz/dx)**2
    j2 += (dxvy/dx-dyvx/dy)**2
    return j2
  # ---
  #data=ath.athdf(fname,quantities=['B1','B2','B3'])
  time,data=ath.athdf(fname,quantities=['Bcc1'])
  bx = data['Bcc1']
  time,data=ath.athdf(fname,quantities=['Bcc2'])
  by = data['Bcc2']
  time,data=ath.athdf(fname,quantities=['Bcc3'])
  bz = data['Bcc3']
  x  = data['x1f'];y  = data['x2f']; z  = data['x3f']
  dx = dz = x[1]-x[0]; dy = y[1]-y[0]
  #j2 = curl(bx,by,bz,dx,dy,dz)
  j2 = curl(bx,by,bz,dx,dy,dz)

  time,data=ath.athdf(fname,quantities=['vel1'])
  bx = data['vel1']
  time,data=ath.athdf(fname,quantities=['vel2'])
  by = data['vel2']
  time,data=ath.athdf(fname,quantities=['vel3'])
  bz = data['vel3']
  time,data=ath.athdf(fname,quantities=['rho'])
  d  = data['rho']
  #w2 = curl(bx,by,bz,dx,dy,dz)*data['rho']
  w2 = curl(bx,by,bz,dx,dy,dz)*d

  return time,x,y,z,j2,w2


def plot_3slice(input_data,xslice=512,yslice=504,zslice=256,dmin=0,dmax=100,cmap='OrRd',figsize_x=6,figsize_y=8,label=None):
    tmp = input_data
    matplotlib.rcParams['figure.figsize'] = (figsize_x, figsize_y)
    matplotlib.rcParams['image.cmap'] = cmap
    gs = gridspec.GridSpec(2, 2,width_ratios=[2, 1],height_ratios=[4, 1])

    ax1 = plt.subplot(gs[0])
    fig1 = plt.imshow(tmp[zslice,:,:],origin='lower',extent =[-1,1,-2,2],vmin=dmin,vmax=dmax)
    plt.ylabel('y',fontsize=15)
    ax2 = plt.subplot(gs[1])
    fig2 = plt.imshow(np.transpose(tmp[:,:,xslice]),origin='lower',extent=[-0.5,0.5,-2,2],vmin=dmin,vmax=dmax)
    ax3 = plt.subplot(gs[2])
    fig3 = plt.imshow(tmp[:,yslice,:],origin='lower',extent=[-1,1,-0.5,0.5],vmin=dmin,vmax=dmax)
    plt.ylabel('z',fontsize=15)
    plt.xlabel('x',fontsize=15)
    ax4 = plt.subplot(gs[3])
    ax4.axis('off')

    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("left", size="8%", pad=0.15)

    ax1.tick_params(labelbottom='off') 
    ax2.tick_params(labelleft='off')

    #ax2.get_yaxis().set_ticks([])
    plt.colorbar(cax=cax)
    if label != None:
      ax3.text(1.8,0.08,label,fontsize=20,rotation=90)

    plt.tight_layout()



def show_contour(x,y,input_data,dmin=0,dmax=100,cmap='OrRd',xspan=(0,0),yspan=(0,0),figsize_x=6,figsize_y=8,xlabel='x',ylabel='y'):
    tmp = input_data
    matplotlib.rcParams['figure.figsize'] = (figsize_x, figsize_y)
    matplotlib.rcParams['image.cmap'] = cmap
    #gs = gridspec.GridSpec(2, 2,width_ratios=[2, 1],height_ratios=[4, 1])
    if xspan[0] == 0 and xspan[1] == 0:
        xmin = -1; xmax=1
    else:
        xmin = xspan[0]; xmax=xspan[1]
    if yspan[0] == 0 and yspan[1] == 0:
        ymin = -2; ymax=2
    else:
        ymin = yspan[0]; ymax=yspan[1]
    
    # reduce the img to required size
    ids=np.abs(x-xmin).argmin()
    ide=np.abs(x-xmax).argmin()
    jds=np.abs(y-ymin).argmin()
    jde=np.abs(y-ymax).argmin()
    ndim = np.shape(tmp)
    jde = jde if jde < ndim[0] else ndim[0]-1
    ide = ide if ide < ndim[1] else ndim[1]-1
    tmp = tmp[jds:jde,ids:ide]

    fig1 = plt.imshow(tmp,origin='lower',extent =[xmin,xmax,ymin,ymax],vmin=dmin,vmax=dmax)
    plt.ylabel(ylabel,fontsize=15)
    plt.xlabel(xlabel,fontsize=15)
    plt.colorbar(pad=0.05,fraction=0.05)




# plot 8-panels in x-z plane (jmag,bx,by,bz,wmag/rho,vx,dvy,vz)
def plot_8panel(t,x,y,j2u,w2u,b1,b2,b3,v1,v2,v3,d,xspan=(-1,1),yspan=(-0.5,0.5),yslice=502,step=10,nslice=1,isrho=False):
    matplotlib.rcParams['figure.figsize'] = (15,40)
    x1b = (x + 0.5*(x[1]-x[0]))[:-1]
    x3b = (y + 0.5*(y[1]-y[0]))[:-1]
    nx,nz = len(x1b),len(x3b)
    count = 1
    jmin,jmax = 0,5
    bmin,bmax = -2,2
    vmin,vmax = -1.5,1.5
    rmin,rmax = 0,2
    #xspan,yspan = (0.3,0.7),(-0.5,-0.1)
    ncol,nrow = 4,22

    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height

    cmap = 'seismic'
    #yslice= 502
    #step = 10

    wu = np.sqrt(w2u)
    wu /= np.average(wu)
    ju = np.sqrt(j2u)
    ju /= np.average(ju)

    for i in (np.arange(nslice)-nslice/2):
      ind = yslice+step*i
      for nvar in np.arange(8): #// j2,vx,bx,by
        cmap = 'seismic'
        if nvar == 0:
          var = ju[:,ind,:] #ysli(np.sqrt(j2u)/np.average(np.sqrt(j2u)))[:,yslice+step*i,:] 
          cmap = 'OrRd'
          dmin,dmax = jmin,jmax
          label = r'$|j|$'
        elif nvar == 1:
          var = b1[:,ind,:]
          dmin,dmax = bmin,bmax
          label = r'$B_x$'                                   
        elif nvar == 2:
          var = b2[:,ind,:]
          dmin,dmax = bmin,bmax
          label = r'$B_y$'                                           
        elif nvar == 3:
          var = b3[:,ind,:]
          dmin,dmax = bmin,bmax
          label = r'$B_z$'  
        elif nvar == 4:
	  if isrho:
            var = d[:,ind,:]
            dmin,dmax = rmin,rmax
            label = r'$\rho$'
	  else:
            var = wu[:,ind,:] #(np.sqrt(w2u)/np.average(np.sqrt(w2u)))[:,yslice+step*i,:] 
            cmap = 'OrRd'
            dmin,dmax = jmin,jmax
            label = r'$|\omega|$'
        elif nvar == 5:
          var = v1[:,ind,:]
          dmin,dmax = vmin,vmax
          label = r'$v_x$'
        elif nvar == 6:
          var = v2[:,ind,:]+1.5*1.0*np.array(np.tile(x1b,(nz,1)))
          dmin,dmax = vmin,vmax
          label = r'$\delta v_y$'
        elif nvar == 7:
          var = v3[:,ind,:]
          dmin,dmax = vmin,vmax
          label = r'$v_z$'                                          
        else:
          print 'out plotting bound'

        plt.subplot(nrow,ncol,count)
        show_contour(x,y,var,dmin=dmin,dmax=dmax,xspan=xspan,yspan=yspan,cmap=cmap,xlabel='x',ylabel='z')
        plt.text(xspan[1]-0.1,yspan[1]-0.05,label,fontsize=15)                                             
        if nvar > 0:
          plt.tick_params(labelbottom='off') 
          plt.tick_params(labelleft='off')
        count += 1

    plt.tight_layout()

# plot 6-panels in x-z plane (jmag,bx,by,wmag/rho,vx,dvy)
def plot_6panel(t,x,y,j2u,w2u,b1,b2,v1,v2,d,xspan=(-1,1),yspan=(-0.5,0.5),yslice=502,nslice=1,isrho=False):
    matplotlib.rcParams['figure.figsize'] = (10,3.5)
    x1b = (x + 0.5*(x[1]-x[0]))[:-1]
    x3b = (y + 0.5*(y[1]-y[0]))[:-1]
    nx,nz = len(x1b),len(x3b)
    count = 1
    jmin,jmax = 0,5
    bmin,bmax = -2,2
    vmin,vmax = -1.5,1.5
    rmin,rmax = 0,2
    #xspan,yspan = (0.3,0.7),(-0.5,-0.1)
    ncol,nrow = 3,2

    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height

    cmap = 'seismic'
    #yslice= 502
    #step = 10

    wu = np.sqrt(w2u)
    wu /= np.average(wu)
    ju = np.sqrt(j2u)
    ju /= np.average(ju)

    ind = yslice
    for nvar in np.arange(6): #// j2,vx,bx,by
      cmap = 'seismic'
      if nvar == 0:
        var = ju[:,ind,:] #ysli(np.sqrt(j2u)/np.average(np.sqrt(j2u)))[:,yslice+step*i,:] 
        cmap = 'OrRd'
        dmin,dmax = jmin,jmax
        label = r'$|j|$'
      elif nvar == 1:
        var = b1[:,ind,:]
        dmin,dmax = bmin,bmax
        label = r'$B_x$'                                   
      elif nvar == 2:
        var = b2[:,ind,:]
        dmin,dmax = bmin,bmax
        label = r'$B_y$'                                           
      elif nvar == 3:
        if isrho:
          var = d[:,ind,:]
          dmin,dmax = rmin,rmax
          label = r'$\rho$'
        else:
          var = wu[:,ind,:] #(np.sqrt(w2u)/np.average(np.sqrt(w2u)))[:,yslice+step*i,:] 
          cmap = 'OrRd'
          dmin,dmax = jmin,jmax
          label = r'$|\omega|$'
      elif nvar == 4:
        var = v1[:,ind,:]
        dmin,dmax = vmin,vmax
        label = r'$v_x$'
      elif nvar == 5:
        var = v2[:,ind,:]+1.5*1.0*np.array(np.tile(x1b,(nz,1)))
        dmin,dmax = vmin,vmax
        label = r'$\delta v_y$'
      else:
        print 'out plotting bound'

      plt.subplot(nrow,ncol,count)
      show_contour(x,y,var,dmin=dmin,dmax=dmax,xspan=xspan,yspan=yspan,cmap=cmap,xlabel='x',ylabel='z')
      plt.text(xspan[1]-0.2,yspan[1]+0.1,label,fontsize=15)                                             
      if nvar > 0:
        plt.tick_params(labelbottom='off') 
        plt.tick_params(labelleft='off')
      count += 1

    plt.tight_layout()


def plot_exx(ncells,jmax,jmaxnorm,diss,dissnorm,size,orient,theta,disscut=4e-5,ntime=21,project=np.array([1,0,0]),cmap='jet',vsheet=False):
    """
    plot the correlation between quantities within current sheet
    """
    nframe = 3
    nplot  = 3
    enlarge = 2
    matplotlib.rcParams['figure.figsize'] = (10*enlarge, enlarge*8.0*nframe/nplot)
    fig = plt.figure()
    fraction=0.046;pad=0.04
    gs = gridspec.GridSpec(nframe, nplot)
    mask = (theta > disscut)
    # calc dissipation rate per sheet for normalization or weighted average
    ncelltot= float(np.sum(np.array(ncells)))
    diss_sheet = np.array(ncells)*np.array(diss)
    disstot= float(np.sum(diss_sheet))
    
    for i in np.arange(9):
    #######################################
    # (1) orientation
    #######################################
      if (i==0):
        vec2d = np.array(orient)[mask,0]#np.cross(np.array(orient)[mask,0],project)
        hdata_x = vec2d[:,2]
        hdata_y = vec2d[:,1]
        xlab = r'$\hat{l}_x$'
        ylab = r'$\hat{l}_y$'
        xmin,xmax = -1,1
        ymin,ymax = -1,1
	vec2d = np.array(orient)[:,0]
	vecnorm = np.average(vec2d,axis=0)
	vecnorm /= np.linalg.norm(vecnorm)
        print 'non-weighted averaged (\hat{l}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(np.array(ncells)*vec2d.T,axis=1)/ncelltot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'vol-weighted averaged (\hat{l}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(diss_sheet*vec2d.T,axis=1)/disstot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'dis-weighted averaged (\hat{l}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
      if (i==1):
        vec2d = np.array(orient)[mask,1]
        hdata_x = vec2d[:,2]
        hdata_y = vec2d[:,1]
        xlab = r'$\hat{\xi}_x$'
        ylab = r'$\hat{\xi}_y$'
        xmin,xmax = -1,1
        ymin,ymax = -1,1
	vec2d = np.array(orient)[:,1]
	vecnorm = np.average(vec2d,axis=0)
	vecnorm /= np.linalg.norm(vecnorm)
        print 'non-weighted averaged (\hat{xi}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(np.array(ncells)*vec2d.T,axis=1)/ncelltot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'vol-weighted averaged (\hat{xi}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(diss_sheet*vec2d.T,axis=1)/disstot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'dis-weighted averaged (\hat{xi}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
      if (i==2):
        vec2d = np.array(orient)[mask,2]
        hdata_x = vec2d[:,2]
        hdata_y = vec2d[:,1]
        xlab = r'$\hat{\lambda}_x$'
        ylab = r'$\hat{\lambda}_y$'
        xmin,xmax = -1,1
        ymin,ymax = -1,1
	vec2d = np.array(orient)[:,2]
	vecnorm = np.average(vec2d,axis=0)
	vecnorm /= np.linalg.norm(vecnorm)
        print 'non-weighted averaged (\hat{\lambda}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(np.array(ncells)*vec2d.T,axis=1)/ncelltot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'vol-weighted averaged (\hat{\lambda}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
	vecnorm = np.average(diss_sheet*vec2d.T,axis=1)/disstot
	vecnorm /= np.linalg.norm(vecnorm)
        print 'dis-weighted averaged (\hat{\lambda}) = (',vecnorm[2],' ,',\
                                           vecnorm[1],' ,',\
                                           vecnorm[0],')'
    #######################################
    # (2) 2d histogram of sheet dimensions
    #######################################
      if (i==3):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(size)[:,1]
        xlab = r'$l$'
        ylab = r'$\xi$'
        xmin,xmax = -2.7,1
        ymin,ymax = -2.7,0
        nbin = 100
      if (i==4):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(size)[:,2]
        xlab = r'$l$'
        ylab = r'$\lambda$'
        xmin,xmax = -2.7,1
        ymin,ymax = -2.7,0
        nbin = 100
      if (i==5):
        hdata_x = np.array(size)[:,1]
        hdata_y = np.array(size)[:,2]
        xlab = r'$\xi$'
        ylab = r'$\lambda$'
        xmin,xmax = -2.7,1
        ymin,ymax = -2.7,0
        nbin = 100
    #######################################
    # (3) dissipation below
    #######################################
      if (i==6):
        hdata_x = np.array(size)[:,0]
        hdata_y = np.array(theta)
        xlab = r'$l$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -2.7,1
        nbin = 100
      if (i==7):
        hdata_x = np.array(size)[:,1]
        hdata_y = np.array(theta)
        xlab = r'$\xi$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -2.7,1
        nbin = 100
      if (i==8):
        hdata_x = np.array(size)[:,2]
        hdata_y = np.array(theta)
        xlab = r'$\lambda$'
	if vsheet:
          ylab = r'$\theta_{\nu}/\mathcal{E}_{\nu}$'
	else:
          ylab = r'$\theta_{\eta}/\mathcal{E}_{\eta}$'
        ymin = np.log10(np.min(hdata_y))
        ymax = np.log10(np.max(hdata_y))
        xmin,xmax = -2.7,1
        nbin = 100
        
      fig.add_subplot(gs[i],aspect='equal')  
     
      if i<3:
        plt.scatter(hdata_x,hdata_y,s=0.02, marker = '.' )
        plt.ylim([ymin,ymax])
        plt.xlim([xmin,xmax])
        plt.xlabel(xlab,size=20)
        plt.ylabel(ylab,size=20)
        if i==0: # add tilt angle measured in Zhdankin paper -17.5degree
            tilt = -17.5/180.*np.pi
            rx,ry = 2.*np.sin(tilt),2.*np.cos(tilt)
            plt.plot([rx,-rx],[ry,-ry],'k:')
      elif i<6:
	#plt.scatter(hdata_x,hdata_y,s=0.01,marker='.')
	#hst are binned in the dissipation values (normalized with total dissipation)
        hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=(np.logspace(ymin,ymax,nbin),\
                                                         np.logspace(xmin,xmax,nbin)),\
                                                         weights=theta) #diss_sheet/disstot)
        xx,yy = np.meshgrid(xe,ye)
	#return hst, xe,ye,xx,yy
        loghst = np.log10(hst+1e-12)
        #plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-5,vmax=-1,cmap=cmap)
	plt.plot(xe,np.log10(np.sum(hst*xx[:-1],axis=0)+1e-12))

        #print 'min/max of hst = ',np.min(loghst),np.max(loghst)
	#print 'sum(hst) = ',np.sum(hst)
        #print 'min/max of xe = ',np.min(xe),np.max(xe)
        #print 'min/max of ye = ',np.min(ye),np.max(ye)
	
        #plt.colorbar(pad=pad,fraction=fraction*0.8)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim([10**xmin,10**xmax])
        #plt.ylim([10**ymin,10**ymax])
        plt.xlabel(xlab,size=20)
        #plt.ylabel(ylab,size=20)
	#if vsheet:
	#  plt.title(r'$\theta_{\nu}/\mathcal{E}_{\nu}$',fontsize=15)
	#else: 
	#  plt.title(r'$\theta_{\eta}/\mathcal{E}_{\eta}$',fontsize=15)
        #if i==3:
        #  index,coef=1.0,0.4
        #  plt.plot(xe,coef*xe**index,'k:')
      #plt.title(r'$\mathrm{Histogram\ of\ d_i:}$')
      else:
        hst,xe,ye = np.histogram2d(hdata_y,hdata_x,bins=[np.logspace(ymin,ymax,nbin),\
                                                         np.logspace(xmin,xmax,nbin)],\
                                                         normed=True)
        xx,yy = np.meshgrid(xe,ye)
        loghst = np.log10(hst+1e-12)
        plt.pcolormesh(yy,xx,np.transpose(loghst),vmin=-3,vmax=9)
        #print 'min/max of hst = ',np.min(loghst),np.max(loghst)
        plt.colorbar(pad=pad,fraction=fraction)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim([10**xmin,10**xmax])
        plt.ylim([10**ymin,10**ymax])
        plt.xlabel(xlab,size=20)
        plt.ylabel(ylab,size=20)
	plt.title(r'$\log$'+'(PDF)',fontsize=15)
        if i==6:
	    index,coef=2.0,1.5e-3
            plt.plot(ye,coef*ye**index,'k:')
        
    plt.tight_layout()
