# Magnetic Reconnection 

  * for 2d reconnection:
      python reconnect_2d.py s1e3z.ideal.np 0 100 3 0.05 &>log.s3e5z &
      python plot_2dcs.py s1e3z.ideal.np  0 100 >& log.s1e5z.ideal &
      use get_phimax.py to get the phi_max
      use notebook two_dimensional_reconnection.ipynb as examples

  * for 3d reconnection: 
      main*.py  identify all current sheets for given frames
      null_finder.py  identify all magnetic null pts for given frames
      get_jsheet*.py  extract useful properties of current sheets 
      use mri_current_sheet.ipynb as an user guide for analyzing data
     
  * for MRI 2d test project: 
      use MRI2d_Test.ipynb and mri2d_test.py to analyze/visualize 
      output from a designed MRI 2D test problem

      




