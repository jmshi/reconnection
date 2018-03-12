#!/bin/bash
#nohup python dump_hist2d.py 'x4y4z1r64pm0.5re3000' &
#nohup python dump_hist2d.py 'x4y4z1r64pm0.5re6000' &
#nohup python dump_hist2d.py 'x4y4z1r64pm1re3000'   &
#nohup python dump_hist2d.py 'x4y4z1r64pm2re3000'   &

#python pspec.py x2y4z1r256ideal 25 51 1 &>log.r256.pspec &
python dump_hist2d.py x2y4z1r256ideal 0 51 >& log.r512.hist &
python cdf_jsheet.py  x2y4z1r256ideal 30 51 1 >& log.r512.jcdf &
python cdf_vsheet2.py  x2y4z1r256ideal 30 51 1 >& log.r512.vcdf &
python main.py x2y4z1r256ideal 30 51 1 >& log.r512.jlist &
python main_vsheet.py x2y4z1r256ideal 30 51 1 >& log.r512.vlist &
python get_jsheet_prop.py x2y4z1r256ideal  30 51 1 1.25e-4,8 >& log.r512.jprop &
python get_vsheet_prop.py x2y4z1r256ideal 30 51 1 1.25e-4,8 >& log.r512.vprop &
