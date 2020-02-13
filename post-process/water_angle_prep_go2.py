# coding: utf-8
import sys
import os.path
import numpy as np
from timeit import default_timer as timer
from datetime import datetime, timedelta

#----------------------------------------------------------------------------------------------
print("Start: ",datetime.now())
start = timer()

list_filename=["00","01","02","03","04","05"]

print("Init Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

for f in list_filename:
    filename = (f+'_AS_water_angle.txt')
    if (os.path.isfile(filename)):
        # Step C9 dist_IS_down dist_IS_up dist_AS TUDon TAcc TDon
        WA_import = np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(0,1,2,3,4,5,6,7,8,9))
        WA_import = np.hstack((np.full((WA_import.shape[0],1),int(f)),WA_import))
        try:
            _ = WA_AS_go2.shape
        except NameError:
            WA_AS_go2 = np.zeros((1,WA_import.shape[1]))
        WA_AS_go2 = np.append(WA_AS_go2, WA_import, axis=0)
        del(WA_import)
WA_AS_go2=np.delete(WA_AS_go2,(0),axis=0)

print("Load AS_water_angle.txt Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

# TrajNum O_id cOE cOH cOA cC cC9 dist_AS_down Angle OH/NAS Angle HH/NAS step
np.savez_compressed('go2_WA_AS.npz', ar1=WA_AS_go2)

print("Write go2_WA_AS.npz Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()
del(WA_AS_go2)

print("Cleanup Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()


for f in list_filename:
    filename = (f+'_IS_water_angle.txt')
    if (os.path.isfile(filename)):
        # Step C9 dist_IS_down dist_IS_up dist_AS TUDon TAcc TDon
        WA_import = np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(0,1,2,3,4,5,6,8,10,12))
        WA_import = np.hstack((np.full((WA_import.shape[0],1),int(f)),WA_import))
        try:
            _ = WA_IS_go2.shape
        except NameError:
            WA_IS_go2 = np.zeros((1,WA_import.shape[1]))
        WA_IS_go2 = np.append(WA_IS_go2, WA_import, axis=0)
        del(WA_import)
WA_IS_go2=np.delete(WA_IS_go2,(0),axis=0)

print("Load IS_water_angle.txt Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

# TrajNum O_id cOE cOH cOA cC cC9 dist_IS_down Angle OH/NIS_down Angle HH/NIS_down step
np.savez_compressed('go2_WA_IS.npz', ar1=WA_IS_go2)

print("Write go2_WA_IS.npz Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
quit()