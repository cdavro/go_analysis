# coding: utf-8
import sys
import os.path
import numpy as np
from datetime import datetime
import time

#----------------------------------------------------------------------------------------------
print("Start: ",datetime.now())
start = time.time()

list_filename=["00","01","02","03","04","05"]

print("Init: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

WA_AS_revPBED3_go4=np.zeros((1,11))
for f in list_filename:
    filename = (f+'_AS_water_angle.txt')
    if (os.path.isfile(filename)):
        # O_id cOE cOH cOA cC cC9 dist_AS_down Angle OH/NAS Angle HH/NAS step
        in_v = np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(0,1,2,3,4,5,6,7,8,9))
        in_c = np.full((in_v.shape[0],1),int(f))
        in_z = np.hstack((in_c,in_v))
        del(in_c,in_v)
        WA_AS_revPBED3_go4 = np.append(WA_AS_revPBED3_go4, in_z, axis=0)
WA_AS_revPBED3_go4=np.delete(WA_AS_revPBED3_go4,(0),axis=0)

print("Read AS: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

# TrajNum O_id cOE cOH cOA cC cC9 dist_AS_down Angle OH/NAS Angle HH/NAS step
np.savez_compressed('go4_WA_AS.npz', ar1=WA_AS_revPBED3_go4)

print("Write_npz  AS: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
del(WA_AS_revPBED3_go4)

#----------------------------------------------------------------------------------------------
print("Start: ",datetime.now())
start = time.time()

list_filename=["00","01","02","03","04","05"]

print("Init: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

WA_IS_revPBED3_go4=np.zeros((1,11))
for f in list_filename:
    filename = (f+'_IS_water_angle.txt')
    if (os.path.isfile(filename)):
        # O_id cOE cOH cOA cC cC9 dist_IS_down Angle OH/NIS_down Angle HH/NIS_down step
        in_v = np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(0,1,2,3,4,5,6,8,10,12))
        in_c = np.full((in_v.shape[0],1),int(f))
        in_z = np.hstack((in_c,in_v))
        del(in_c,in_v)
        WA_IS_revPBED3_go4 = np.append(WA_IS_revPBED3_go4, in_z, axis=0)
WA_IS_revPBED3_go4=np.delete(WA_IS_revPBED3_go4,(0),axis=0)

print("Read IS: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

# TrajNum O_id cOE cOH cOA cC cC9 dist_IS_down Angle OH/NIS_down Angle HH/NIS_down step
np.savez_compressed('go4_WA_IS.npz', ar1=WA_IS_revPBED3_go4)

print("Write_npz  IS: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
quit()