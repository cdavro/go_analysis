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

hbonds_2=np.zeros((1,6))
for f in list_filename:
    filename = (f+'_OH_Hbonds.txt')
    if (os.path.isfile(filename)):
        # Step C9 dist_IS_down dist_IS_up dist_AS
        in_v = np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(3,11,15,16,17))
        in_c = np.full((in_v.shape[0],1),int(f))
        in_z = np.hstack((in_c,in_v))
        del(in_c,in_v)
        hbonds_2 = np.append(hbonds_2, in_z, axis=0)
        del(in_z)
hbonds_2=np.delete(hbonds_2,(0),axis=0)

print("hbonds_2: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

hbonds_1=np.zeros((1,14))
for f in list_filename:
    filename = (f+'_OH_Hbonds_dist.txt')
    if (os.path.isfile(filename)):
        # OH_id O_type O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist
        hbonds_1 = np.append(hbonds_1, np.genfromtxt(filename,dtype='float64',skip_header=1),axis=0)
hbonds_1=np.delete(hbonds_1,(0),axis=0)

print("hbonds_1: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

hbonds=np.hstack((hbonds_2,hbonds_1))
del(hbonds_2,hbonds_1)

print("hbonds: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

# TrajNum Step C9 dist_IS_down dist_IS_up dist_AS OH_id O_type O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist
np.savez_compressed('go4_HB_TD_full.npz', ar1=hbonds)

print("hbonds: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

hbonds = hbonds[::10,:]

print("hbonds_modulo: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

hbonds_0 = hbonds[:,:8]
hbonds_1 = hbonds[:,8:11]
hbonds_2 = hbonds[:,11:14]
hbonds_3 = hbonds[:,14:17]
hbonds_4 = hbonds[:,17:20]
hbonds_5=np.hstack((hbonds_0,hbonds_1))
hbonds_6=np.hstack((hbonds_0,hbonds_2))
hbonds_7=np.hstack((hbonds_0,hbonds_3))
hbonds_8=np.hstack((hbonds_0,hbonds_4))
del(hbonds_0,hbonds_1,hbonds_2,hbonds_3,hbonds_4)
hbonds_v=np.concatenate((hbonds_5,hbonds_6,hbonds_7,hbonds_8))
del(hbonds_5,hbonds_6,hbonds_7,hbonds_8)
hbonds_n = hbonds_v[hbonds_v[:,9] !=0 ]
del(hbonds_v)

print(hbonds_n.shape)
print("hbonds_n: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

# TrajNum Step C9 dist_IS_down dist_IS_up dist_AS OH_id O_type O_id O_Type Dist
# TrajNum Step C9 dist_IS_down dist_IS_up dist_AS OH_id O_type O_id O_Type Dist
# TrajNum Step C9 dist_IS_down dist_IS_up dist_AS OH_id O_type O_id O_Type Dist
# TrajNum Step C9 dist_IS_down dist_IS_up dist_AS OH_id O_type O_id O_Type Dist
np.savez_compressed('go4_HB_TD_inline_mod10.npz', ar1=hbonds_n)

print("write: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
quit()