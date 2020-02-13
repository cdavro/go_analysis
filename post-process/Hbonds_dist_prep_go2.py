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
    filename = (f+'_OH_Hbonds.txt')
    if (os.path.isfile(filename)):
        # Step C9 dist_IS_down dist_IS_up dist_AS TUDon TAcc TDon
        H1_import = np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(3,11,15,16,17,14,12,13))
        H1_import = np.hstack((np.full((H1_import.shape[0],1),int(f)),H1_import))
        try:
            _ = H1.shape
        except NameError:
            H1 = np.zeros((1,H1_import.shape[1]))
        H1 = np.append(H1, H1_import, axis=0)
        del(H1_import)
H1=np.delete(H1,(0),axis=0)

print("Load OH_Hbonds.txt Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

for f in list_filename:
    filename = (f+'_OH_Hbonds_dist.txt')
    if (os.path.isfile(filename)):
        # OH_id O_type O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist
        H2_import = np.genfromtxt(filename,dtype='float64',skip_header=1)
        try:
            _ = H2.shape
        except NameError:
            H2 = np.zeros((1,H2_import.shape[1]))
        H2 = np.append(H2, H2_import, axis=0)
        del(H2_import)
H2=np.delete(H2,(0),axis=0)

print("Load OH_Hbonds_dist.txt Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

HStacked=np.hstack((H1,H2))
del(H1,H2)

print("HStack Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

# TrajNum Step C9 dist_IS_down dist_IS_up dist_AS TUDon TAcc TDon OH_id O_type O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist O_id O_Type Dist
np.savez_compressed('go2_HB_TD_full.npz', ar1=HStacked)

print("Write go2_HB_TD_full.npz Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

HStacked = HStacked[(HStacked[:,1] % 10) == 0]

print("Modulo Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

H_in=np.concatenate((HStacked[:,:14],
np.hstack((HStacked[HStacked[:,14] !=0 ][:,:11],HStacked[HStacked[:,14] !=0 ][:,14:17])),
np.hstack((HStacked[HStacked[:,17] !=0 ][:,:11],HStacked[HStacked[:,17] !=0 ][:,17:20])),
np.hstack((HStacked[HStacked[:,20] !=0 ][:,:11],HStacked[HStacked[:,20] !=0 ][:,20:23]))))

del(HStacked)

print("Inline Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

# TrajNum Step C9 dist_IS_down dist_IS_up dist_AS TUDon TAcc TDon OH_id O_type O_id O_Type Dist
np.savez_compressed('go2_HB_TD_inline_mod10.npz', ar1=H_in)

print("Write go2_HB_TD_inline_mod10.npz Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
print("End: ",datetime.now())
quit()