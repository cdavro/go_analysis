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

for f in list_filename:
    filename = (f+'_O_Hbonds.txt')
    if (os.path.isfile(filename)):
        # Step O_type Acc UDon C9 dist_IS_down
        H1_import = np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(1,2,3,5,10,11))
        try:
            _ = H1.shape
        except NameError:
            H1 = np.zeros((1,H1_import.shape[1]))
        H1 = np.append(H1, H1_import, axis=0)
H1=np.delete(H1,(0),axis=0)

print("Load O_Hbonds.txt Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

# Step O_type Acc UDon C9 dist_IS_down
np.savez_compressed('go4_HB_TADCIS.npz', ar1=H1)

print("Write go4_HB_TADCIS.npz Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
quit()