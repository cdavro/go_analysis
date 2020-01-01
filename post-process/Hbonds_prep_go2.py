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

hbonds=np.zeros((1,6))
for f in list_filename:
    filename = (f+'_O_Hbonds.txt')
    if (os.path.isfile(filename)):
        hbonds = np.append(hbonds, np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(1,2,3,5,10,11)), axis=0)

hbonds=np.delete(hbonds,(0),axis=0)

print("Read: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

np.savez_compressed('Hbonds_revPBED3_go2.npz', ar1=hbonds)

print("write: ",datetime.now())
print("Timings: ", time.time()-start)


#----------------------------------------------------------------------------------------------
quit()