# coding: utf-8
import sys
import os.path
import numpy as np
from timeit import default_timer as timer
from datetime import datetime, timedelta

#----------------------------------------------------------------------------------------------
print("Start: ",datetime.now())
start = timer()

hbonds_revPBED3_go4 = np.load('go4_HB_TD_inline_mod10.npz')['ar1']

print("Load go4_HB_TD_inline_mod10.npz Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()

BIGONE=np.zeros((1,16))
donerrors=0
accerrors=0
herrors=0
dherrors=0
tj=np.unique(np.sort(hbonds_revPBED3_go4[:,0],axis=0)).astype(int)
j0=0
for j in tj:
    # 0TrajNum 1Step 2C9 3dist_IS_down 4dist_IS_up 5dist_AS 6TUDon 7TAcc 8TDon 9OH_id 10H_id 11O_type 12O_id 13O_Type Dist
    H1=hbonds_revPBED3_go4[hbonds_revPBED3_go4[:,0] == j]
    ts=np.unique(np.sort(H1[:,1],axis=0)).astype(int)
    s0=0
    for s in ts:
        H2=H1[H1[:,1] == s]
        list = np.sort(H2[:,9],axis=0).astype(int)
        listu = np.unique(list).astype(int)
        count = list.shape[0]
        GOP=np.hstack((np.expand_dims(listu,1),np.zeros((listu.shape[0],15))))
        del(listu,list)
        # 0-O_ID, 1-H_ID1, 2-H_ID2, 3-O_TYPE, 4-TUDon, 5-TDon, 6-DOG1, 7-DOG2, 8-DWAT1, 9-DWAT2, 10-TAcc, 11-AOH, 12-AWAT, 13-DIST_IS, 14-TRAJ, 15-TIMESTEP
        for n in np.arange(count):
            locOS=np.where(GOP[:,0] == H2[n,9]) #ID_SOURCE
            locOT=np.where(GOP[:,0] == H2[n,12]) #ID_TARGET

            GOP[locOS,3] = H2[n,11] # O_Type
            GOP[locOS,4] = H2[n,6] # TUDON
            GOP[locOS,5] = H2[n,8] # TDON
            GOP[locOS,10] = H2[n,7] # TACC
            GOP[locOS,13] = H2[n,3]
            GOP[locOS,14] = tj[j0]
            GOP[locOS,15] = ts[s0]

            if (H2[n,10] == GOP[locOS,2]):
                GOP[locOS,2] = H2[n,10]
            elif (H2[n,10] == GOP[locOS,1]):
                GOP[locOS,1] = H2[n,10]
            elif (GOP[locOS,1] != 0) & (H2[n,10] != GOP[locOS,1]):
                GOP[locOS,2] = H2[n,10]
            elif (GOP[locOS,1] == 0):
                GOP[locOS,1] = H2[n,10]
            else:
                herrors=herrors+1

            if (H2[n,13] == 0):
                #print('EXIT')
                continue
            if GOP[locOS,1] == H2[n,10]:
                if (H2[n,13] != 13):
                    GOP[locOS,6] = GOP[locOS,6] + 1
                else:
                    GOP[locOS,8] = GOP[locOS,8] + 1
            elif GOP[locOS,2] == H2[n,10]:
                if (H2[n,13] != 13):
                    GOP[locOS,7] = GOP[locOS,7] + 1
                else:
                    GOP[locOS,9] = GOP[locOS,9] + 1
            else:
                dherrors=dherrors+1

            if (H2[n,11] != 13):
                GOP[locOT,11] = GOP[locOT,11] + 1
            else:
                GOP[locOT,12] = GOP[locOT,12] + 1

        for n in np.arange(GOP.shape[0]):
            if GOP[n,5] != (GOP[n,6]+GOP[n,7]+GOP[n,8]+GOP[n,9]):
                donerrors=donerrors+1
            if GOP[n,10] != (GOP[n,11]+GOP[n,12]):
                accerrors=accerrors+1

        BIGONE = np.append(BIGONE, GOP, axis=0)
        del(GOP,H2)
        s0=s0+1
    j0=j0+1
BIGONE=np.delete(BIGONE,(0),axis=0)

print("donerrors",donerrors)
print("accerrors",accerrors)
print("herrors",herrors)
print("dherrors",dherrors)

print("Array transformation Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
start = timer()
np.savez_compressed('go4_HB_TDATSIS_mod10.npz', ar1=BIGONE)

print("Write go4_HB_TDATSIS_mod10.npz Timings: ", timedelta(seconds=timer()-start))

#----------------------------------------------------------------------------------------------
print("End: ",datetime.now())
quit()