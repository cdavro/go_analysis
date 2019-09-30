# coding: utf-8
import sys
import os.path
import numpy as np

if ((sys.argv[1]) == ''):
    ext='txt'
else:
    ext=(sys.argv[1])

list_filename=["00","01","02","03","04","05"]

list_type=["all","0","1N","1S","2S","3S","2D","3D","4D","5D"]
list_donnor=[-1,0,0,1,1,1,2,2,2,2]
list_acceptor=[-1,0,1,0,1,2,0,1,2,3]

list_L=["all","L1","L2","L3"]
#list_Ld=[-20.0,-20.0,2.75,5.25]
#list_Lu=[20.0,2.75,5.25,20.0]
list_Ld=[-20.0,-20.0,3.25,6.50]
list_Lu=[20.0,3.25,6.50,20.0]

hbonds=np.zeros((1,6))
for f in list_filename:
    filename = (f+'_O_hbonds.'+ext)
    if (os.path.isfile(filename)):
        hbonds = np.append(hbonds, np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(1,3,4,9,10,11)), axis=0)
hbonds=np.delete(hbonds,(0),axis=0)

oh={}
for n,o,p in zip(list_type,list_donnor,list_acceptor):
    oh[n]={}
    for k,l,m in zip(list_L,list_Ld,list_Lu):
        oh[n][k]=0
        if n == 'all':
            oh[n][k] = oh[n][k] + np.count_nonzero(hbonds[(hbonds[:,4]==13)
                & (hbonds[:,3]==1)
                & (hbonds[:,5]>l)
                & (hbonds[:,5]<=m)
                ][:,0])
        else:
            oh[n][k] = oh[n][k] + np.count_nonzero(hbonds[(hbonds[:,4]==13)
                & (hbonds[:,3]==1)
                & (hbonds[:,5]>l)
                & (hbonds[:,5]<=m)
                & (hbonds[:,1]==o)
                & (hbonds[:,2]==p)
                ][:,0])
print(oh)
for n,o,p in zip(list_type,list_donnor,list_acceptor):
    i = 0
    for k,l,m in zip(list_L,list_Ld,list_Lu):
        print(n,list(oh[n].keys())[i])
        if oh[n]['all'] != 0:
            print(oh[n][k],oh[n][k]/oh[n]['all']*100)
        else:
            print(oh[n][k],0)
        i = i+1
        print('----')

oh={}
for n,o,p in zip(list_L,list_Ld,list_Lu):
    oh[n]={}
    for k,l,m in zip(list_type,list_donnor,list_acceptor):
        oh[n][k]=0
        if k == 'all':
            oh[n][k] = oh[n][k] + np.count_nonzero(hbonds[(hbonds[:,4]==13)
                & (hbonds[:,3]==1)
                & (hbonds[:,5]>o)
                & (hbonds[:,5]<=p)
                ][:,0])
        else:
            oh[n][k] = oh[n][k] + np.count_nonzero(hbonds[(hbonds[:,4]==13)
                & (hbonds[:,3]==1)
                & (hbonds[:,5]>o)
                & (hbonds[:,5]<=p)
                & (hbonds[:,1]==l)
                & (hbonds[:,2]==m)
                ][:,0])
print(oh)

for n,o,p in zip(list_L,list_Ld,list_Lu):
    i = 0
    for k,l,m in zip(list_type,list_donnor,list_acceptor):
        print(n,list(oh[n].keys())[i])
        if oh[n]['all'] != 0:
            print(oh[n][k],oh[n][k]/oh[n]['all']*100)
        else:
            print(oh[n][k],0)
        i = i+1
        print('----')

quit()