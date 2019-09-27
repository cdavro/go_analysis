import numpy as np

hbonds_revPBED3={}

list_filename=["01","02","03","04","05"]

list_type=["all","0","1N","1S","2S","3S","2D","3D","4D","5D"]
list_donnor=[-1,0,0,1,1,1,2,2,2,2]
list_acceptor=[-1,0,1,0,1,2,0,1,2,3]


list_L=["all","L1","L2","L3"]
list_Ld=[-20.0,-20.0,2.75,5.25]
list_Lu=[20.0,2.75,5.25,20.0]

for i in list_filename:
    with open (i+'_O_hbonds.dat') as f:
        hbonds_revPBED3[i] = np.genfromtxt(f,dtype='float64',skip_header=1)

oh={}
for n,o,p in zip(list_type,list_donnor,list_acceptor):
    oh[n]={}
    for k,l,m in zip(list_L,list_Ld,list_Lu):

        oh[n][k]=0
        for i in list_filename:
            if n == 'all':
                oh[n][k] = oh[n][k] + np.count_nonzero(hbonds_revPBED3[i][(hbonds_revPBED3[i][:,10]==13)
                    & (hbonds_revPBED3[i][:,9]==1)
                    & (hbonds_revPBED3[i][:,11]>l)
                    & (hbonds_revPBED3[i][:,11]<=m)
                    ][:,1])
            else:
                oh[n][k] = oh[n][k] + np.count_nonzero(hbonds_revPBED3[i][(hbonds_revPBED3[i][:,10]==13)
                    & (hbonds_revPBED3[i][:,9]==1)
                    & (hbonds_revPBED3[i][:,11]>l)
                    & (hbonds_revPBED3[i][:,11]<=m)
                    & (hbonds_revPBED3[i][:,3]==o)
                    & (hbonds_revPBED3[i][:,4]==p)
                    ][:,1])
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
        for i in list_filename:
            if k == 'all':
                oh[n][k] = oh[n][k] + np.count_nonzero(hbonds_revPBED3[i][(hbonds_revPBED3[i][:,10]==13)
                    & (hbonds_revPBED3[i][:,9]==1)
                    & (hbonds_revPBED3[i][:,11]>o)
                    & (hbonds_revPBED3[i][:,11]<=p)
                    ][:,1])
            else:
                oh[n][k] = oh[n][k] + np.count_nonzero(hbonds_revPBED3[i][(hbonds_revPBED3[i][:,10]==13)
                    & (hbonds_revPBED3[i][:,9]==1)
                    & (hbonds_revPBED3[i][:,11]>o)
                    & (hbonds_revPBED3[i][:,11]<=p)
                    & (hbonds_revPBED3[i][:,3]==l)
                    & (hbonds_revPBED3[i][:,4]==m)
                    ][:,1])
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