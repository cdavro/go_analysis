# coding: utf-8
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

if (len(sys.argv) == 1):
    ext='txt'
else:
    ext=(sys.argv[1])

list_filename=["00","01","02","03","04","05"]

#list_Ld=[-20.0,-20.0,2.75,5.25]
#list_Lu=[20.0,2.75,5.25,20.0]
list_Ld=[-20.0,-20.0,3.25,6.50]
list_Lu=[20.0,3.25,6.50,20.0]

# Bins de merde (1° width)
bin_0_180_01=[1,0.99985,0.99939,0.99863,0.99756,0.99619,0.99452,0.99255,0.99027,0.98769,0.98481,0.98163,0.97815,0.97437,0.9703,0.96593,0.96126,0.9563,0.95106,0.94552,0.93969,0.93358,0.92718,0.9205,0.91355,0.90631,0.89879,0.89101,0.88295,0.87462,0.86603,0.85717,0.84805,0.83867,0.82904,0.81915,0.80902,0.79864,0.78801,0.77715,0.76604,0.75471,0.74314,0.73135,0.71934,0.70711,0.69466,0.682,0.66913,0.65606,0.64279,0.62932,0.61566,0.60182,0.58779,0.57358,0.55919,0.54464,0.52992,0.51504,0.5,0.48481,0.46947,0.45399,0.43837,0.42262,0.40674,0.39073,0.37461,0.35837,0.34202,0.32557,0.30902,0.29237,0.27564,0.25882,0.24192,0.22495,0.20791,0.19081,0.17365,0.15643,0.13917,0.12187,0.10453,0.08716,0.06976,0.05234,0.0349,0.01745,0,-0.01745,-0.0349,-0.05234,-0.06976,-0.08716,-0.10453,-0.12187,-0.13917,-0.15643,-0.17365,-0.19081,-0.20791,-0.22495,-0.24192,-0.25882,-0.27564,-0.29237,-0.30902,-0.32557,-0.34202,-0.35837,-0.37461,-0.39073,-0.40674,-0.42262,-0.43837,-0.45399,-0.46947,-0.48481,-0.5,-0.51504,-0.52992,-0.54464,-0.55919,-0.57358,-0.58779,-0.60182,-0.61566,-0.62932,-0.64279,-0.65606,-0.66913,-0.682,-0.69466,-0.70711,-0.71934,-0.73135,-0.74314,-0.75471,-0.76604,-0.77715,-0.78801,-0.79864,-0.80902,-0.81915,-0.82904,-0.83867,-0.84805,-0.85717,-0.86603,-0.87462,-0.88295,-0.89101,-0.89879,-0.90631,-0.91355,-0.9205,-0.92718,-0.93358,-0.93969,-0.94552,-0.95106,-0.9563,-0.96126,-0.96593,-0.9703,-0.97437,-0.97815,-0.98163,-0.98481,-0.98769,-0.99027,-0.99255,-0.99452,-0.99619,-0.99756,-0.99863,-0.99939,-0.99985,-1]
bin_0_180_01=np.flip(np.asarray(bin_0_180_01))

WA_revPBED3=np.zeros((1,4))
for f in list_filename:
    filename = (f+'_water-angle.'+ext)
    if (os.path.isfile(filename)):
        WA_revPBED3 = np.append(WA_revPBED3, np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(1,2,4,6)), axis=0)
WA_revPBED3=np.delete(WA_revPBED3,(0),axis=0)

plt.rcParams["figure.figsize"] = [14,7]
fig1, ax1 = plt.subplots(dpi=300)
sns.distplot(np.cos(WA_revPBED3[(WA_revPBED3[:,0]==1) & (WA_revPBED3[:,1]>list_Ld[1]) &(WA_revPBED3[:,1]<=list_Lu[1])][:,2])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1) , "linestyle":"-","linewidth":"3.0","color":"gold"}, label=r'$L1$')
sns.distplot(np.cos(WA_revPBED3[(WA_revPBED3[:,0]==1) & (WA_revPBED3[:,1]>list_Ld[2]) &(WA_revPBED3[:,1]<=list_Lu[2])][:,2])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1) , "linestyle":"--","linewidth":"3.0","color":"red"}, label=r'$L2$')
sns.distplot(np.cos(WA_revPBED3[(WA_revPBED3[:,0]==1) & (WA_revPBED3[:,1]>list_Ld[3]) &(WA_revPBED3[:,1]<=list_Lu[3])][:,2])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1), "linestyle":":","linewidth":"3.0","color":"black"}, label=r'$L3$')
plt.ylabel(r'P(cos(θ$_{DW}$))', fontsize=28)
plt.xlabel(r'cos(θ$_{DW}$)', fontsize=28)
ax1.legend(fontsize=20,framealpha=0.3,fancybox=True)
plt.yticks(fontweight='bold',fontsize=20)
plt.xlim(-1, 1)
#plt.ylim(0, 1.2)
plt.xticks(np.arange(-1, 1.25, 0.25), fontweight='bold',fontsize=20)
#plt.yticks(np.arange(0,1.4,0.2), fontweight='bold',fontsize=12)
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
#ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
plt.setp(ax1.spines.values(), linewidth=4)
ax1.axes.tick_params(which='minor',size=6,width=2)
ax1.axes.tick_params(which='major',size=8,width=4)
ax1.grid(axis='x',which='minor', alpha=0.3)
ax1.grid(axis='both',which='major', alpha=0.6)
ax1.xaxis.labelpad = 15
ax1.yaxis.labelpad = 15
plt.tight_layout()
plt.savefig('figure_cos_DW_revPBED3', bbox_inches='tight', transparent = True, dpi = 300)

quit()