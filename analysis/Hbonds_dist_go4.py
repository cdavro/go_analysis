# coding: utf-8
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import time
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#----------------------------------------------------------------------------------------------
print("Start: ",datetime.now())
start = time.time()

loaded = np.load('Hbonds_dist_revPBED3_go4.npz')
hbonds_revPBED3_go4 = loaded['hbonds_dist_revPBED3_go4']
print(hbonds_revPBED3_go4.shape)

print("Load hbonds_revPBED3_go4: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

plt.rcParams["figure.figsize"] = [14,10]
fig1, ax1 = plt.subplots(dpi=300)

plt.xticks(np.arange(0.5, 3.5, 0.5), fontweight='bold',fontsize=20)
plt.yticks(fontweight='bold',fontsize=20)
ax1.xaxis.set_minor_locator(MultipleLocator(0.125))
ax1.yaxis.set_minor_locator(MultipleLocator(0.125))
plt.grid()
#plt.ylim(-0.25,6.00)
plt.xlim(0.5, 3.0)
plt.setp(ax1.spines.values(), linewidth=2)

ax1.axes.tick_params(which='minor',size=4,width=2)
ax1.axes.tick_params(which='major',size=5,width=3)

ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"b"}, label=u'$H_{ALC} - O_{WAT}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==11) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"r"}, label=u'$H_{ALC} - O_{ALC}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==10) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"g"}, label=u'$H_{ALC} - O_{EPO}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==12) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"k"}, label=u'$H_{ALC} - O_{ALK}$')

ax1.legend(fontsize=20)
plt.xlabel(r'r$_{H_D-O_A}$ (Å)', fontsize=28)
plt.ylabel(r'P (r$_{H_D-O_A}$)', fontsize=28)

fig1.savefig("Hbonds_dist_ALC.png", dpi=300,bbox_inches='tight', transparent = True)
plt.clf()
plt.cla()
plt.close()

print("Hbonds_dist_ALC: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

plt.rcParams["figure.figsize"] = [14,10]
fig1, ax1 = plt.subplots(dpi=300)

plt.xticks(np.arange(0.5, 3.5, 0.5), fontweight='bold',fontsize=20)
plt.yticks(fontweight='bold',fontsize=20)
ax1.xaxis.set_minor_locator(MultipleLocator(0.125))
ax1.yaxis.set_minor_locator(MultipleLocator(0.125))
plt.grid()
#plt.ylim(-0.25,6.00)
plt.xlim(0.5, 3.0)
plt.setp(ax1.spines.values(), linewidth=2)

ax1.axes.tick_params(which='minor',size=4,width=2)
ax1.axes.tick_params(which='major',size=5,width=3)

ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"b"}, label=u'$H_{WAT} - O_{WAT}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==11) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"r"}, label=u'$H_{WAT} - O_{ALC}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==10) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"g"}, label=u'$H_{WAT} - O_{EPO}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==12) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"k"}, label=u'$H_{WAT} - O_{ALK}$')

ax1.legend(fontsize=20)
plt.xlabel(r'r$_{H_D-O_A}$ (Å)', fontsize=28)
plt.ylabel(r'P (r$_{H_D-O_A}$)', fontsize=28)

fig1.savefig("Hbonds_dist_WAT.png", dpi=300,bbox_inches='tight', transparent = True)
plt.clf()
plt.cla()
plt.close()

print("Hbonds_dist_WAT: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

plt.rcParams["figure.figsize"] = [14,10]
fig1, ax1 = plt.subplots(dpi=300)

plt.xticks(np.arange(0.5, 3.5, 0.5), fontweight='bold',fontsize=20)
plt.yticks(fontweight='bold',fontsize=20)
ax1.xaxis.set_minor_locator(MultipleLocator(0.125))
ax1.yaxis.set_minor_locator(MultipleLocator(0.125))
plt.grid()
#plt.ylim(-0.25,6.00)
plt.xlim(0.5, 3.0)
plt.setp(ax1.spines.values(), linewidth=2)

ax1.axes.tick_params(which='minor',size=4,width=2)
ax1.axes.tick_params(which='major',size=5,width=3)

#list_Ld=[-20.0,-20.0,2.75,5.75]
#List_Lu=[20.0,2.75,5.75,20.0]

#ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=-0.25) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
#                   hist_kws={'range': (0.5,3.0)},
#                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"b"}, label=u'$L0 - H_{WAT} - O_{WAT}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=2.75) & (hbonds_revPBED3_go4[:,1]>-9.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"r"}, label=u'$L1 - H_{WAT} - O_{WAT}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=5.75) & (hbonds_revPBED3_go4[:,1]>2.75) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"g"}, label=u'$L2 - H_{WAT} - O_{WAT}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==13) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>5.75) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"k"}, label=u'$L3 - H_{WAT} - O_{WAT}$')

ax1.legend(fontsize=20)
plt.xlabel(r'r$_{H_D-O_A}$ (Å)', fontsize=28)
plt.ylabel(r'P (r$_{H_D-O_A}$)', fontsize=28)

fig1.savefig("Hbonds_dist_WAT_WAT_L.png", dpi=300,bbox_inches='tight', transparent = True)
plt.clf()
plt.cla()
plt.close()

print("Hbonds_dist_WAT: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

plt.rcParams["figure.figsize"] = [14,10]
fig1, ax1 = plt.subplots(dpi=300)

plt.xticks(np.arange(0.5, 3.5, 0.5), fontweight='bold',fontsize=20)
plt.yticks(fontweight='bold',fontsize=20)
ax1.xaxis.set_minor_locator(MultipleLocator(0.125))
ax1.yaxis.set_minor_locator(MultipleLocator(0.125))
plt.grid()
#plt.ylim(-0.25,6.00)
plt.xlim(0.5, 3.0)
plt.setp(ax1.spines.values(), linewidth=2)

ax1.axes.tick_params(which='minor',size=4,width=2)
ax1.axes.tick_params(which='major',size=5,width=3)

ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]<=0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"b"}, label=u'$H_{ALC} - O_{WAT}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==11) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]<=0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"r"}, label=u'$H_{ALC} - O_{ALC}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==10) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]<=0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"g"}, label=u'$H_{ALC} - O_{EPO}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==12) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]<=0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"k"}, label=u'$H_{ALC} - O_{ALK}$')

ax1.legend(fontsize=20)
plt.xlabel(r'r$_{H_D-O_A}$ (Å)', fontsize=28)
plt.ylabel(r'P (r$_{H_D-O_A}$)', fontsize=28)

fig1.savefig("Hbonds_dist_ALC_DOWN.png", dpi=300,bbox_inches='tight', transparent = True)
plt.clf()
plt.cla()
plt.close()

print("Hbonds_dist_ALC_DOWN: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
start = time.time()

plt.rcParams["figure.figsize"] = [14,10]
fig1, ax1 = plt.subplots(dpi=300)

plt.xticks(np.arange(0.5, 3.5, 0.5), fontweight='bold',fontsize=20)
plt.yticks(fontweight='bold',fontsize=20)
ax1.xaxis.set_minor_locator(MultipleLocator(0.125))
ax1.yaxis.set_minor_locator(MultipleLocator(0.125))
plt.grid()
#plt.ylim(-0.25,10.00)
plt.xlim(0.5, 3.0)
plt.setp(ax1.spines.values(), linewidth=2)

ax1.axes.tick_params(which='minor',size=4,width=2)
ax1.axes.tick_params(which='major',size=5,width=3)

ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==13) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]>0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"b"}, label=u'$H_{ALC} - O_{WAT}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==11) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]>0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"r"}, label=u'$H_{ALC} - O_{ALC}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==10) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]>0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"g"}, label=u'$H_{ALC} - O_{EPO}$')
ax1 = sns.distplot(hbonds_revPBED3_go4[(hbonds_revPBED3_go4[:,5]==11) & (hbonds_revPBED3_go4[:,6]==12) & (hbonds_revPBED3_go4[:,1]<=9.0) & (hbonds_revPBED3_go4[:,1]>-9.0) & (hbonds_revPBED3_go4[:,3]>0.0) ][:,7],hist=False,bins=200,
                   hist_kws={'range': (0.5,3.0)},
                   kde_kws={'clip': (0.5,3.0), "linestyle":"-","linewidth":"3.0","color":"k"}, label=u'$H_{ALC} - O_{ALK}$')

ax1.legend(fontsize=20)
plt.xlabel(r'r$_{H_D-O_A}$ (Å)', fontsize=28)
plt.ylabel(r'P (r$_{H_D-O_A}$)', fontsize=28)

fig1.savefig("Hbonds_dist_ALC_UP.png", dpi=300,bbox_inches='tight', transparent = True)
plt.clf()
plt.cla()
plt.close()

print("Hbonds_dist_ALC_UP: ",datetime.now())
print("Timings: ", time.time()-start)

#----------------------------------------------------------------------------------------------
quit()