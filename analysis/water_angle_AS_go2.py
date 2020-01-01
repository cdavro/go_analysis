# coding: utf-8
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter, AutoMinorLocator

if (len(sys.argv) == 1):
    ext='txt'
else:
    ext=(sys.argv[1])

list_filename=["00","01","02","03","04","05"]

bin_0_180_01=[]
for i in np.arange(0,181,1):
    bin_0_180_01.append(np.cos(np.radians(i)))
bin_0_180_01=np.flip(np.asarray(bin_0_180_01))

bin_90_180_01=[]
for i in np.arange(90,181,1):
    bin_90_180_01.append(np.cos(np.radians(i)))
bin_0_180_01=np.flip(np.asarray(bin_90_180_01))

bin_0_180=[]
for i in np.arange(0,185,5):
    bin_0_180.append(np.cos(np.radians(i)))
bin_0_180=np.flip(np.asarray(bin_0_180))

bin_90_180=[]
for i in np.arange(90,185,5):
    bin_90_180.append(np.cos(np.radians(i)))
bin_90_180=np.flip(np.asarray(bin_90_180))

WA_revPBED3=np.zeros((1,7))
for f in list_filename:
    filename = (f+'_AS_water_angle.'+ext)
    if (os.path.isfile(filename)):
        WA_revPBED3 = np.append(WA_revPBED3, np.genfromtxt(filename,dtype='float64',skip_header=1,usecols=(1,2,3,4,5,7,8)), axis=0)
WA_revPBED3=np.delete(WA_revPBED3,(0),axis=0)
print(WA_revPBED3)

# ---------------------------------------------------------------------------------------------------------------------------------------
plt.rcParams["figure.figsize"] = [14,7]
fig1, ax1 = plt.subplots(dpi=300)
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,5])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1) , "linestyle":"-","linewidth":"3.0","color":"gold"}, label=r'$gol$')
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,5])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1) , "linestyle":"--","linewidth":"3.0","color":"red"}, label=r'$gl$')
sns.distplot(np.cos(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,5])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1), "linestyle":":","linewidth":"3.0","color":"black"}, label=r'$bl$')

plt.setp(ax1.spines.values(), linewidth=4)
ax1.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.axes.tick_params(which='minor',size=6,width=2)
ax1.axes.tick_params(which='major',size=8,width=4)
ax1.grid(axis='both',which='minor', alpha=0.2)
ax1.grid(axis='both',which='major', alpha=0.6)
ax1.xaxis.labelpad = 15
ax1.yaxis.labelpad = 15

ax1.set_xlim(-1,1)
#ax1.set_ylim(-1,1)

ax1.axes.set_xticks(np.arange(-1, 1.25, 0.25))
ax1.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 1.25, 0.25)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_xlabel(r'cos(θ$_{DW}$)',{'fontsize':24,'fontweight':'normal'})

ax1.axes.set_yticks(ax1.axes.get_yticks())
ax1.axes.set_yticklabels(['{0:.2f}'.format(x) for x in ax1.axes.get_yticks()],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_ylabel(r'P(cos(θ$_{DW}$))',{'fontsize':25,'fontweight':'normal'})

plt.savefig('figure_AS_cos_DW_T_revPBED3', bbox_inches='tight', transparent = True, dpi = 300)
plt.clf()
plt.cla()
plt.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
plt.rcParams["figure.figsize"] = [14,7]
fig1, ax1 = plt.subplots(dpi=300)
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,6])
    ,hist=False,bins=bin_90_180_01,kde_kws={'clip': (-1,0) , "linestyle":"-","linewidth":"3.0","color":"gold"}, label=r'$gol$')
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,6])
    ,hist=False,bins=bin_90_180_01,kde_kws={'clip': (-1,0) , "linestyle":"--","linewidth":"3.0","color":"red"}, label=r'$gl$')
sns.distplot(np.cos(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,6])
    ,hist=False,bins=bin_90_180_01,kde_kws={'clip': (-1,0), "linestyle":":","linewidth":"3.0","color":"black"}, label=r'$bl$')

plt.setp(ax1.spines.values(), linewidth=4)
ax1.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.axes.tick_params(which='minor',size=6,width=2)
ax1.axes.tick_params(which='major',size=8,width=4)
ax1.grid(axis='both',which='minor', alpha=0.2)
ax1.grid(axis='both',which='major', alpha=0.6)
ax1.xaxis.labelpad = 15
ax1.yaxis.labelpad = 15

ax1.set_xlim(-1,1)
#ax1.set_ylim(-1,1)

ax1.axes.set_xticks(np.arange(-1, 0.25, 0.25))
ax1.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 0.25, 0.25)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_xlabel(r'cos(θ$_{HH}$)',{'fontsize':24,'fontweight':'normal'})

ax1.axes.set_yticks(ax1.axes.get_yticks())
ax1.axes.set_yticklabels(['{0:.2f}'.format(x) for x in ax1.axes.get_yticks()],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_ylabel(r'P(cos(θ$_{HH}$))',{'fontsize':25,'fontweight':'normal'})

plt.savefig('figure_AS_cos_HH_T_revPBED3', bbox_inches='tight', transparent = True, dpi = 300)
plt.clf()
plt.cla()
plt.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
plt.rcParams["figure.figsize"] = [14,14]
fig1, ax = plt.subplots(dpi=300)

ax1 = plt.subplot(211)
sns.distplot(np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,5])
    ,hist=False,bins=180,kde_kws={'clip': (0,180) , "linestyle":"-","linewidth":"3.0","color":"gold"}, label=r'$gol$')
sns.distplot(np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,5])
    ,hist=False,bins=180,kde_kws={'clip': (0,180) , "linestyle":"--","linewidth":"3.0","color":"red"}, label=r'$gl$')
sns.distplot(np.degrees(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,5])
    ,hist=False,bins=180,kde_kws={'clip': (0,180), "linestyle":":","linewidth":"3.0","color":"black"}, label=r'$bl$')

plt.setp(ax1.spines.values(), linewidth=4)
ax1.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.axes.tick_params(which='minor',size=6,width=2)
ax1.axes.tick_params(which='major',size=8,width=4)
ax1.grid(axis='both',which='minor', alpha=0.2)
ax1.grid(axis='both',which='major', alpha=0.6)
ax1.xaxis.labelpad = 15
ax1.yaxis.labelpad = 15

ax1.set_xlim(0,180)
#ax1.set_ylim(-1,1)

ax1.axes.set_xticks(np.arange(0, 200, 20))
ax1.axes.set_xticklabels(['{0:.0f}'.format(x) for x in np.arange(0, 200, 20)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax1.axes.set_yticks(ax1.axes.get_yticks())
ax1.axes.set_yticklabels(['{0:.3f}'.format(x) for x in ax1.axes.get_yticks()],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_ylabel(r'P(θ$_{DW}$)',{'fontsize':25,'fontweight':'normal'})


ax2 = plt.subplot(212)
sns.distplot(np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,6])
    ,hist=False,bins=90,kde_kws={'clip': (90,180) , "linestyle":"-","linewidth":"3.0","color":"gold"}, label=r'$gol$')
sns.distplot(np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,6])
    ,hist=False,bins=90,kde_kws={'clip': (90,180) , "linestyle":"--","linewidth":"3.0","color":"red"}, label=r'$gl$')
sns.distplot(np.degrees(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,6])
    ,hist=False,bins=90,kde_kws={'clip': (90,180), "linestyle":":","linewidth":"3.0","color":"black"}, label=r'$bl$')

plt.setp(ax2.spines.values(), linewidth=4)
ax2.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.axes.tick_params(which='minor',size=6,width=2)
ax2.axes.tick_params(which='major',size=8,width=4)
ax2.grid(axis='both',which='minor', alpha=0.2)
ax2.grid(axis='both',which='major', alpha=0.6)
ax2.xaxis.labelpad = 15
ax2.yaxis.labelpad = 15

ax2.set_xlim(0,180)
#ax1.set_ylim(-1,1)

ax2.axes.set_xticks(np.arange(0, 200, 20))
ax2.axes.set_xticklabels(['{0:.0f}'.format(x) for x in np.arange(0, 200, 20)],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_xlabel(r'θ$_{HH}$',{'fontsize':24,'fontweight':'normal'})

ax2.axes.set_yticks(ax2.axes.get_yticks())
ax2.axes.set_yticklabels(['{0:.3f}'.format(x) for x in ax2.axes.get_yticks()],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})

plt.savefig('figure_AS_DW_HH_T_revPBED3', bbox_inches='tight', transparent = True, dpi = 300)
plt.clf()
plt.cla()
plt.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
plt.rcParams["figure.figsize"] = [14,14]
fig1, ax = plt.subplots(dpi=300)

ax1 = plt.subplot(211)
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,5])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1) , "linestyle":"-","linewidth":"3.0","color":"gold"}, label=r'$gol$')
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,5])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1) , "linestyle":"--","linewidth":"3.0","color":"red"}, label=r'$gl$')
sns.distplot(np.cos(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,5])
    ,hist=False,bins=bin_0_180_01,kde_kws={'clip': (-1,1), "linestyle":":","linewidth":"3.0","color":"black"}, label=r'$bl$')

plt.setp(ax1.spines.values(), linewidth=4)
ax1.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.axes.tick_params(which='minor',size=6,width=2)
ax1.axes.tick_params(which='major',size=8,width=4)
ax1.grid(axis='both',which='minor', alpha=0.2)
ax1.grid(axis='both',which='major', alpha=0.6)
ax1.xaxis.labelpad = 15
ax1.yaxis.labelpad = 15

ax1.set_xlim(-1,1)
#ax1.set_ylim(-1,1)

ax1.axes.set_xticks(np.arange(-1, 1.25, 0.25))
ax1.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 1.25, 0.25)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_xlabel(r'cos(θ$_{DW}$)',{'fontsize':24,'fontweight':'normal'})

ax1.axes.set_yticks(ax1.axes.get_yticks())
ax1.axes.set_yticklabels(['{0:.3f}'.format(x) for x in ax1.axes.get_yticks()],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_ylabel(r'P(cos(θ$_{DW}$))',{'fontsize':25,'fontweight':'normal'})


ax2 = plt.subplot(212)
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,6])
    ,hist=False,bins=bin_90_180_01,kde_kws={'clip': (-1,0) , "linestyle":"-","linewidth":"3.0","color":"gold"}, label=r'$gol$')
sns.distplot(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,6])
    ,hist=False,bins=bin_90_180_01,kde_kws={'clip': (-1,0) , "linestyle":"--","linewidth":"3.0","color":"red"}, label=r'$gl$')
sns.distplot(np.cos(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,6])
    ,hist=False,bins=bin_90_180_01,kde_kws={'clip': (-1,0), "linestyle":":","linewidth":"3.0","color":"black"}, label=r'$bl$')

plt.setp(ax2.spines.values(), linewidth=4)
ax2.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.axes.tick_params(which='minor',size=6,width=2)
ax2.axes.tick_params(which='major',size=8,width=4)
ax2.grid(axis='both',which='minor', alpha=0.2)
ax2.grid(axis='both',which='major', alpha=0.6)
ax2.xaxis.labelpad = 15
ax2.yaxis.labelpad = 15

ax2.set_xlim(-1,1)
#ax1.set_ylim(-1,1)

ax2.axes.set_xticks(np.arange(-1, 1.25, 0.25))
ax2.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 1.25, 0.25)],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_xlabel(r'cos(θ$_{HH}$)',{'fontsize':24,'fontweight':'normal'})

ax2.axes.set_yticks(ax2.axes.get_yticks())
ax2.axes.set_yticklabels(['{0:.3f}'.format(x) for x in ax2.axes.get_yticks()],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_ylabel(r'P(cos(θ$_{HH}$))',{'fontsize':25,'fontweight':'normal'})

plt.savefig('figure_AS_cos_DW_HH_T_revPBED3', bbox_inches='tight', transparent = True, dpi = 300)
plt.clf()
plt.cla()
plt.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
plt.rcParams["figure.figsize"] = [20,12]
fig1, ax = plt.subplots(dpi=300)

ax1 = plt.subplot(221)
plt.hist2d(np.degrees(WA_revPBED3[(WA_revPBED3[:,4] == 1)][:,5]),
    np.degrees(WA_revPBED3[(WA_revPBED3[:,4] == 1)][:,6]),
    range=((0,180),(90,180)),bins=[36,18],cmap="Blues")
cb1 = plt.colorbar()
cb1.set_label("Occurences",size=20)
cb1.outline.set_linewidth(4)
cb1.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax1.spines.values(), linewidth=4)
#ax1.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.axes.tick_params(which='minor',size=6,width=2)
ax1.axes.tick_params(which='major',size=8,width=4)
ax1.grid(axis='both',which='minor', alpha=0.2)
ax1.grid(axis='both',which='major', alpha=0.6)
ax1.xaxis.labelpad = 15
ax1.yaxis.labelpad = 15

ax1.set_xlim(0,180)
ax1.set_ylim(90,180)
ax1.axes.set_xticks(np.arange(0,200,20))
ax1.axes.set_xticklabels(['{0:.0f}'.format(x) for x in np.arange(0,200,20)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax1.axes.set_yticks(np.arange(90,190,10))
ax1.axes.set_yticklabels(['{0:.0f}'.format(x) for x in np.arange(90,190,10)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax1.text(0.02, 0.92, "(revPBE-D3-All)", transform=ax1.transAxes, fontsize=14, fontweight='bold')


ax2 = plt.subplot(222)
plt.hist2d(np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,5]),
    np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,6]),
    range=((0,180),(90,180)),bins=[36,18],cmap="Blues")
cb2 = plt.colorbar()
cb2.set_label("Occurences",size=20)
cb2.outline.set_linewidth(4)
cb2.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax2.spines.values(), linewidth=4)
#ax2.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.axes.tick_params(which='minor',size=6,width=2)
ax2.axes.tick_params(which='major',size=8,width=4)
ax2.grid(axis='both',which='minor', alpha=0.2)
ax2.grid(axis='both',which='major', alpha=0.6)
ax2.xaxis.labelpad = 15
ax2.yaxis.labelpad = 15

ax2.set_xlim(0,180)
ax2.set_ylim(90,180)
ax2.axes.set_xticks(np.arange(0,200,20))
ax2.axes.set_xticklabels(['{0:.0f}'.format(x) for x in np.arange(0,200,20)],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax2.axes.set_yticks(np.arange(90,190,10))
ax2.axes.set_yticklabels(['{0:.0f}'.format(x) for x in np.arange(90,190,10)],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax2.text(0.02, 0.92, "(revPBE-D3-gol)", transform=ax2.transAxes, fontsize=14, fontweight='bold')


ax3 = plt.subplot(223)
plt.hist2d(np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,5]),
    np.degrees(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,6]),
    range=((0,180),(90,180)),bins=[36,18],cmap="Blues")
cb3 = plt.colorbar()
cb3.set_label("Occurences",size=20)
cb3.outline.set_linewidth(4)
cb3.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax3.spines.values(), linewidth=4)
#ax3.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax3.xaxis.set_minor_locator(AutoMinorLocator(2))
ax3.yaxis.set_minor_locator(AutoMinorLocator(2))
ax3.axes.tick_params(which='minor',size=6,width=2)
ax3.axes.tick_params(which='major',size=8,width=4)
ax3.grid(axis='both',which='minor', alpha=0.2)
ax3.grid(axis='both',which='major', alpha=0.6)
ax3.xaxis.labelpad = 15
ax3.yaxis.labelpad = 15

ax3.set_xlim(0,180)
ax3.set_ylim(90,180)
ax3.axes.set_xticks(np.arange(0,200,20))
ax3.axes.set_xticklabels(['{0:.0f}'.format(x) for x in np.arange(0,200,20)],{'fontsize':16,'fontweight':'bold'})
ax3.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax3.axes.set_yticks(np.arange(90,190,10))
ax3.axes.set_yticklabels(['{0:.0f}'.format(x) for x in np.arange(90,190,10)],{'fontsize':16,'fontweight':'bold'})
ax3.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax3.text(0.02, 0.92, "(revPBE-D3-gl)", transform=ax3.transAxes, fontsize=14, fontweight='bold')

ax4 = plt.subplot(224)
plt.hist2d(np.degrees(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,5]),
    np.degrees(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,6]),
    range=((0,180),(90,180)),bins=[36,18],cmap="Blues")
cb4 = plt.colorbar()
cb4.set_label("Occurences",size=20)
cb4.outline.set_linewidth(4)
cb4.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax4.spines.values(), linewidth=4)
#ax4.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax4.xaxis.set_minor_locator(AutoMinorLocator(2))
ax4.yaxis.set_minor_locator(AutoMinorLocator(2))
ax4.axes.tick_params(which='minor',size=6,width=2)
ax4.axes.tick_params(which='major',size=8,width=4)
ax4.grid(axis='both',which='minor', alpha=0.2)
ax4.grid(axis='both',which='major', alpha=0.6)
ax4.xaxis.labelpad = 15
ax4.yaxis.labelpad = 15

ax4.set_xlim(0,180)
ax4.set_ylim(90,180)
ax4.axes.set_xticks(np.arange(0,200,20))
ax4.axes.set_xticklabels(['{0:.0f}'.format(x) for x in np.arange(0,200,20)],{'fontsize':16,'fontweight':'bold'})
ax4.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax4.axes.set_yticks(np.arange(90,190,10))
ax4.axes.set_yticklabels(['{0:.0f}'.format(x) for x in np.arange(90,190,10)],{'fontsize':16,'fontweight':'bold'})
ax4.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax4.text(0.02, 0.92, "(revPBE-D3-bl)", transform=ax4.transAxes, fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figure_AS_DW-HH_T_revPBED3', bbox_inches='tight', transparent = True, dpi = 300)
plt.clf()
plt.cla()
plt.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
plt.rcParams["figure.figsize"] = [20,12]
fig1, ax = plt.subplots(dpi=300)

ax1 = plt.subplot(221)
plt.hist2d(np.cos(WA_revPBED3[(WA_revPBED3[:,4] == 1)][:,5]),
    np.cos(WA_revPBED3[(WA_revPBED3[:,4] == 1)][:,6]),
    range=((-1,1),(-1,0)),bins=[bin_0_180,bin_90_180],cmap="Blues")

cb1 = plt.colorbar()
cb1.set_label("Occurences",size=20)
cb1.outline.set_linewidth(4)
cb1.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax1.spines.values(), linewidth=4)
#ax1.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.axes.tick_params(which='minor',size=6,width=2)
ax1.axes.tick_params(which='major',size=8,width=4)
#ax1.grid(axis='both',which='minor', alpha=0.2)
#ax1.grid(axis='both',which='major', alpha=0.6)
ax1.xaxis.labelpad = 15
ax1.yaxis.labelpad = 15

ax1.set_xlim(-1,1)
ax1.set_ylim(-1,0)
ax1.axes.set_xticks(np.arange(-1, 1.2, 0.2))
ax1.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 1.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax1.axes.set_yticks(np.arange(-1, 0.2, 0.2))
ax1.axes.set_yticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 0.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax1.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax1.text(0.02, 0.02, "(revPBE-D3-All)", transform=ax1.transAxes, fontsize=14, fontweight='bold')


ax2 = plt.subplot(222)
plt.hist2d(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,5]),
    np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 1) | (WA_revPBED3[:,1] == 1) | (WA_revPBED3[:,2] == 1) ) & (WA_revPBED3[:,3] == 1) )][:,6]),
    range=((-1,1),(-1,0)),bins=[bin_0_180,bin_90_180],cmap="Blues")
cb2 = plt.colorbar()
cb2.set_label("Occurences",size=20)
cb2.outline.set_linewidth(4)
cb2.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax2.spines.values(), linewidth=4)
#ax2.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.axes.tick_params(which='minor',size=6,width=2)
ax2.axes.tick_params(which='major',size=8,width=4)
#ax2.grid(axis='both',which='minor', alpha=0.2)
#ax2.grid(axis='both',which='major', alpha=0.6)
ax2.xaxis.labelpad = 15
ax2.yaxis.labelpad = 15

ax2.set_xlim(-1,1)
ax2.set_ylim(-1,0)
ax2.axes.set_xticks(np.arange(-1, 1.2, 0.2))
ax2.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 1.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax2.axes.set_yticks(np.arange(-1, 0.2, 0.2))
ax2.axes.set_yticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 0.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax2.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax2.text(0.02, 0.02, "(revPBE-D3-gol)", transform=ax2.transAxes, fontsize=14, fontweight='bold')


ax3 = plt.subplot(223)
plt.hist2d(np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,5]),
    np.cos(WA_revPBED3[( ( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) ) & (WA_revPBED3[:,3] == 1) )][:,6]),
    range=((-1,1),(-1,0)),bins=[bin_0_180,bin_90_180],cmap="Blues")
cb3 = plt.colorbar()
cb3.set_label("Occurences",size=20)
cb3.outline.set_linewidth(4)
cb3.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax3.spines.values(), linewidth=4)
#ax3.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax3.xaxis.set_minor_locator(AutoMinorLocator(2))
ax3.yaxis.set_minor_locator(AutoMinorLocator(2))
ax3.axes.tick_params(which='minor',size=6,width=2)
ax3.axes.tick_params(which='major',size=8,width=4)
#ax3.grid(axis='both',which='minor', alpha=0.2)
#ax3.grid(axis='both',which='major', alpha=0.6)
ax3.xaxis.labelpad = 15
ax3.yaxis.labelpad = 15

ax3.set_xlim(-1,1)
ax3.set_ylim(-1,0)
ax3.axes.set_xticks(np.arange(-1, 1.2, 0.2))
ax3.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 1.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax3.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax3.axes.set_yticks(np.arange(-1, 0.2, 0.2))
ax3.axes.set_yticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 0.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax3.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax3.text(0.02, 0.02, "(revPBE-D3-gl)", transform=ax3.transAxes, fontsize=14, fontweight='bold')

ax4 = plt.subplot(224)
plt.hist2d(np.cos(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,5]),
    np.cos(WA_revPBED3[( (WA_revPBED3[:,0] == 0) & (WA_revPBED3[:,1] == 0) & (WA_revPBED3[:,2] == 0) & (WA_revPBED3[:,3] == 0) & (WA_revPBED3[:,4] == 1) )][:,6]),
    range=((-1,1),(-1,0)),bins=[bin_0_180,bin_90_180],cmap="Blues")
cb4 = plt.colorbar()
cb4.set_label("Occurences",size=20)
cb4.outline.set_linewidth(4)
cb4.ax.tick_params(which='major',size=6,width=4,labelsize=16)

plt.setp(ax4.spines.values(), linewidth=4)
#ax4.legend(fontsize=20,framealpha=0.3,fancybox=True)
ax4.xaxis.set_minor_locator(AutoMinorLocator(2))
ax4.yaxis.set_minor_locator(AutoMinorLocator(2))
ax4.axes.tick_params(which='minor',size=6,width=2)
ax4.axes.tick_params(which='major',size=8,width=4)
#ax4.grid(axis='both',which='minor', alpha=0.2)
#ax4.grid(axis='both',which='major', alpha=0.6)
ax4.xaxis.labelpad = 15
ax4.yaxis.labelpad = 15

ax4.set_xlim(-1,1)
ax4.set_ylim(-1,0)
ax4.axes.set_xticks(np.arange(-1, 1.2, 0.2))
ax4.axes.set_xticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 1.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax4.axes.set_xlabel(r'θ$_{DW}$',{'fontsize':24,'fontweight':'normal'})

ax4.axes.set_yticks(np.arange(-1, 0.2, 0.2))
ax4.axes.set_yticklabels(['{0:.2f}'.format(x) for x in np.arange(-1, 0.2, 0.2)],{'fontsize':16,'fontweight':'bold'})
ax4.axes.set_ylabel(r'P(θ$_{HH}$)',{'fontsize':25,'fontweight':'normal'})
ax4.text(0.02, 0.02, "(revPBE-D3-bl)", transform=ax4.transAxes, fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figure_AS_cos_DW-HH_T_revPBED3', bbox_inches='tight', transparent = True, dpi = 300)
plt.clf()
plt.cla()
plt.close()

# ---------------------------------------------------------------------------------------------------------------------------------------
quit()
