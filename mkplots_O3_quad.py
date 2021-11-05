import numpy as np
from pylab import plt
import sys
import snronpairs_quad_dbl as sn
import contfunc as cf

## O3 
det="O3"


def contplot(kk,xdata,ydata,gwt,gwm,gwmerr,xdata_un,ydata_un,lab,extent_im):
   # Extract x and y
    ax1=fig.add_subplot(3,2,kk)
    xxl,yyl,ffl=cf.contfunc(xdata,ydata)

    ax1.plot(0,0,label=lab,color="g")
    levels_l=cf.contlevs(ffl) 
    cset=ax1.contour(np.rot90(ffl), levels_l, colors='g', origin='upper', extent=extent_im)

    xxu,yyu,ffu=cf.contfunc(xdata_un,ydata_un)
    levels=cf.contlevs(ffu) 
    cset=ax1.contour(np.rot90(ffu), levels, colors='k', origin='upper', extent=extent_im)

    levels_n=cf.contlevs(ffu, returnfirst=True) 
    cfset = ax1.contourf(xxu, yyu,  ffu, levels_n,cmap='YlOrBr')
    if(kk==5):
        ax1.errorbar(np.log10(gwt),np.log10(gwm),label="O3-GW Events",color="blue",yerr=gwmerr/gwm/np.log(10.),alpha=0.7,fmt="o") 
    else:
        ax1.errorbar(np.log10(gwt),np.log10(gwm),color="blue",yerr=gwmerr/gwm/np.log(10.),alpha=0.7,fmt="o") 
   
    if(kk<5): 
        ax1.tick_params(labelbottom = False, bottom = False)
    if(kk%2==0 ): 
        ax1.tick_params(left = False, labelleft = False )
    if(kk==3):
        ax1.set_ylabel("Log (Relative magnification)")
    ax1.set_xlim(extent_im[0],extent_im[1])
    ax1.set_ylim(extent_im[2],extent_im[3])
    plt.legend(loc="lower left") 


### Lensed pairs
magr31,tdel31,magr32,tdel32,magr41,tdel41,magr42,tdel42,magr21,tdel21,magr43,tdel43,dbmg21,dbtd21=sn.getsnr_forpairs(det)

print(magr31.size,magr32.size,magr41.size,magr42.size,magr21.size,magr43.size,dbmg21.size)


#### Unlensed pairs
tdu,magu=np.loadtxt("unlensedpairs/unlensedpairs_tdmag_%s.txt"%(det),unpack=1) 
idxu=(tdu>1e-3) & (magu>1e-3)

## GW events difference of gps times; in seconds
gwtdel0,gwrmag0,gwrmag0err=np.loadtxt("gwevents_O3/gwevents_phdf0.txt",usecols=(2,3,4),unpack=1)
gwtdel90,gwrmag90,gwrmag90err=np.loadtxt("gwevents_O3/gwevents_phdf90.txt",usecols=(2,3,4),unpack=1)

#####
fig = plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_xlabel("Log (Time delay) (days)")

## Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
xmin =np.log10(1e-2) 
xmax =np.log10(5e2)  
ymin =np.log10(1e-2)    
ymax =np.log10(1e2)  
extent_im=[xmin, xmax, ymin, ymax]


## Plotting

contplot(1,tdel31,magr31,gwtdel90/3600./24,gwrmag90,gwrmag90err,tdu[idxu],magu[idxu],"O3-31",extent_im)
contplot(2,tdel32,magr32,gwtdel90/3600./24,gwrmag90,gwrmag90err,tdu[idxu],magu[idxu],"O3-32",extent_im)
contplot(3,tdel41,magr41,gwtdel90/3600./24,gwrmag90,gwrmag90err,tdu[idxu],magu[idxu],"O3-41",extent_im)
contplot(4,tdel42,magr42,gwtdel90/3600./24,gwrmag90,gwrmag90err,tdu[idxu],magu[idxu],"O3-42",extent_im)
contplot(5,tdel21,magr21, gwtdel0/3600./24,gwrmag0 ,gwrmag0err,tdu[idxu],magu[idxu],"O3-21",extent_im)
contplot(6,tdel43,magr43, gwtdel0/3600./24,gwrmag0 ,gwrmag0err,tdu[idxu],magu[idxu],"O3-43",extent_im)

plt.subplots_adjust(hspace=0.05,wspace=0.01)
plt.savefig("O3_quad.pdf")

