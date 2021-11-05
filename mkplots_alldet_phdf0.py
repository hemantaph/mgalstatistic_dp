import numpy as np
from pylab import plt
import scipy.stats as st
import sys
import snronpairs_phdf0 as sn
from scipy.interpolate import interp1d
import contfunc as cf

## O3, Ap, CE or ET
det=np.array(["AL", "Ap", "CE","ET" ])
#det=sys.argv[1]

fig = plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_ylabel("Log (Relative magnification)")
ax.set_xlabel("Log (Time delay) (days)")

## Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

def make_plot(kk,tdel_comb,mag_comb,tdu,magu,extent_im):    
    
    ax1=fig.add_subplot(2,2,kk+1)
     
    ##cfset = ax1.contourf(xx, yy, f, cmap='YlGnBu')
    ##cset = ax1.contour(xx, yy, f, colors='g')
     
    ###cfset = ax1.contourf(xx, yy, fu, cmap='YlOrBr')
    ###ax1.imshow(np.rot90(fu), cmap='YlGnBu', extent=[xmin, xmax, ymin, ymax],alpha=0.1,aspect=(xmin-xmax)/(ymin-ymax))
    ###cset = ax1.contour(xx, yy, fu, colors='k')
 
    xxl,yyl,ffl=cf.contfunc(tdel_comb,mag_comb)
    ax1.plot(0,0,label=det[kk],color="g")
    levels_l=cf.contlevs(ffl) 
    cset=ax1.contour(np.rot90(ffl), levels_l, colors='g', origin='upper', extent=extent_im)

    xxu,yyu,ffu=cf.contfunc(tdu,magu)
    levels=cf.contlevs(ffu) 
    cset=ax1.contour(np.rot90(ffu), levels, colors='k', origin='upper', extent=extent_im)
      
    levels_n=cf.contlevs(ffu, returnfirst=True) 
    cfset = ax1.contourf(xxu, yyu,  ffu, levels_n,cmap='YlOrBr')
    
    if((kk+1)<2): 
        ax1.tick_params(labelbottom = False, bottom = False,labelleft=True)
    if((kk+1)%2==0 ): 
        ax1.tick_params(left = False, labelleft = False )
    plt.legend(loc="lower left") 


xmin =np.log10(1e-2) 
xmax =np.log10(5e2)  
ymin =np.log10(1e-2)    
ymax =np.log10(1e2)  
extent_im=[xmin, xmax, ymin, ymax]

for kk in range(det.size):
    ### lensed pairs
    mag_comb,tdel_comb=sn.getsnr_forpairs(det[kk]) 
    idx=(tdel_comb>1e-3) & (mag_comb>1e-3)
        
    #### Unlensed pairs
    tdu,magu=np.loadtxt("unlensedpairs/unlensedpairs_tdmag_%s.txt"%(det[kk]),unpack=1) 
    idxu=(tdu>1e-3) & (magu>1e-3)
        
    make_plot(kk,tdel_comb[idx],mag_comb[idx],tdu[idxu],magu[idxu],extent_im)


plt.subplots_adjust(hspace=0.,wspace=0.02)
plt.savefig("alldet_phdf0.pdf")
