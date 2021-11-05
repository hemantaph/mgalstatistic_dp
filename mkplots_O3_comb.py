import numpy as np
from pylab import plt
import sys
import snronpairs_quad_dbl as sn
import contfunc as cf

## O3 
det="O3"

fig = plt.figure()
lab=np.array(["O3, $\Delta\phi$=0", "O3, $\Delta\phi$=90"])


def makeplot(kk,mag_comb, tdel_comb, gwm,gwt,gwmerr,magu,tdu,lab):
    
    ax1=fig.add_subplot(2,2,kk)

    #######
    
    xxl,yyl,ffl=cf.contfunc(tdel_comb,mag_comb)
    
     
    ax1.plot(0,0,label=lab,color="g")
    levels_l=cf.contlevs(ffl) 
    cset=ax1.contour(np.rot90(ffl), levels_l, colors='g', origin='upper', extent=extent_im)
    
    
    xxu,yyu,ffu=cf.contfunc(tdu,magu)
    levels=cf.contlevs(ffu) 
    cset=ax1.contour(np.rot90(ffu), levels, colors='k', origin='upper', extent=extent_im)
     
    levels_n=cf.contlevs(ffu, returnfirst=True) 
    cfset = ax1.contourf(xxu, yyu,  ffu, levels_n,cmap='YlOrBr')
        
   
    if(kk==1):
        ax1.errorbar(np.log10(gwt),np.log10(gwm),color="blue",yerr=gwmerr/gwm/np.log(10.),fmt="o",alpha=0.7) 
        ax1.set_ylabel("Log (Relative magnification)")
    if(kk==2): 
        ax1.tick_params(left = False, labelleft = False )
        ax1.errorbar(np.log10(gwt),np.log10(gwm),label="O3-GW Events",color="blue",yerr=gwmerr/gwm/np.log(10.),alpha=0.7,fmt="o") 
    
    ax1.set_xlabel("Log (Time delay) (days)")
    plt.legend(loc="lower left",fontsize=11)



## Extract the magnification and time delay distribution for O3 - quads and doubles
mag31,tdel31,mag32,tdel32,mag41,tdel41,mag42,tdel42,mag21,tdel21,mag43,tdel43,dbmg21,dbtd21= sn.getsnr_forpairs(det)

print(mag31.size,mag32.size,mag41.size,mag42.size, mag21.size,mag43.size, dbmg21.size)

## phase diff of 0
mag1=np.concatenate((mag21,mag43))
tdel1=np.concatenate((tdel21,tdel43))
idx1=(tdel1>1e-3) & (mag1>5e-2)

## phase diff of 90
mag2=np.concatenate([mag31,mag32,mag41,mag42,dbmg21])
tdel2=np.concatenate([tdel31,tdel32,tdel41,tdel42,dbtd21])
idx2=(tdel2>1e-3) & (mag2>5e-2)


#### Unlensed pairs
tdu,magu=np.loadtxt("unlensedpairs/unlensedpairs_tdmag_%s.txt"%(det),unpack=1) 
idxu=(tdu>1e-3) & (magu>1e-3)

## GW events difference of gps times; in seconds
gwtdel0,gwrmag0,gwmerr0=np.loadtxt("gwevents_O3/gwevents_phdf0.txt",usecols=(2,3,4),unpack=1)
gwtdel90,gwrmag90,gwmerr90=np.loadtxt("gwevents_O3/gwevents_phdf90.txt",usecols=(2,3,4),unpack=1)

xmin =np.log10(1e-2) 
xmax =np.log10(5e2)  
ymin =np.log10(1e-2)    
ymax =np.log10(1e2)  
extent_im=[xmin, xmax, ymin, ymax]


makeplot(1,mag1[idx1],tdel1[idx1],gwrmag0,gwtdel0/3600/24,gwmerr0,magu[idxu],tdu[idxu],lab[0])
makeplot(2,mag2[idx2],tdel2[idx2],gwrmag90,gwtdel90/3600/24,gwmerr90,magu[idxu],tdu[idxu],lab[1])

plt.subplots_adjust(hspace=0.,wspace=0.02)
plt.savefig("O3_comb.pdf")

