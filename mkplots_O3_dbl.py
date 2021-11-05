import numpy as np
from pylab import plt
import scipy.stats as st
import sys
from scipy.interpolate import interp1d
import contfunc as cf

det="O3"


## Double SIE
dbid,dbimx,dbimy,dbmag,dbtdel,dblparity,id2=np.loadtxt("out_fincat_%s/dblall_imgprop.txt"%(det),unpack=True)

fig = plt.figure()
for ii in range(0,dbtdel.size,4):
    indx1=[ii,ii+1,ii+2,ii+3]

indx1=np.arange(1,dbtdel.size+1,1)
dbetdel2=dbtdel[indx1%4==2]

dbemag1=dbmag[indx1%4==1]
dbermag2= -dbmag[indx1%4==2]/dbemag1

### Double SIS
#dbid,dbimx,dbimy,dbmag,dbtdel,dblparity,id2=np.loadtxt("out_fincat_%ssis/dblall_imgprop.txt"%(det),unpack=True)
#
#for ii in range(0,dbtdel.size,4):
#    indx1=[ii,ii+1,ii+2,ii+3]
#
#indx1=np.arange(1,dbtdel.size+1,1)
#dbstdel2=dbtdel[indx1%4==2]
#
#dbsmag1=dbmag[indx1%4==1]
#dbsrmag2= -dbmag[indx1%4==2]/dbsmag1
#


def contplot(kk,xdata,ydata,gwt,gwm,gwmerr,xdata_un,ydata_un,lab,extent_im):
   # Extract x and y
    ax1=fig.add_subplot(2,2,kk)
    xxl,yyl,ffl=cf.contfunc(xdata,ydata)
 
    ax1.plot(0,0,label=lab,color="g")
    ##cset = ax1.contour(xxl, yyl, ffl, colors='r') ##linewidths=1.0)

    levels_l=cf.contlevs(ffl) 
    cset=ax1.contour(np.rot90(ffl), levels_l, colors='g', origin='upper', extent=extent_im)

    xxu,yyu,ffu=cf.contfunc(xdata_un,ydata_un)
    levels=cf.contlevs(ffu) 
    cset=ax1.contour(np.rot90(ffu), levels, colors='k', origin='upper', extent=extent_im)
      
    levels_n=cf.contlevs(ffu, returnfirst=True) 
    cfset = ax1.contourf(xxu, yyu,  ffu, levels_n,cmap='YlOrBr')
   #ax1.imshow(np.rot90(ffu), cmap='YlGnBu', extent=extent_im,alpha=0.1) #aspect=(extent_im[0]-extent_im[1])/(extent_im[2]-extent_im[3]))
    if(kk==1):
        ax1.errorbar(np.log10(gwt),np.log10(gwm),label="O3-GW Events",color="blue",yerr=(gwmerr/gwm/np.log(10.)),alpha=0.7,fmt="o") 
    else:
        ax1.errorbar(np.log10(gwt),np.log10(gwm),color="blue",yerr=(gwmerr/gwm/np.log(10.)),alpha=0.7,fmt="o") 

    #ax1.clabel(cset, inline=1, fontsize=10)
   #cset = ax1.contour(xxu, yyu, ffu, colors='k')
   
    ax1.set_xlabel("Log (Time delay) (days)")
    if(kk==1):
        ax1.set_ylabel("Log (Relative magnification)")
    if(kk==2):
   #     ax1.text(-3.2,0.3,"Log (Relative magnification)",fontsize=11,rotation=90)
        ax1.tick_params(left = False, labelleft = False )
      # ax1.tick_params(labelcolor='w', left=False)

    ax1.set_xlim(extent_im[0],extent_im[1])
    ax1.set_ylim(extent_im[2],extent_im[3])
    plt.legend(loc="lower left") 



#### Unlensed pairs
tdu,magu=np.loadtxt("unlensedpairs/unlensedpairs_tdmag_%s.txt"%(det),unpack=1) 
idxu=(tdu>1e-3) & (magu>1e-3)

## GW events difference of gps times; in seconds
##gwtdel0,gwrmag0=np.loadtxt("gwevents_O3/gwevents_phdf0.txt",usecols=(2,3),unpack=1)
gwtdel90,gwrmag90,gwrmag90err=np.loadtxt("gwevents_O3/gwevents_phdf90.txt",usecols=(2,3,4),unpack=1)

######

xmin =np.log10(1e-2) 
xmax =np.log10(5e2)  
ymin =np.log10(1e-2)    
ymax =np.log10(1e2)  
extent_im=[xmin, xmax, ymin, ymax]


### Plotting

contplot(1,dbetdel2,dbermag2,gwtdel90/3600/24,gwrmag90,gwrmag90err,tdu[idxu],magu[idxu],"O3-dbl-SIE",extent_im)
#contplot(2,dbstdel2,dbsrmag2,gwtdel90/3600/24,gwrmag90,gwrmag90err,tdu[idxu],magu[idxu],"O3-dbl-SIS",extent_im)

plt.subplots_adjust(hspace=0.,wspace=0.02)
plt.savefig("O3_dbl.pdf")

