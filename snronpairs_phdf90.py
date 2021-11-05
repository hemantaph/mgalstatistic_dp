import numpy as np
import sys
import pylab as plt
import scipy.stats as st
import contfunc as cf

## For phase diff. of 0 i.e pairs 4-3, 2-1

det=np.array(["AL", "Ap","CE","ET" ])
SNRthresh=8

fig = plt.figure()

def getsnr_forpairs(det):
   lid,uSNR=np.loadtxt("out_fincat_%s/all_lensprop.txt"%(det),usecols=(0,6),unpack=True)
   
   qdid,qdmag,qdtdel=np.loadtxt("out_fincat_%s/quadall_imgprop.txt"%(det),usecols=(0,3,4),unpack=True)
   dbid,dbmag,dbtdel=np.loadtxt("out_fincat_%s/dblall_imgprop.txt"%(det),usecols=(0,3,4),unpack=True)
       
   
   mag_rat = []
   tdel_rat = []
   
   
   for ii in range(lid.size):
       ## get indices of the matched IDs 
       indx = qdid==np.int(lid[ii])
       indx1 = dbid==np.int(lid[ii])

       if(np.sum(indx)==0):
          ## double 
          dbtdel2= dbtdel[indx1][1] 
          dbmag1 = dbmag[indx1][0]         
          dbmag2 = dbmag[indx1][1]
       else:
          ## quad
          qdtdel2=qdtdel[indx][1] 
          qdtdel3=qdtdel[indx][2] 
          qdtdel4=qdtdel[indx][3] 
        
          qdmag1 = qdmag[indx][0]         
          qdmag2= qdmag[indx][1]  
          qdmag3= qdmag[indx][2]  
          qdmag4= qdmag[indx][3]  
      

          ########
          qdmin42 = min(np.absolute(qdmag4),np.absolute(qdmag2)) 
          qdmin32 = min(np.absolute(qdmag3),np.absolute(qdmag2)) 
          qdmin41 = min(np.absolute(qdmag4),np.absolute(qdmag1)) 
          qdmin31 = min(np.absolute(qdmag3),np.absolute(qdmag1)) 
          dbmin21 = min(np.absolute(dbmag2),np.absolute(dbmag1)) 


          if( np.sqrt(qdmin42)*uSNR[ii] > SNRthresh ):
             mag_rat.append( -qdmag4/qdmag2 )
             tdel_rat.append( qdtdel4-qdtdel2 )
     
          if( np.sqrt(qdmin32)*uSNR[ii] > SNRthresh ):
             mag_rat.append( -qdmag3/qdmag2 )
             tdel_rat.append( qdtdel3-qdtdel2 )
      
          if( np.sqrt(qdmin41)*uSNR[ii] > SNRthresh ):
             mag_rat.append( -qdmag4/qdmag1 )
             tdel_rat.append( qdtdel4 )
      
          if( np.sqrt(qdmin31)*uSNR[ii] > SNRthresh ):
             mag_rat.append( -qdmag3/qdmag1 )
             tdel_rat.append( qdtdel3 )
   
          if( np.sqrt(dbmin21)*uSNR[ii] > SNRthresh ):
             mag_rat.append(-dbmag2/dbmag1 )
             tdel_rat.append( dbtdel2 )
          
   mag_rat=np.array(mag_rat)
   tdel_rat=np.array(tdel_rat)
   
   return mag_rat,tdel_rat



if __name__== "__main__":

   for kk in range(det.size):
   
      mag_rat,tdel_rat=getsnr_forpairs(det[kk])
   
      idx=(tdel_rat>1e-3) & (mag_rat>5e-2)
   
      xmin =np.log10(1e-2) 
      xmax =np.log10(5e2)  
      ymin =np.log10(1e-2)    
      ymax =np.log10(1e2)  
      extent_im=[xmin, xmax, ymin, ymax]
      
      ax1=fig.add_subplot(2,2,kk+1)
      ax1.plot(0,0,label=det[kk],color="g")
      xxl,yyl,ffl=cf.contfunc(tdel_rat[idx],mag_rat[idx])
      ax1.plot(0,0,label=det[kk],color="g")
      levels_l=cf.contlevs(ffl) 
      cset=ax1.contour(np.rot90(ffl), levels_l, colors='g', origin='upper', extent=extent_im)
      plt.legend(loc="lower left") 
   
   plt.subplots_adjust(hspace=0.,wspace=0.02)
   plt.savefig("sample_phdf90.pdf")
   
