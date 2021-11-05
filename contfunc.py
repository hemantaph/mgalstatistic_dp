from scipy.interpolate import interp1d
import scipy.stats as st
import numpy as np


def contfunc(xdat,ydat):
    # Extract x and y
    xu=np.log10(xdat)
    yu=np.log10(ydat)

    xmin = np.log10(1e-2)
    xmax = np.log10(5e2) 
    ymin = np.log10(1e-2)
    ymax = np.log10(1e2) 
    
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([xu, yu])
    kernel = st.gaussian_kde(values)
    ff = np.reshape(kernel(positions).T, xx.shape)
    return xx,yy,ff


def contlevs(ff, returnfirst=False):
    zsort=-np.sort(-ff.flatten())
    
    cumz = np.cumsum(zsort)/np.sum(zsort)
    spl = interp1d(cumz, zsort)
    level1 = spl(0.10)
    level2 = spl(0.40)
    level3 = spl(0.68)
    level4 = spl(0.95)
    if returnfirst:
        level5 = zsort[0]
        levels= np.array([level5, level1,level2, level3, level4])[::-1]
    else:
        levels= np.array([level1,level2, level3, level4])[::-1]
    return levels
