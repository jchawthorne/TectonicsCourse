# a collection of functions useful for the lunar gravity lab

import numpy as np
import pyshtools as pysh
import matplotlib.pyplot as plt

def filter_sph(data,llow=0,lhigh=None,n=2):
    """
    filter the spherical harmonics to reduce signals at 
    degrees smaller than llow or higher than lhigh
    :param        data: a pysh object with coefficients
    :param        llow: the lowest degree of interest
    :param       lhigh: the highest degree of interest
    :param           n: amplitude reduced by a factor of (l-lref)^(-n)
    """

    
    # initial data and degrees
    cfs=data.coeffs
    l=data.degrees().astype(float)

    # pick references
    if lhigh is None:
        lhigh=np.max(l)

    # amplitude to multiply by
    amp=np.ones(l.size,dtype=float)
    ilow=l<llow
    amp[ilow]=np.power(llow/l[ilow],-n)
    ihigh=l>lhigh
    amp[ihigh]=np.power(l[ihigh]/lhigh,-n)

    amp=amp.reshape([1,amp.size,1])

    # amply filter
    cfs=np.multiply(cfs,amp)

    # and make new SH object
    dataf=pysh.SHCoeffs.from_array(cfs,normalization=data.normalization)

    return dataf
    
