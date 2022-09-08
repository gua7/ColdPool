#!/usr/bin/env python
# coding: utf-8

# In[1]:

from scipy.ndimage import gaussian_filter as Gaussian
import netCDF4 as nc
import xarray as xr
import numpy as np
from clusterpy import cluster
import numpy.ma as ma
import scipy as sp


# In[6]:


fn ='the0001g5.nc'
#fn2 ='qt.nc'

ds = nc.Dataset(fn)
ds2 = nc.Dataset('/net/shared-public/les-precip/192x192_final/reference/r1/output/cross/0001/combined/crossxy.0001.qtxy.100.nc')
ds3 = nc.Dataset('/net/shared-public/les-precip/192x192_final/reference/r1/output/cross/0001/combined/crossxy.0001.uxy.100.nc')
ds4 = nc.Dataset('/net/shared-public/les-precip/192x192_final/reference/r1/output/cross/0001/combined/crossxy.0001.vxy.100.nc')
#qt = qtf.variables['qtxy'][:]


# In[28]:


times=ds['time'][:]



# In[10]:


# In[28]:

numbers = np.zeros((len(times)))
mean_the = np.nan*np.ones((len(times),1000))
min_the = np.nan*np.ones((len(times),1000))
size = np.nan*np.ones((len(times),1000))
mean_qt = np.nan*np.ones((len(times),1000))
mean_w = np.nan*np.ones((len(times),1000))
infile = fn
ds_disk = xr.open_dataset(infile, decode_times=False)

for t in range(len(times)):
    
    the = ds['thexy'][t,:,:]
    qt = ds2['qtxy'][t,:,:]
    u = ds3['uxy'][t,:,:]
    v = ds4['vxy'][t,:,:]
    wind = np.sqrt(u**2+v**2)
    wind = wind - np.mean(wind)
    wind = Gaussian(wind,sigma=5,mode='wrap')
    qt = Gaussian(qt,sigma=5,mode='wrap')
    thres = -1
    tslice = ds_disk.isel(time=t, x=range(ds_disk.sizes['x']), y=range(ds_disk.sizes['y'])).to_array()
    cldata = xr.where(tslice > -999999999, -1, -1)
    cells = cluster.ClusterArray(-tslice,-thres).get_clusterarray()
    cells = cells[0,:,:]  
    cells_masked = ma.masked_less(cells,0)
    the_masked = ma.masked_array(the,mask = cells_masked.mask)
    qt_masked = ma.masked_array(qt,mask = cells_masked.mask)
    wind_masked = ma.masked_array(wind,mask = cells_masked.mask)
    cp_num = int(np.max(cells+1))
    #numbers[t] = cp_num
    for i in range(cp_num):
        cell = ma.masked_not_equal(cells,i)
        cp = ma.masked_array(the_masked,mask=cell.mask)
        humi = ma.masked_array(qt_masked,mask=cell.mask)
        w = ma.masked_array(wind_masked,mask=cell.mask)
        mean_the[t,i] = np.mean(cp)
        min_the[t,i] = np.min(cp)
        size[t,i] = np.sum(cells_masked==i)
        mean_qt[t,i] = np.mean(humi)
        mean_w[t,i] = np.mean(w)
        

# In[2]:


np.savetxt('mean_the',mean_the)
np.savetxt('min_the',min_the)
np.savetxt('size_',size)
np.savetxt('mean_wind',mean_w)
np.savetxt('mean_qt',mean_qt)



ds.close()
ds2.close()




