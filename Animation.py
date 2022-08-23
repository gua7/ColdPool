## Thanks to Pouriya Alinaghi for the code ##
## Code of -4 K condition

import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.ndimage import gaussian_filter
import matplotlib.animation as animation
from skimage.transform import rescale
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches


ds0 = nc.Dataset('nc/01ref/crossxy.0001.thlxy.100.nc') 
ds3 = nc.Dataset('nc/-4/crossxy.0003.wxy.100.nc') #gust front
ds5 = nc.Dataset('nc/-4/cape.lwp.100.nc') #cloud

ds1 = nc.Dataset('nc/-4/crossxy.0001.thlxy.100.nc')
ds2 = nc.Dataset('nc/-4/crossxy.0001.qtxy.100.nc')

ds4 = nc.Dataset('nc/-4/cape.surfprec.100.nc') #precipitation

ds6 = nc.Dataset('nc/-4/crossxy.0001.uxy.100.nc')
ds7 = nc.Dataset('nc/-4/crossxy.0001.vxy.100.nc') #horizontal wind



def update_plot(t):
    ax1.clear()
    ax2.clear()

    CL  = ds5['lwp'][t,:,:]
    CL = gaussian_filter(CL, sigma=2)
    thrCL= 0.05
    CL  = 1 * (CL > thrCL) + 0 * (CL <= thrCL)
    ax1.contourf(CL, colors = 'gray',levels = [0.99,1],alpha = 0.3)
 
    thl = ds1['thlxy'][t,:,:]
    qt = ds2['qtxy'][t,:,:]
    C = thl + 2500 * qt
    Cm = np.mean(C)
    C  = C - Cm * np.ones(np.shape(C))
    C = gaussian_filter(C, sigma=2)
    thrC = 1
    C  = 1 * (C <= -thrC) + 0 * (C > -thrC)
    ax1.contourf(C, colors = 'gold',levels = [0.99,1], alpha = 0.75)
    #C = np.round_(C , decimals = 3)
    #ax2.contour(C, colors = 'black', levels=[1])
    thl = ds1['thlxy'][t,:,:]
    qt = ds2['qtxy'][t,:,:]
    C2 = thl + 2500 * qt 
    C2m = np.mean(C2)
    C2  = C2 - C2m * np.ones(np.shape(C2))
    C2 = gaussian_filter(C2, sigma=2)
    thrC2= 2
    C2  = 1 * (C2 <= -thrC2) + 0 * (C2 > -thrC2)
    ax1.contourf(C2, colors = 'orange',levels = [0.99,1],alpha=0.75)
    
    
    w = ds3['wxy'][t,:,:]
    w = gaussian_filter(w, sigma=2)
    thrw = 0.1
    w1 = 1 * (w >= thrw) + 0 * (w < thrw)
    w2 = 1 * (w <= -thrw) + 0 * (w > -thrw)
    ax1.contourf(w1, colors = 'olive',levels = [0.99,1],alpha=0.7)
   

    R  = ds4['surfprec'][t,:,:]
    thrR1= 1e-5
    R1  = 1 * (R > thrR1) + 0 * (R <= thrR1)
    ax1.contourf(R1, colors = 'lightskyblue',levels = [0.99,1],alpha = 0.3)
    R2 = ds4['surfprec'][t,:,:]
    thrR2= 1e-4
    R2 = 1 * (R > thrR2) + 0 * (R <= thrR2)
    ax1.contourf(R2, colors = 'deepskyblue',levels = [0.99,1],alpha = 0.7)
  
    ax1.set_title(f'Cold pools with cloud, rain and gust front,time = {t + 1}')


    cloudcol = mpatches.Patch(color='gray', label='Clouds')
    raincol = mpatches.Patch(color='lightskyblue', label='Rain(1e-5)')
    raincol2 = mpatches.Patch(color='deepskyblue', label='Rain(1e-4)')
    cdpcol = mpatches.Patch(color='gold', label='Cold Pool(-1K)')
    cdp2col = mpatches.Patch(color='orange', label='Cold Pool(-2K)')
    wcol = mpatches.Patch(color='olive', label='Gust Front')
    wcol2 = mpatches.Patch(color='lime', label='Downward Wind')
    ax1.legend(handles=[cloudcol,raincol,raincol2,cdpcol,cdp2col,wcol])

    u = ds6['uxy'][t,:,:]
    u = gaussian_filter(u, sigma=2)
    

    v = ds7['vxy'][t,:,:]
    v = gaussian_filter(v, sigma=2)
   # u = rescale(u, 0.04)
   # v= rescale(v, 0.04)

    wind = (u**2 + v**2) ** 0.5
    Wm = np.mean(wind)
    wind = wind-Wm
    wind[wind > 5] = 5
    wind[wind < -5] = -5
    levels = np.arange(-5, 5.1, 0.5)
 
    ax2.contourf(wind, levels=levels, cmap = 'bwr')
    im = ax2.contourf(wind, levels=levels, cmap = 'bwr')
    ax2.contourf(w1, colors = 'olive',levels = [0.99,1],alpha=0.7)
    ax2.contourf(w2, colors = 'lime',levels = [0.99,1],alpha=0.7)
    #ax2.quiver(u, v)
    ax2.set_title(f'Horizontal Wind Anomaly with Gust Front and Downward Wind,time = {t + 1}')
    fig.colorbar(im,cax=cax)
    ax2.legend(handles=[wcol,wcol2])
    print(t)


fig = plt.figure()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize =(24,10))
div = make_axes_locatable(ax2)
cax = div.append_axes('right', '5%', '5%')




anim = FuncAnimation(fig, update_plot,save_count=1800)

f = r"m4.mp4" 
writervideo = animation.FFMpegWriter(fps=30) 
anim.save(f, writer=writervideo,dpi=200)
