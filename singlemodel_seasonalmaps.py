#########################################################################
#
#  OPTICAL PROPERTIES: SEASONAL MAPS
#
#  For a given simulation: maps of AOD, AAOD, AE
#  Sim, obs, sim-obs for each season.
#
#########################################################################

import utilities.datareaders as readers
from utilities.regridder import regridder_global as regrid

import numpy as np
import xarray as xr
import xesmf as xe

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

import yaml
import glob
import calendar
from datetime import datetime

import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings('ignore',category=ShapelyDeprecationWarning)
warnings.filterwarnings('ignore',message='Input array is not C_CONTIGUOUS')


#########################################################################
#  MAIN FIGURE
#########################################################################

def plot_maps_comparisons(model,runid,season,y0,yf,varlist,simdata,obsdata,pdf_plots):

    # -- figure setup ---------------------------------------------------

    f,ax = plt.subplots(len(varlist),3,figsize=(12,2*len(varlist)),
                        subplot_kw={'projection':ccrs.PlateCarree()})
    title = '%s %s aerosol properties, %s %d-%d mean'%(model,runid,season,y0,yf)
    f.suptitle(title+'\n(all maps regridded to 2x2deg)')

    # -- subplot labels ---------------------------------------------------

    for i,var in enumerate(varlist):
        for j,src in enumerate([runid,'MISR','%s - MISR'%runid]):
            ax[i,j].set_title('%s %s'%(src,var.upper()))

    # -- colormap normalization -------------------------------------------

    cmaps = {'aod':'cividis', 'aaod':'afmhot_r', 'ae':'BuPu', 'dif':'coolwarm'}
    vmaxs = {'aod':0.5, 'aaod':0.1, 'ae':2}

    # -- plot maps -------------------------------------------------------

    for i,var in enumerate(varlist):

        sim = simdata[var]
        obs = obsdata[var]
        sim = sim.sel(time=sim.time.dt.season==season).mean('time')
        obs = obs.sel(time=obs.time.dt.season==season).mean('time')

        im0 = ax[i,0].pcolormesh(sim.lon,sim.lat,sim[var],cmap=cmaps[var],
                                 vmin=0,vmax=vmaxs[var],shading='auto',rasterized=True)
        im1 = ax[i,1].pcolormesh(obs.lon,obs.lat,obs[var],cmap=cmaps[var],
                                 vmin=0,vmax=vmaxs[var],shading='auto',rasterized=True)
        im2 = ax[i,2].pcolormesh(obs.lon,obs.lat,sim[var]-obs[var],cmap=cmaps['dif'],
                                 vmin=-vmaxs[var]/2,vmax=vmaxs[var]/2,shading='auto',rasterized=True)

        for j,im in enumerate([im0,im1,im2]):
            ax[i,j].coastlines()
            f.colorbar(im,ax=ax[i,j])

    plt.savefig(pdf_plots,format='pdf')

    return



#########################################################################
#  MAIN CALLS
#########################################################################

def call_seasonal_maps(model,runid,varlist,obsdicts,y0,yf):

    print('\nseasonal maps for %s %s'%(model,runid))
    pdf_plots = PdfPages('plots/seasonalmaps_%s.pdf'%runid)

    print('...reading+regridding model data')
    simdata = {var: regrid(readers.read_model_singlerlzn(model,runid,var,y0,yf),res=2)
               for var in varlist}

    print('...regridding obs data')
    obsdata = {var: regrid(obsdicts[var]['MISR'],res=2) for var in varlist}
        
    print('...plotting')
    for season in ['DJF','MAM','JJA','SON']:
        plot_maps_comparisons(model,runid,season,y0,yf,varlist,simdata,obsdata,pdf_plots)

    print('...cleaning up')
    for ds in simdata+obsdata: ds.close()
    pdf_plots.close()

    return