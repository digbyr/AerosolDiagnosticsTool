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

def plot_maps_comparisons(model,tag,simtype,season,y0,yf,varlist,simdata,obsdata,pdf_plots):

    # -- figure setup ---------------------------------------------------

    f,ax = plt.subplots(2*len(varlist),3,figsize=(12,4*len(varlist)),
                        subplot_kw={'projection':ccrs.PlateCarree()})
    ensnote = ' (1st realization)' if simtype=='ens' else ''
    title = '%s %s%s aerosol properties, %s %d-%d mean'%(model,tag,ensnote,season,y0,yf)
    f.suptitle(title+'\n(all maps regridded to 2x2deg)')

    satyrs = {'MISR':(max(y0,2001),min(2020,yf)),
              'POLDER-GRASP':(max(y0,2006),min(2012,yf)),
              'MODIS Aqua':(max(y0,2003),min(2020,yf)),
              'ACROS':(max(y0,2007),min(2020,yf)),
              'MIDAS':(max(y0,2003),min(2017,yf)),
              'CALIOP':(max(y0,2007),min(2020,yf))}

    # -- subplot labels ---------------------------------------------------

    simlabel = model if simtype=='ens' else tag
    for i,var in enumerate(varlist):

        # sim: plot every other row (1x/var)
        ax[2*i,0].set_title('%s %s'%(simlabel,var.upper()),fontsize=8)

        # obs: plot every row (2 sats/var)
        for ii,sat in enumerate(obsdata[var].keys()):
            if sat=='POLDER-GRASP': yrs = '(2006-2012 climatology)'
            elif satyrs[sat][0]>satyrs[sat][1]: yrs = '(no obs data available)'
            else: yrs = '(%d-%d)'%(satyrs[sat][0],satyrs[sat][1])
            for j,src in enumerate([sat,'%s - %s'%(simlabel,sat)]):
                ax[2*i+ii,j+1].set_title('%s %s\n%s'%(src,var.upper(),yrs),fontsize=8)

    # -- colormap normalization -------------------------------------------

    cmaps = {'aod':'cividis', 'aaod':'afmhot_r', 'dod':'hot_r',
             'ae':'BuPu', 'dif':'coolwarm'}
    vmaxs = {'aod':0.5, 'aaod':0.03, 'dod':0.5, 'ae':2}

    # -- plot maps -------------------------------------------------------

    for i,var in enumerate(varlist):

        sim = simdata[var]
        simfull = sim.sel(time=sim.time.dt.season==season).mean('time')
        
        im0 = ax[2*i,0].pcolormesh(sim.lon,sim.lat,simfull[var],cmap=cmaps[var],
                                   vmin=0,vmax=vmaxs[var],shading='auto',rasterized=True)
        f.colorbar(im0,ax=ax[2*i,0],shrink=0.9)
        ax[2*i,0].coastlines()
        ax[2*i+1,0].axis('off')
        
        for ii,sat in enumerate(obsdata[var]):
            if sat=='POLDER-GRASP':
                obstrim = obsdata[var][sat]
            else:
                yi0,yif = str(satyrs[sat][0]),str(satyrs[sat][1])
                obstrim = obsdata[var][sat].sel(time=slice(yi0,yif))
            simtrim = sim.sel(time=slice(yi0,yif))
            obstrim = obstrim.sel(time=obstrim.time.dt.season==season).mean('time')
            simtrim = simtrim.sel(time=simtrim.time.dt.season==season).mean('time')

            im1 = ax[2*i+ii,1].pcolormesh(obstrim.lon,obstrim.lat,obstrim[var],cmap=cmaps[var],
                                          vmin=0,vmax=vmaxs[var],shading='auto',rasterized=True)
            im2 = ax[2*i+ii,2].pcolormesh(obstrim.lon,obstrim.lat,simtrim[var]-obstrim[var],cmap=cmaps['dif'],
                                          vmin=-vmaxs[var]/2,vmax=vmaxs[var]/2,shading='auto',rasterized=True)

            for j,im in enumerate([im1,im2]): 
                f.colorbar(im,ax=ax[2*i+ii,j+1],shrink=0.9)
                ax[2*i+ii,j+1].coastlines()

    plt.subplots_adjust(top=0.92,bottom=0.02)
    plt.savefig(pdf_plots,format='pdf')

    return



#########################################################################
#  MAIN CALLS
#########################################################################

def call_seasonal_maps_2obs(model,tag,simtype,varlist,obsdicts,y0,yf):

    print('\nseasonal maps for %s %s'%(model,tag))
    simlabel = model if simtype=='ens' else tag
    pdf_plots = PdfPages('plots/seasonalmaps_%s.pdf'%simlabel)

    print('...reading+regridding model data')
    if simtype=='ens':
        tmpdata = [readers.read_model_ensemble(model,tag,var,y0,yf) for var in varlist]
        run = tmpdata[0].run.values[0]
        simdata = {var: regrid(tmp.sel(run=run).squeeze(),res=2) for var,tmp in zip(varlist,tmpdata)}
        for tmp in tmpdata: tmp.close()
    elif simtype=='ind':
        simdata = {var: regrid(readers.read_model_singlerlzn(model,tag,var,y0,yf),res=2)
                   for var in varlist}

    print('...regridding obs data')
    mapobs = {'aod':['MISR','POLDER-GRASP'],
              'aaod':['MISR','POLDER-GRASP'],
              'dod':['ACROS','MIDAS'],
              'ae':['MISR','POLDER-GRASP']}
    obsdata = {var: {sat: regrid(obsdicts[var][sat],res=2) for sat in mapobs[var]}
               for var in varlist }

    print('...plotting')
    for season in ['DJF','MAM','JJA','SON']:
        plot_maps_comparisons(model,tag,simtype,season,y0,yf,varlist,simdata,obsdata,pdf_plots)

    print('...cleaning up')
    for var,ds in simdata.items(): ds.close()
    for var in obsdata.keys():
        for sat,ds in obsdata[var].items(): ds.close()
    pdf_plots.close()

    return
