#########################################################################
#
#  Land vs Ocean Figs
#
#  One massive function that makes all the plots simultaneously.
#  Maybe not as elegant as the other scripts, but it gets the job done.
#
#  There is code here, currently commented out, for plotting maps 
#  of the total, land, and ocean fields in a subset of model+obs
#  datasets. This was really just for debugging purposes so is 
#  not included in the default output, but can be uncommented if desired.
#
#########################################################################

import utilities.datareaders as readers
from utilities.regridder import regridder_global as regrid
from utilities.taylorDiagram import TaylorDiagram
from analysis_scripts.comparemodels_globalproperties import plot_taylor_point

import numpy as np
import xarray as xr
import xesmf as xe
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

import yaml
import glob
import calendar
from tqdm import tqdm
from datetime import datetime

import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings('ignore',category=ShapelyDeprecationWarning)
warnings.filterwarnings('ignore',message='Input array is not C_CONTIGUOUS')
warnings.filterwarnings('ignore',message='The handle <matplotlib.lines.Line2D object')



#########################################################################
#  ONE FUNCTION TO BRUTE FORCE THEM ALL
#########################################################################

def land_ocn_figs(var,ens_models,ens_expts,ind_models,ind_runids,
                  obsdata,y0,yf,cdict,mkdict,lsdict,pdf_plots):
    
    # -- set up all the figures --------------------------------------------------
    
    print('......setting up')

    #fmap,axmap = plt.subplots(3,4,figsize=(12,4),subplot_kw={'projection':ccrs.PlateCarree()})
    fts,axts = plt.subplots(3,1,figsize=(8,8),sharey=True)
    ftay,axtay = plt.figure(figsize=(12,4)), []

    #fmap.suptitle('%s sample masking (2x2deg, limited to 70S-70N) '%var.upper())
    fts.suptitle('%s mean, 2x2deg, 70S-70N'%var.upper())
    for i,loc in enumerate(['All','Land Only','Ocean Only']): 
        axts[i].set_title(loc,fontweight='bold')
        axts[i].set_ylabel(var.upper())
    ftay.suptitle('%s, 70S-70N, regridded to 2x2deg'%var.upper())

    # -- landmasks to use for all obs (1x1) and for all data in taylor diagram (2x2) ----------

    if len(ind_models)>0: 
        lndmsk = readers.read_model_singlerlzn(ind_models[0],ind_runids[0],'lndmsk',None,None,aux=True)
    else: 
        lndmsk = readers.read_model_ensemble(ens_models[0],ens_expts[0],'lndmsk',None,None,aux=True)
    if 'bnds' in lndmsk.dims: lndmsk = lndmsk.drop_dims('bnds')

    lndmsk2x2 = regrid(lndmsk,res=2).sel(lat=slice(-70,70))


    # -- set up axes for taylor diagram -------------------------------------------
    
    # obs reference: mean of time-mean fields
    satfields = [regrid(ds.mean('time'),res=2).sel(lat=slice(-70,70))
                 for sat,ds in obsdata.items()]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore','Mean of empty slice')
        satmean = xr.ones_like(satfields[0][var])
        satmean.values = np.nanmean([ds[var].values for ds in satfields],axis=0)
    
    # 3 reference fields: all, land, ocean
    satrefs = [satmean,
               satmean.where(lndmsk2x2.lndmsk.values),
               satmean.where(lndmsk2x2.lndmsk.values==0)]

    # build taylor diagram axes
    axtay = [TaylorDiagram(np.nanstd(satref),fig=ftay,rect=rect,
                           label='obs mean',srange=(0,2.5)) 
             for satref,rect in zip(satrefs,[131,132,133])]
    for i,loc in enumerate(['All','Land Only','Ocean Only']): 
        axtay[i]._ax.set_title(loc,fontweight='bold')
    
    # -- iterate through sim data ------------------------------------------------------------

    print('......plotting sims')

    simtypes = ['ens']*len(ens_models) + ['ind']*len(ind_models)
    for m,(model,tag,simtype) in enumerate(zip(ens_models+ind_models,
                                               ens_expts+ind_runids,
                                               simtypes)):
        if simtype=='ens': 
            ds = readers.read_model_ensemble(model,tag,var,y0,yf)
            if 'run' in ds.dims: ds = ds.chunk(dict(run=-1))
            c,ls,mk = cdict[model], lsdict[model], mkdict[model]
            label = '%s (%s)'%(model,ds.dims['run'])

        elif simtype=='ind':
            ds = readers.read_model_singlerlzn(model,tag,var,y0,yf)
            c,ls,mk = cdict[tag], lsdict[tag], mkdict[tag]
            label = '%s (%s)'%(tag,model)

        # (lat,lon,time) -> (time,lat,lon) to allow masking with (lat,lon)
        ds = regrid(ds,res=2).transpose('time',...).sel(lat=slice(-70,70))
        
        # do masking
        ds_all = ds[var]
        ds_lnd = ds[var].where(lndmsk2x2.lndmsk.values)
        ds_ocn = ds[var].where(lndmsk2x2.lndmsk.values==0)

        '''
        # first two models: plot maps to check masking
        vmax = {'aod':0.5, 'aaod':0.03, 'dod':0.5, 'ae':1.5}[var]
        if m<2:
            axmap[0,m].set_title(label.split('(')[0])
            for i,dsi in enumerate([ds_all,ds_lnd,ds_ocn]):
                im = axmap[i,m].pcolormesh(dsi.lon,dsi.lat,dsi.mean('time'),vmin=0,vmax=vmax,
                                           shading='auto',rasterized=True)
                axmap[i,m].coastlines()
        '''
        
        # plot timeseries
        weights = np.cos(np.deg2rad(ds.lat))
        time = [np.datetime64(t) for t in ds.time.values]
        for i,dsi in enumerate([ds_all,ds_lnd,ds_ocn]):
            tsi = dsi.weighted(weights).mean(('lat','lon'))
            if 'run' in ds.dims:
                axts[i].plot(time,tsi.mean('run'),c=c,ls=ls,lw=3,label=label)
                axts[i].fill_between(time,tsi.quantile(0.05,'run'),
                                     tsi.quantile(0.95,'run'),
                                     fc=c,alpha=0.5)
            else:
                axts[i].plot(time,tsi,c=c,ls=ls,lw=3,label=label)    
            tsi.close()

        # plot taylor fig
        for a,dsi,satref in zip(axtay,[ds_all,ds_lnd,ds_ocn],satrefs):
            dsi = dsi.mean('time')
            if 'run' in dsi.dims:
                for run in dsi.run.values: 
                    plot_taylor_point(dsi.sel(run=run).values,
                                      satref.values,a,c,mk,label,'minor')
                plot_taylor_point(dsi.mean('run').values,
                                  satref,a,c,mk,label,'major')
            else: 
                plot_taylor_point(dsi.values,satref.values,
                                  a,c,mk,label,'major')    
            dsi.close()
        
        
    # -- iterate over obs ------------------------------------------------------

    print('......plotting obs')

    for s,sat in enumerate(obsdata.keys()):
        
        # regrid to standard 1x1 to allow masking with sim lndmsk
        ds = regrid(obsdata[sat],res=2).sel(lat=slice(-70,70)).transpose('time',...)
        c,ls,mk,label = cdict[sat],lsdict[sat],mkdict[sat],sat
        
        # do masking
        ds_all = ds[var]
        ds_lnd = ds[var].where(lndmsk2x2.lndmsk.values)
        ds_ocn = ds[var].where(lndmsk2x2.lndmsk.values==0)

        '''
        # first two sats: plot maps to check masking
        vmax = {'aod':0.5, 'aaod':0.03, 'dod':0.5, 'ae':2}[var]
        if s<2:
            axmap[0,s+2].set_title(label)
            for i,dsi in enumerate([ds_all,ds_lnd,ds_ocn]):
                axmap[i,s+2].pcolormesh(dsi.lon,dsi.lat,dsi.mean('time'),vmin=0,vmax=vmax,
                                        shading='auto',rasterized=True)
                axmap[i,s+2].coastlines()
        '''
        
        # plot timeseries
        if not ((sat=='POLDER-GRASP')*(y0>2012)):
            weights = np.cos(np.deg2rad(ds.lat))
            time = [np.datetime64(t) for t in ds.time.values]
            for i,dsi in enumerate([ds_all,ds_lnd,ds_ocn]):
                tsi = dsi.weighted(weights).mean(('lat','lon'))
                axts[i].plot(time,tsi,c=c,ls=ls,lw=3,label=label)    
                tsi.close()

        # plot taylor fig
        for a,dsi,satref in zip(axtay,[ds_all,ds_lnd,ds_ocn],satrefs):
            dsi = dsi.mean('time')
            plot_taylor_point(dsi.values,satref.values,a,c,mk,label,'major')    
            dsi.close()

        
    # -- finish figs ------------------------------------------------------------
    
    print('......finishing')
    
    # maps
    #fmap.colorbar(im,ax=axmap.ravel().tolist(),label=var.upper())

    # timeseries
    axts[0].legend(fontsize=8)
    fts.tight_layout()

    # taylor
    for a in axtay: 
        a.add_grid()
        a.add_contours(levels=5,colors='0.9')
    ftay.legend(axtay[0].samplePoints,[p.get_label() for p in axtay[0].samplePoints],
                numpoints=1,loc='upper left',bbox_to_anchor=(1.1,1),
                bbox_transform=axtay[-1]._ax.transAxes,fontsize=8)
    ftay.subplots_adjust(wspace=0.4,left=0.07,right=0.75,bottom=0.02,top=0.99)

    # save figs
    #for fig in [fmap,fts,ftay]: fig.savefig(pdf_plots,format='pdf')
    for fig in [fts,ftay]: fig.savefig(pdf_plots,format='pdf')

    return




#########################################################################
#  MAIN CALLS
#########################################################################

def call_compare_landvsocn(figtag,ens_models,ens_expts,ind_models,ind_runids,
                           varlist,obsdicts,y0,yf):

    print('land-vs-ocn comparisons for %s and %s'\
          %(', '.join(ens_models),', '.join(ind_runids)))
    pdf_plots = PdfPages('plots/landvsocn_%s.pdf'%figtag)

    simtags = ens_models+ind_runids
    sats = ['MODIS Aqua','MISR','CALIOP','MIDAS','ACROS','POLDER-GRASP']

    # set colors to be consistent from plot to plot
    cdict_sim = {tag:plt.cm.hsv(i/len(simtags)) for i,tag in enumerate(simtags)}
    cdict_obs = {sat:plt.cm.gist_earth(i/len(sats)) for i,sat in enumerate(sats)}
    cdict = {**cdict_sim,**cdict_obs}

    # set linestyles to be consistent from plot to plot
    lsdict_sim = {tag:'-' for tag in simtags}
    lsdict_obs = {'MODIS Aqua':'--', 'MISR':':', 'CALIOP':'-.', 
                  'MIDAS':'--', 'ACROS':'-.','POLDER-GRASP':(0,(4,1,1,1,1,1,1,1,1,1))}
    lsdict = {**lsdict_sim,**lsdict_obs}

    # set markerstyles to be consistent from plot to plot (currently only used in Taylor fig)
    mkdict_sim = {tag:'o' for tag in simtags}
    mkdict_obs = {sat:'s' for sat in sats}
    mkdict = {**mkdict_sim,**mkdict_obs}

    for var in varlist:
        
        print('...plotting %s'%var.upper())
        land_ocn_figs(var,ens_models,ens_expts,ind_models,ind_runids,
                      obsdicts[var],y0,yf,cdict,mkdict,lsdict,pdf_plots)    
    
    pdf_plots.close()

    return



