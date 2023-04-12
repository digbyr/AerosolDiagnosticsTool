
#################################################################
#
#  Regional Analyses
#
#  For now: use land mask to make land-vs-ocean plots.
#  Eventually: other region subsetting (eg lat-lon box).
#  May incorporate these options into the existing analyses
#  (comparemodels_globalproperties and singlemodel_seasonalmaps) 
#  or may keep as a separate script.
#
#################################################################

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
#  LAND vs OCEAN MAPS
#########################################################################

def land_ocn_maps(var,model,tag,simtype,y0,yf,pdf_plots):

    f,ax = plt.subplots(3,1,figsize=(8,8),subplot_kw={'projection':ccrs.PlateCarree()})
    ax[0].set_title('Global %s (70S-70N)'%var.upper())
    ax[1].set_title('Land Only')
    ax[2].set_title('Ocean Only')

    if simtype=='ens': 
        ds = readers.read_model_ensemble(model,tag,var,y0,yf).mean(('run','time'))
        if 'run' in ds.dims: ds = ds.chunk(dict(run=-1))
        lndmsk = readers.read_model_ensemble(model,tag,'lndmsk',None,None,aux=True)
        
    elif simtype=='ind':
        ds = readers.read_model_singlerlzn(model,tag,var,y0,yf).mean('time')
        lndmsk = readers.read_model_singlerlzn(model,tag,'lndmsk',None,None,aux=True)
            
    ds_all = ds[var]
    ds_lnd = ds[var].where(lndmsk.lndmsk.values)
    ds_ocn = ds[var].where(lndmsk.lndmsk.values==0)

    vmax = {'aod':0.5, 'aaod':0.03, 'dod':0.5, 'ae':2}[var]

    for a,dsi in zip(ax.flat,[ds_all,ds_lnd,ds_ocn]):
        a.pcolormesh(dsi.lon,dsi.lat,dsi,vmin=0,vmax=vmax,
                     shading='auto',rasterized=True)
        a.coastlines()

    plt.savefig(pdf_plots,format='pdf')
    
    return



    
#########################################################################
#  LAND vs OCEAN TIMESERIES
#########################################################################

def land_ocn_ts(var,ens_models,ens_expts,ind_models,ind_runids,
                obsdata,y0,yf,cdict,lsdict,pdf_plots):

    f,ax = plt.subplots(3,1,figsize=(8,8),sharey=True)
    ax[0].set_title('Global %s (70S-70N)'%var.upper())
    ax[1].set_title('Land Only')
    ax[2].set_title('Ocean Only')

    # -- plot sim data ------------------------------------------------------------
    
    simtypes = ['ens']*len(ens_models) + ['ind']*len(ind_models)
    for model,tag,simtype in zip(ens_models+ind_models,
                                 ens_expts+ind_runids,
                                 simtypes):
        
        if simtype=='ens': 
            ds = readers.read_model_ensemble(model,tag,var,y0,yf)
            if 'run' in ds.dims: ds = ds.chunk(dict(run=-1))
            lndmsk = readers.read_model_ensemble(model,tag,'lndmsk',None,None,aux=True)
            c,ls = cdict[model], lsdict[model]
            label = '%s\n(%s)'%(model,ds.dims['run'])

        elif simtype=='ind':
            ds = readers.read_model_singlerlzn(model,tag,var,y0,yf)
            lndmsk = readers.read_model_singlerlzn(model,tag,'lndmsk',None,None,aux=True)
            c,ls = cdict[tag], lsdict[tag]
            label = '%s\n(%s)'%(tag,model)
            
        if 'bnds' in lndmsk.dims: lndmsk = lndmsk.drop_dims('bnds')
        ds = ds.transpose('time','lon','lat')
        lndmsk = lndmsk.transpose('lon','lat')

        weights = np.cos(np.deg2rad(ds.lat))
        ts_all = ds[var].weighted(weights).mean(('lat','lon'))
        ts_lnd = ds[var].where(lndmsk.lndmsk.values).weighted(weights).mean(('lat','lon'))
        ts_ocn = ds[var].where(lndmsk.lndmsk.values==0).weighted(weights).mean(('lat','lon'))
        time = [np.datetime64(t) for t in ts_all.time.values]
        
        for i,ts in enumerate([ts_all,ts_lnd,ts_ocn]):
            if 'run' in ds.dims:
                ax[i].plot(time,ts.mean('run'),c=c,ls=ls,lw=3,label=label)
                ax[i].fill_between(time,ts.quantile(0.05,'run'),
                                   ts.quantile(0.95,'run'),
                                   fc=c,alpha=0.5)
            else:
                ax[i].plot(time,ts,c=c,ls=ls,lw=3,label=label)    
            ts.close()

        for dsi in [ds,lndmsk]: dsi.close()


    # -- plot obs data ------------------------------------------------------------
    
    # need to use the sim lndmsk for the obs -- regrid both to 1x1
    lndmsk = regrid(lndmsk,res=1).transpose('lon','lat')

    for sat in obsdata.keys():
        
        if (sat=='POLDER-GRASP') & (y0>2012): continue

        ds = regrid(obsdata[sat],res=1).transpose('time','lon','lat')
        weights = np.cos(np.deg2rad(ds.lat))
        ts_all = ds[var].weighted(weights).mean(('lat','lon'))
        ts_lnd = ds[var].where(lndmsk.lndmsk.values).weighted(weights).mean(('lat','lon'))
        ts_ocn = ds[var].where(lndmsk.lndmsk.values==0).weighted(weights).mean(('lat','lon'))
        time = [np.datetime64(t) for t in ts_all.time.values]
            
        for i,ts in enumerate([ts_all,ts_lnd,ts_ocn]):
            ax[i].plot(time,ts,c=cdict[sat],ls=lsdict[sat],lw=3,label=sat)
            ts.close()
    
    # -- finish fig ------------------------------------------------------------
    
    ax[0].set_ylabel(var)
    ax[0].legend()
    plt.tight_layout()
    
    plt.savefig(pdf_plots,format='pdf')
    
    return

    

#########################################################################
#  LAND vs OCEAN TAYLOR DIAGRAMS
#########################################################################

def land_ocn_taylor(var,ens_models,ens_expts,ind_models,ind_runids,
                    obsdata,y0,yf,cdict,mkdict,pdf_plots):

    f = plt.figure(figsize=(12,4))
    f.suptitle('%d-%d mean %s, all maps regridded to 2x2deg'%(y0,yf,var.upper()))

    # read in one reference landmask, regrid to 2x2
    if len(ind_models)>0: 
        lndmsk = readers.read_model_singlerlzn(ind_models[0],ind_runids[0],
                                               'lndmsk',None,None,aux=True)
    else: 
        lndmsk = readers.read_model_ensemble(ens_models[0],ens_expts[0],
                                             'lndmsk',None,None,aux=True)
    lndmsk = regrid(lndmsk,res=2).transpose('lon','lat')

    # set up reference (obs) data: time-mean maps, regridded by 2x2
    satfields = [regrid(ds.mean('time'),res=2).transpose('lon','lat') 
                 for sat,ds in obsdata.items()]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore','Mean of empty slice')
        satmean = xr.ones_like(satfields[0][var])
        satmean.values = np.nanmean([ds[var].values for ds in satfields],axis=0)
    
    # 3 reference fields: all, land, ocean
    satrefs = [satmean,
               satmean.where(lndmsk.lndmsk.values),
               satmean.where(lndmsk.lndmsk.values==0)]

    # build taylor diagram axes
    ax = [TaylorDiagram(np.nanstd(satref),fig=f,rect=rect,
                        label='obs mean',srange=(0,3)) 
          for satref,rect in zip(satrefs,[131,132,133])]
    
    ax[0]._ax.set_title('Global %s (70S-70N)'%var.upper())
    ax[1]._ax.set_title('Land Only')
    ax[2]._ax.set_title('Ocean Only')

    # -- iterate over sims ------------------------------------------------------

    simtypes = ['ens']*len(ens_models) + ['ind']*len(ind_models)
    for model,tag,simtype in zip(ens_models+ind_models,
                                 ens_expts+ind_runids,
                                 simtypes):
        
        if simtype=='ens': 
            ds = readers.read_model_ensemble(model,tag,var,y0,yf)
            if 'run' in ds.dims: ds = ds.chunk(dict(run=-1))
            c,mk = cdict[model], mkdict[model]
            label = '%s\n(%s)'%(model,ds.dims['run'])

        elif simtype=='ind':
            ds = readers.read_model_singlerlzn(model,tag,var,y0,yf)
            c,mk = cdict[tag], mkdict[tag]
            label = '%s\n(%s)'%(tag,model)

        ds = regrid(ds.mean('time'),res=2).transpose('lon','lat')
        ds_lnd = ds.where(lndmsk.lndmsk.values)
        ds_ocn = ds.where(lndmsk.lndmsk.values==0)

        # plot points
        for a,dsi,satref in zip(ax,[ds,ds_lnd,ds_ocn],satrefs):

            if 'run' in dsi.dims:
                for run in dsi.run.values: 
                    plot_taylor_point(dsi[var].sel(run=run).values,
                                      satref.values,a,c,mk,label,'minor')
                plot_taylor_point(dsi[var].mean('run').values,
                                  satref,a,c,mk,label,'major')
            else: 
                plot_taylor_point(dsi[var].values,satref.values,
                                  a,c,mk,label,'major')    
            dsi.close()


    # -- iterate over obs ------------------------------------------------------

    for sat,ds in obsdata.items():
        ds = regrid(ds.mean('time'),res=2).transpose('lon','lat')
        ds_lnd = ds.where(lndmsk.lndmsk.values)
        ds_ocn = ds.where(lndmsk.lndmsk.values==0)
        c,mk,label = cdict[sat], mkdict[sat], sat
        for a,dsi,satref in zip(ax,[ds,ds_lnd,ds_ocn],satrefs):
            plot_taylor_point(dsi[var].values,satref.values,
                              a,c,mk,label,'major')  


    # -- finish figure ------------------------------------------------------
        
    for a in ax:
        a.add_grid()
        a.add_contours(levels=5,colors='0.9')
        
    f.legend(ax[0].samplePoints,[p.get_label() for p in ax[0].samplePoints],
             numpoints=1,loc='upper left',bbox_to_anchor=(1.1,1),
             bbox_transform=ax[-1]._ax.transAxes)

    plt.subplots_adjust(wspace=0.4,right=0.8,bottom=0.02,top=0.99)
    plt.savefig(pdf_plots,format='pdf')

    return





#########################################################################
#  MAIN CALLS
#########################################################################

def call_compare_landvsocn(figtag,ens_models,ens_expts,ind_models,ind_runids,
                           varlist,obsdicts,y0,yf):

    print('\land-vs-ocn comparisons for %s and %s'\
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

    print('...plotting')
    for var in varlist:
        
        print('......%s land-ocean maps'%var)
        land_ocn_maps(var,ind_models[-1],ind_runids[-1],'ind',y0,yf,pdf_plots)
        
        print('......%s land-ocean timeseries'%var)
        land_ocn_ts(var,ens_models,ens_expts,
                    ind_models,ind_runids,
                    obsdicts[var],y0,yf,
                    cdict,lsdict,pdf_plots)
        
        if var!='ae':
            print('......%s land-ocean taylors'%var)
            land_ocn_taylor(var,ens_models,ens_expts,
                            ind_models,ind_runids,
                            obsdicts[var],y0,yf,
                            cdict,mkdict,pdf_plots)
        
    
    pdf_plots.close()

    return

