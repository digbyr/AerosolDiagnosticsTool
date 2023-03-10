#########################################################################
#
#  COMPARE MAPS: GLOBAL PROPERTIES
#
#  For a selection of models (could be ensembles, individual realizations,
#  or a mix thereof) plot globally-averaged properties.
#  Includes time series, seasonal cycle, zonal means, taylor diagram.
#  "Global" actually means (70S, 70N) to allow comparison with obs.
#
#########################################################################

import utilities.datareaders as readers
from utilities.regridder import regridder_global as regrid
from utilities.taylorDiagram import TaylorDiagram

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
#  PLOTS: TIMESERIES
#########################################################################

def plot_timeseries_seasons_combined(var,ens_models,ens_expts,
                                     ind_models,ind_runids,
                                     obsdata,y0,yf,
                                     cdict,lsdict,pdf_plots):


    f = plt.figure(figsize=(16,4))
    gs = GridSpec(1,2,width_ratios=[2,1])
    ax = np.array([f.add_subplot(gs[0]),f.add_subplot(gs[1])])
    
    ax[0].set_title('%s, %d-%d (70S-70N)'%(var.upper(),y0,yf))
    ax[1].set_title('mean seasonal cycle')
    
    months = np.arange(1,13,1)
    
    # -- plot sim data ------------------------------------------------------------
    
    simtypes = ['ens']*len(ens_models) + ['ind']*len(ind_models)

    for model,tag,simtype in zip(ens_models+ind_models,
                                 ens_expts+ind_runids,
                                 simtypes):
        
        if simtype=='ens': 
            ds = readers.read_model_ensemble(model,tag,var,y0,yf)
            if 'run' in ds.dims: ds = ds.chunk(dict(run=-1))
            c,ls = cdict[model], lsdict[model]
            label = '%s\n(%s)'%(model,ds.dims['run'])

        elif simtype=='ind':
            ds = readers.read_model_singlerlzn(model,tag,var,y0,yf)
            c,ls = cdict[tag], lsdict[tag]
            label = '%s\n(%s)'%(tag,model)
            
        ts = ds[var].weighted(np.cos(np.deg2rad(ds.lat))).mean(('lat','lon'))   
        ts_season = ts.groupby(ts.time.dt.month).mean('time')
        time = [np.datetime64(t) for t in ts.time.values]
        
        if 'run' in ds.dims:

            ax[0].plot(time,ts.mean('run'),c=c,ls=ls,lw=3,label=label)
            ax[1].plot(months,ts_season.median('run'),c=c,ls=ls,lw=3,label=label)

            ax[0].fill_between(time,ts.quantile(0.05,'run'),
                               ts.quantile(0.95,'run'),
                               fc=c,alpha=0.5)
            ax[1].fill_between(months,ts_season.quantile(0.05,'run'),
                               ts_season.quantile(0.95,'run'),
                               fc=c,alpha=0.5)            
        else:
            
            ax[0].plot(time,ts.values,c=c,ls=ls,lw=3,label=label)
            ax[1].plot(months,ts_season.values,c=c,ls=ls,lw=3,label=label)

        for dsi in [ds,ts,ts_season]: dsi.close()
        
        
    # -- plot obs data ------------------------------------------------------------
    
    for sat in obsdata.keys():
        
        ds = obsdata[sat]
        ts = ds[var].weighted(np.cos(np.deg2rad(ds.lat))).mean(('lat','lon'))   
        ts_season = ts.groupby(ts.time.dt.month).mean('time')
        time = [np.datetime64(t) for t in ts.time.values]
            
        ax[0].plot(time,ts,c=cdict[sat],ls=lsdict[sat],lw=3,label=sat)
        ax[1].plot(months,ts_season,c=cdict[sat],ls=lsdict[sat],lw=3,label=sat)    
    
    
    # -- finish fig ------------------------------------------------------------
    
    ax[0].set_ylabel(var)
    ax[1].set_xticks(months[::3])
    ax[1].set_xticklabels([calendar.month_abbr[m] for m in months[::3]])
        
    # force shared y axis
    ylim = ax[0].get_ylim()
    for i in range(2): 
        ax[i].set_ylim(ylim)
    ax[1].set_yticklabels('')
    
    plt.subplots_adjust(left=0.05,right=0.85,wspace=0.05)
    ax[0].legend(loc='upper left',bbox_to_anchor=(1.05,0.95),
                 bbox_transform=ax[1].transAxes)
    
    plt.savefig(pdf_plots,format='pdf')
    
    return


####################################################################
#  PLOTS: ZONAL MEANS
####################################################################

def plot_zonal_mean(var,ens_models,ens_expts,
                    ind_models,ind_runids,
                    obsdata,y0,yf,
                    cdict,lsdict,pdf_plots):
    
    f,ax = plt.subplots(1,figsize=(8,8))
    ax.set_title('Zonal Mean %s (%d-%d mean)'%(var.upper(),y0,yf))
    
    # -- plot sim data ------------------------------------------------------------
    
    simtypes = ['ens']*len(ens_models) + ['ind']*len(ind_models)

    for model,tag,simtype in zip(ens_models+ind_models,
                                 ens_expts+ind_runids,
                                 simtypes):
        
        if simtype=='ens': 
            ds = readers.read_model_ensemble(model,tag,var,y0,yf)
            if 'run' in ds.dims: ds = ds.chunk(dict(run=-1))
            c,ls = cdict[model], lsdict[model]
            label = '%s (%s)'%(model,ds.dims['run'])

        elif simtype=='ind':
            ds = readers.read_model_singlerlzn(model,tag,var,y0,yf)
            c,ls = cdict[tag], lsdict[tag]
            label = '%s (%s)'%(tag,model)
                
        zm = ds[var].mean(('lon','time'))
            
        if 'run' in ds.dims: 
            ax.plot(zm.lat,zm.mean('run'),c=c,ls=ls,lw=3,label=label)
            ax.fill_between(zm.lat,zm.quantile(0.05,'run'),
                            zm.quantile(0.95,'run'),
                            fc=c,alpha=0.5)
        else: 
            ax.plot(zm.lat,zm,c=c,ls=ls,lw=3,label=label)

        for dsi in [ds,zm]: dsi.close()
    

    # -- plot obs data ------------------------------------------------------------
    
    for sat in obsdata.keys():
        
        ds = obsdata[sat]
        zm = ds[var].mean(('lon','time'))
        ax.plot(zm.lat,zm,c=cdict[sat],ls=lsdict[sat],lw=3,label=sat)
    
    
    # -- finish fig ------------------------------------------------------------
    
    ax.legend(loc='upper left')
    ax.set_xlabel('latitude')
    ax.set_ylabel(var.upper())
    
    plt.savefig(pdf_plots,format='pdf')

    return




####################################################################
#  PLOTS: TAYLOR DIAGRAM
####################################################################

def plot_taylor_diagram(varlist,ens_models,ens_expts,
                        ind_models,ind_runids,
                        obsdicts,y0,yf,
                        cdict,mkdict,pdf_plots):

    # -- figure setup ------------------------------------------------------
    
    f = plt.figure(figsize=(8,6))
    f.suptitle('%d-%d mean, all maps regridded to 2x2deg'%(y0,yf))
    
    # set up reference (obs) data: time-mean maps, regridded by 2x2
    satrefs = {}
    for var in varlist:
        satfields = [regrid(ds.mean('time'),res=2) for sat,ds in obsdicts[var].items()]
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore','Mean of empty slice')
            satrefs[var] = np.nanmean([ds[var].values for ds in satfields],axis=0)
    
    # build taylor diagram axes
    ax = []
    rects = {1:[111], 2:[121,122], 3:[221,222,223], 4:[221,222,223,224]}
    for var,rect in zip(varlist,rects[len(varlist)]):
        ax.append(TaylorDiagram(np.nanstd(satrefs[var]),fig=f,rect=rect,
                                label='obs mean',srange=(0,3)))

    
    # -- iterate over sims ------------------------------------------------------

    simtypes = ['ens']*len(ens_models) + ['ind']*len(ind_models)
    for model,tag,simtype in zip(ens_models+ind_models,
                                 ens_expts+ind_runids,
                                 simtypes):
        
        for i,var in enumerate(varlist):

            # read in data: ensembles
            if simtype=='ens': 
                ds = readers.read_model_ensemble(model,tag,var,y0,yf)
                if 'run' in ds.dims: ds = ds.chunk(dict(run=-1))
                c,mk = cdict[model], mkdict[model]
                label = '%s\n(%s)'%(model,ds.dims['run'])

            # read in data: single realizations
            elif simtype=='ind':
                ds = readers.read_model_singlerlzn(model,tag,var,y0,yf)
                c,mk = cdict[tag], mkdict[tag]
                label = '%s\n(%s)'%(tag,model)
                    
            # process data
            ds = regrid(ds.mean('time'),res=2)
            
            # plot point(s)
            if 'run' in ds.dims:
                for run in ds.run.values: 
                    plot_taylor_point(ds[var].sel(run=run).values,satrefs[var],
                                      ax[i],c,mk,label,'minor')
                plot_taylor_point(ds[var].mean('run').values,satrefs[var],
                                  ax[i],c,mk,label,'major')
            else: 
                plot_taylor_point(ds[var].values,satrefs[var],
                                  ax[i],c,mk,label,'major')    

            ds.close()
        

    # -- iterate over obs ------------------------------------------------------

    for i,var in enumerate(varlist):
        for sat,ds in obsdicts[var].items():
            ds = regrid(ds.mean('time'),res=2)
            c,mk,label = cdict[sat], mkdict[sat], sat
            plot_taylor_point(ds[var].values,satrefs[var],
                              ax[i],c,mk,label,'major')    


    # -- finish figure ------------------------------------------------------
        
    for i,var in enumerate(varlist):
        ax[i]._ax.set_title(var,fontweight='bold')
        ax[i].add_grid()
        ax[i].add_contours(levels=5,colors='0.9')
        ax[i]._ax.tick_params(axis='x',rotation=45)

    # one legend for whole fig  
    pts,labs = ax[0].samplePoints, [p.get_label() for p in ax[0].samplePoints]
    for a in ax[1:]: 
        ptsi = a.samplePoints
        pts.extend([pt for pt in ptsi if pt.get_label() not in labs])
        labs.extend([pt.get_label() for pt in ptsi if pt.get_label() not in labs])
    anchor = ax[1] if len(varlist)==4 else ax[-1]
    f.legend(pts,labs,numpoints=1,loc='upper left',bbox_to_anchor=(1.1,1),
             bbox_transform=anchor._ax.transAxes)
    
    right = 0.9 if len(varlist)==3 else 0.7
    plt.subplots_adjust(hspace=0.3,wspace=0.3,right=right)
    plt.savefig(pdf_plots,format='pdf')

    return


def plot_taylor_point(dat,ref,ax,c,mk,label,ptype):

    # calculate Pearson correlation coefficient (ie how tight is scatter)
    msk = np.isfinite(dat)*np.isfinite(ref)
    std = np.nanstd(dat)
    cor = stats.pearsonr(ref[msk],dat[msk])[0]

    # ptype (point type) = major or minor
    # major is used for ensemble medians and single-realization products (eg obs)
    # minor is used for individual realizations in an ensemble
    
    if ptype=='major':
        ax.add_sample(std,cor,c=c,marker=mk,mec='k',ms=10,ls='',label=label)
    elif ptype=='minor':
        ax.add_sample(std,cor,c=c,marker='.',alpha=0.5,ms=5,ls='',label='')

    return




#########################################################################
#  MAIN CALLS
#########################################################################

def call_compare_global(figtag,ens_models,ens_expts,ind_models,ind_runids,
                        varlist,obsdicts,y0,yf):

    print('\nglobal comparisons for %s and %s'\
          %(', '.join(ens_models),', '.join(ind_runids)))
    pdf_plots = PdfPages('plots/globalproperties_%s.pdf'%figtag)

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

        print('......%s timeseries'%var)
        plot_timeseries_seasons_combined(var,ens_models,ens_expts,
                                         ind_models,ind_runids,
                                         obsdicts[var],y0,yf,
                                         cdict,lsdict,pdf_plots)
        
        print('......%s zonal mean'%var)
        plot_zonal_mean(var,ens_models,ens_expts,
                        ind_models,ind_runids,
                        obsdicts[var],y0,yf,
                        cdict,lsdict,pdf_plots)
    
    print('......multi-var taylor diagram')
    plot_taylor_diagram(varlist,ens_models,ens_expts,
                        ind_models,ind_runids,
                        obsdicts,y0,yf,
                        cdict,mkdict,pdf_plots)

    pdf_plots.close()

    return

