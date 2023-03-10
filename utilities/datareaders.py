########################################################################### 
#
#  DATA READER
#  Co-developed by Ruth Digby and Yi Ren
#
#  A script to read and standardize data from various observation and 
#  simulation datasets. Uses config files to set data paths. 
#
#  Functions accept t0,tf (start and end times) as int or str. 
#  This assumes that an int will be the year only (eg 2020)
#  And that str will be either 'YYYY' or 'YYYY-MM'. 
#  Other formats are not accepted (except I guess 'YYYY-MM-DD').
#
#  If there is no data available for some or all of the requested
#  period, we return "dummy" data (same resoln but all nan's).
#  Variable convention: 
#    - ymin,ymax = years that have data available
#    - years = years for which data have been requested
#
#  All output is in (lat,lon,time) with time as monthly means.
#
#  Fields get reordered onto lon = (-180,180), lat = (-90,90).
#  There are strong arguments to be made for (0,360) but personally I 
#  dislike plotting it that way, so (-180,180) it is. 
#
###########################################################################


import numpy as np
import xarray as xr
import pandas as pd
from pyhdf.SD import SD,SDC
from netCDF4 import Dataset as ncDataset

import os
import glob
import yaml
import calendar
from datetime import datetime

from tqdm import tqdm  # if calling from command line or .py
#from tqdm.notebook import tqdm  # if calling from ipynb

import warnings
warnings.filterwarnings(action='ignore',message='invalid value encountered in reduce')
warnings.filterwarnings(action='ignore',message='Converting a CFTimeIndex with dates from a non-standard calendar')

try:
    paths = yaml.safe_load(open('datapaths.yaml','r'))
    modelvars = yaml.safe_load(open('modelvars.yaml','r'))
except FileNotFoundError:
    paths = yaml.safe_load(open('../datapaths.yaml','r'))
    modelvars = yaml.safe_load(open('../modelvars.yaml','r'))


#-------------------------------------------------------------------------
#  Interpret input dates and generate list of years
#-------------------------------------------------------------------------

def interpret_t0tf(t0,tf):
    
    if type(t0)!=str: 
        years = np.arange(t0,tf+1)
        t0,tf = str(t0),str(tf)
    else: 
        years = np.arange(int(t0[0:4]),int(tf[0:4])+1)
    months = np.arange(1,13,1)
    
    return t0,tf,years,months


#-------------------------------------------------------------------------
#  Standardize grid: reorder to (-180,180,-90,90)
#-------------------------------------------------------------------------

def standardize_grid(ds): 
    
    # arguably should use (0,360) since that's the model default
    # but we hates it
    
    # if lon includes both 0 and 360, trim ds to exclude 360
    if (ds.lon[0]==0.)*(ds.lon[-1]==360.):
        ds = ds.sel(lon=slice(0,ds.lon[-2]))

    # if lon is on (0,360), rotate to (-180,180)
    if np.abs(ds.lon[0])<10:
        ds.coords['lon'] = (ds.coords['lon'] + 180.) % 360. - 180.
        ds = ds.sortby(ds.lon)
        
    # if lat is on (90,-90), reflect to (-90,90)
    if ds.lat[0]>0:
        ds.coords['lat'] = ds.coords['lat'][::-1]
        for var in ds.data_vars:
            ds[var].values = ds[var].values[:,::-1,:]

    return ds


#-------------------------------------------------------------------------
#  Standardize time coordinate: make into datetime if not
#-------------------------------------------------------------------------

def standardize_calendar(ds,years,months):
    
    # get a reference time value to assess
    d0 = ds.time.values[0]
        
    # want everything to be np.datetime64
    if type(d0)==np.datetime64:
        return ds
    
    # some models: time coord is just months (arbitrarily assign to y=yf)
    if np.array_equal(ds.time.values,months): 
        newtime = np.array([datetime(years[-1],m,15) for m in months])
        ds = ds.assign_coords({'time':newtime})
        
    # some models: time is YYYYMM as str/int/float
    elif isinstance(d0,(str,int,float)):
        if len(str(int(d0)))==6:
            dstrs = [str(date) for date in ds.time.values]
            newtime = np.array([datetime(int(dstr[0:4]),int(dstr[4:6]),15) for dstr in dstrs])
            ds = ds.assign_coords({'time':newtime})

    # otherwise it's probably some other datetime-esque format
    ds['time'] = ds.indexes['time'].to_datetimeindex()
    newtime = np.array([datetime(t.dt.year,t.dt.month,16,0) for t in ds.time])
    ds = ds.assign_coords({'time':newtime.astype('datetime64[ns]')})
    
    return ds


#-----------------------------------------------------------------------------
#  Build dummy ds's for y/m combos that don't have data
#-----------------------------------------------------------------------------

def extend_dslist(dslist,ymin,ymax,years,var):

    # ymin,ymax = years that have data available
    # years = years for which data have been requested
    # if there is no overlap, len(dslist)=0 and we build from scratch
    # (this is the only case in which we need var)

    if len(dslist)==0: 
        ref = xr.Dataset(data_vars={var:(['lat','lon'], np.ones((90,180))*np.nan)},
                         coords={'lat':np.arange(-90,90,2),'lon':np.arange(-180,180,2),
                                 'time':datetime(2000,1,15)})
    else: ref = dslist[0]

    if years[0]<ymin:
        prelist = [xr.full_like(ref,np.nan).assign_coords({'time':datetime(y,m,15)})
                   for y in np.arange(years[0],ymin+1) for m in np.arange(1,13)]
    else: prelist = []

    if years[-1]>ymax:
        postlist = [xr.full_like(ref,np.nan).assign_coords({'time':datetime(y,m,15)})
                   for y in np.arange(ymax,years[-1]+1) for m in np.arange(1,13)]
    else: postlist = []

    extended_dslist = prelist + dslist + postlist

    return extended_dslist


def build_single_dummy(ref,y,m):

    # in some cases I just need a single fill, so extend_dslist doesn't help
    dummy = xr.full_like(ref,np.nan)
    dummy = dummy.assign_coords({'time':datetime(y,m,15)})

    return dummy



#-----------------------------------------------------------------------------
#  Read Model Ensemble Data (assuming dir structure /model/experiment/var/)
#-----------------------------------------------------------------------------

def read_model_ensemble(model,expt,var,t0,tf):

    try: 
        dpath = paths['model_ensembles']+'%s/%s/%s/'%(model,expt,var)
        if not os.path.isdir(dpath): raise FileNotFoundError()
    except FileNotFoundError:
        localvar = modelvars[var]
        dpath = paths['model_ensembles']+'%s/%s/%s/'%(model,expt,localvar)
        if not os.path.isdir(dpath): 
            raise FileNotFoundError('path not found; tried %s and /%s/'%(dpath,var))
        
    flist = glob.glob(dpath+'*nc')
    t0,tf,years,months = interpret_t0tf(t0,tf)
    
    dslist = [standardize_calendar(xr.open_dataset(f),years,months) for f in flist]
        
    runs = np.array([ds.variant_label for ds in dslist])
    if len(np.unique(runs))!=len(runs):
        dslist = [xr.concat([dsi for ri,dsi in zip(runs,dslist) if ri==run],'time')
                  for run in np.unique(runs)]
                               
    ds = xr.concat(dslist,pd.Index(np.unique(runs),name='run'))
    for dsi in dslist: dsi.close()
    ds = ds.sel(time=slice(t0,tf))
                
    inv_keys = {v: k for k, v in modelvars.items() if v in list(ds.variables)}
    ds = ds.rename(inv_keys)
    ds = standardize_grid(ds)

    return ds


#-----------------------------------------------------------------------------
#  Read Model Data - Single Rlzn  (assuming dir structure /model/runid/<all files>)
#-----------------------------------------------------------------------------

def read_model_singlerlzn(model,runid,var,t0,tf):
    
    dpath = paths['model_singles']+'%s/%s/'%(model,runid)
    try: 
        ds = xr.open_mfdataset(glob.glob(dpath+'%s*nc'%var))
    except:
        ds = xr.open_mfdataset(glob.glob(dpath+'%s*nc'%modelvars[var]))
        ds = ds.rename({modelvars[var]:var})

    t0,tf,years,months = interpret_t0tf(t0,tf)
    ds = standardize_calendar(ds,years,months).sel(time=slice(t0,tf))
    inv_keys = {v: k for k, v in modelvars.items() if v in list(ds.variables)}
    ds = ds.rename(inv_keys)
    ds = standardize_grid(ds)
    
    return ds


#-----------------------------------------------------------------------------
#  Read MISR
#-----------------------------------------------------------------------------

def read_misr(var,t0,tf):
    
    print('reading misr %s'%var)

    dpath = paths['misr']
    ymin,ymax = 2001,2020
    t0,tf,years,months = interpret_t0tf(t0,tf)
    
    if var=='aod': longvar = 'Aerosol_Optical_Depth'
    elif var=='aaod': longvar = 'Absorbing_Optical_Depth'
    elif var=='ae': longvar = 'Angstrom_Exponent_550_860'

    dslist = []
    for y in years:
        for m in months: 
            if (y<ymin) or (y>ymax): continue
            fname = 'MISR_AM1_CGAS_%s_%d_F15_0032.nc'%(calendar.month_abbr[m].upper(),y)
            dsraw = xr.open_dataset(dpath+fname,group='Aerosol_Parameter_Average')
            vals = dsraw.sel(Band='green_558nm',Optical_Depth_Range='all')[longvar].values
            dsym = xr.Dataset(data_vars={var:(['lat','lon'],vals)},
                                         coords={'lat':dsraw.Latitude.values,
                              'lon':dsraw.Longitude.values,
                              'time':datetime(y,m,15)})
            dslist.extend([dsym])
            dsraw.close()

    dslist = extend_dslist(dslist,ymin,ymax,years,var)
    ds = xr.concat(dslist,dim='time').sel(time=slice(t0,tf))
    for dsi in dslist: dsi.close()
        
    ds = standardize_grid(ds)
    
    return ds




#-----------------------------------------------------------------------------
#  Read MODIS (comes with a dedicated function for opening files)
#-----------------------------------------------------------------------------


def read_modis(sat,var,t0,tf):
    
    print('reading MODIS-%s %s'%(sat,var))
    
    dslist = []
    dpath = paths['modis-%s'%sat.lower()]
    ymin= 2001 if sat=='Terra' else 2003
    ymax = 2020
    t0,tf,years,months = interpret_t0tf(t0,tf)
    
    for y in years:
        if (y<ymin) or (y>ymax): continue
        fpaths = sorted(glob.glob(dpath+'*M3.A%d*'%y))
        if len(fpaths)!=12: 
            print('found %d fpaths for %d'%(len(fpaths),y))
            print('searching for %s*M3.A%d*'%(dpath,y))
            print('paths found:',fpaths)
            raise FileNotFoundError
        for m,fpath in zip(months,fpaths): 
            dslist.extend([open_modis_file(y,m,fpath,var)])

    dslist = extend_dslist(dslist,ymin,ymax,years,var)
    ds = xr.concat(dslist,dim='time').sel(time=slice(t0,tf))
    for dsi in dslist: dsi.close()
        
    ds = standardize_grid(ds)
    
    return ds


def open_modis_file(y,m,fpath,var):

    if var=='aod':
        longvar = 'AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean'

    elif var=='ae':
        longvar = 'Deep_Blue_Angstrom_Exponent_Land_Mean_Mean'

    f = SD(fpath,SDC.READ)
    sds_obj = f.select(longvar)
    vals = sds_obj.get().astype(float)
    
    # restrict to w/in valid range (set outside vals to nan)
    # handles missing data at the same time (fill=-9999)
    valid_range = sds_obj.attributes()['valid_range']
    vals[(vals<valid_range[0])|(vals>valid_range[1])] = np.nan

    # unpack with scale_factor and add_offset
    scale_factor = sds_obj.attributes()['scale_factor']
    add_offset = sds_obj.attributes()['add_offset']
    vals = (vals - add_offset) * scale_factor

    # get lats (89.5,-89.5,1) and lons (-179.5,179.5,1)
    lats = f.select('YDim').get().astype(float)
    lons = f.select('XDim').get().astype(float)

    # make mini ds
    dsym = xr.Dataset(data_vars={var:(['lat','lon'],vals)},
                      coords={'lat':lats,'lon':lons,'time':datetime(y,m,15)})
    
    return dsym
    




#-----------------------------------------------------------------------------
#  Read MIDAS 1x1 degree product
#-----------------------------------------------------------------------------

def read_midas1x1(t0,tf):
    
    print('reading midas dod at 1x1deg')

    dpath = paths['midas1x1']
    ymin,ymax = 2015,2020
    t0,tf,years,months = interpret_t0tf(t0,tf)
    
    dslist = []
    fbase = 'MODIS-AQUA_AOD-and-DOD-GRID_RESOLUTION_1.0'
    
    for y in years:
        for m in months: 

            if (y<ymin) or (y>ymax): continue
        
            flist = glob.glob(dpath+'%d/%s-%d%02d*nc'%(y,fbase,y,m))
            dsi = xr.open_mfdataset(flist,drop_variables=['Latitude','Longitude']).mean('Time')
        
            nc = ncDataset(flist[0])
            coords = {'Longitude':nc['Longitude'][0,:].data,
                      'Latitude':nc['Latitude'][:,0].data,
                      'time':np.datetime64('%d-%02d-15'%(y,m))}
            
            dsi = dsi.assign_coords(coords)[['Modis-total-dust-optical-depth-at-550nm']]
            dsi = dsi.rename({'Latitude':'lat','Longitude':'lon',
                              'Modis-total-dust-optical-depth-at-550nm':'dod'})
            dslist.extend([dsi])

    dslist = extend_dslist(dslist,ymin,ymax,years,'dod')
    ds = xr.concat(dslist,dim='time').sel(time=slice(t0,tf))
    for dsi in dslist: dsi.close()
    ds = standardize_grid(ds)
    
    return ds





#-----------------------------------------------------------------------------
#  Read MIDAS 0.1x0.1 degree product
#-----------------------------------------------------------------------------

def read_midas01x01(t0,tf):
    
    print('reading midas dod at 0.1x0.1deg')
    
    dpath = paths['midas01x01']
    ymin,ymax = 2003,2017
    t0,tf,years,months = interpret_t0tf(t0,tf)
    fbase = 'MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-'

    dslist = []
    for y in years:
        ddir = dpath+'%d/GRID_RESOLUTION_0.1/'%y
        for m in months: 
            
            if (y<ymin) or (y>ymax): continue

            flist = glob.glob(ddir+fbase+'%d%02d*nc'%(y,m))
            dsi = xr.open_mfdataset(flist,drop_variables=['Latitude','Longitude']).mean('Time')
    
            nc = ncDataset(flist[0])
            coords = {'Longitude':nc['Longitude'][0,:].data,
                      'Latitude':nc['Latitude'][:,0].data,
                      'time':np.datetime64('%d-%02d-15'%(y,m))}
        
            dsi = dsi.assign_coords(coords)[['Modis-total-dust-optical-depth-at-550nm']]
            dsi = dsi.rename({'Latitude':'lat','Longitude':'lon',
                              'Modis-total-dust-optical-depth-at-550nm':'dod'})
            dslist.extend([dsi])

    dslist = extend_dslist(dslist,ymin,ymax,years,'dod')
    ds = xr.concat(dslist,dim='time')
    for dsi in dslist: dsi.close()
    ds = ds.sel(time=slice(t0,tf))
    ds = standardize_grid(ds)
        
    return ds





#-----------------------------------------------------------------------------
#  Read CALIOP (comes with a dedicated function for opening files)
#-----------------------------------------------------------------------------

def read_caliop(sky,daynight,var,t0,tf):
    
    print('reading caliop %s %s %s'%(sky,daynight,var))

    ymin,ymax = 2007,2020
    t0,tf,years,months = interpret_t0tf(t0,tf)
    dslist = []
    
    if len(tf)==7: maxmonth = int(tf[5:])
    elif years[-1]==2020: maxmonth = 7 # latest data avail
    else: maxmonth = 12

    for y in years: 
        if (y<ymin) or (y>ymax): continue
        for m in months: 
            if (y==2016)&(m==2):  
                dslist.append(build_single_dummy(dslist[-1],y,m))
            elif (y==2019)&(m==11)*(sky=='AllSky')*(daynight=='Day'):
                dslist.append(build_single_dummy(dslist[-1],y,m))
            elif (y==years[-1])&(m>=maxmonth): 
                dslist.append(build_single_dummy(dslist[-1],y,m))
            else: 
                dslist.append(open_caliop_file(sky,daynight,var,y,m))
                
    dslist = extend_dslist(dslist,ymin,ymax,years,var)
    ds = xr.concat(dslist,dim='time').sel(time=slice(t0,tf))
    for dsi in dslist: dsi.close()
        
    ds = standardize_grid(ds)
    
    return ds


def open_caliop_file(sky,daynight,var,y,m):
    
    dpath = paths['caliop']+'Aerosol_Tropos_%s_%s/'%(sky,daynight)
    prefix = 'CAL_LID_L3_Tropospheric_APro'
    fname = prefix+'_%s-Standard-V4-20.%d-%02d%s.hdf'%(sky,y,m,daynight[0])

    longvar = {'aod':'AOD_Mean','dod':'AOD_Mean_Dust'}[var]
    
    try: 
        f = SD(dpath + fname, SDC.READ)
    except: 
        print('\nno file found for %d %02d'%(y,m))
        f = SD(dpath + fname, SDC.READ)
        
    vals = f.select(longvar).get().astype(float)
    valid_range = f.select(longvar).attributes()['valid_range'].split('...')
    vals[vals<float(valid_range[0])] = np.nan
    vals[vals>float(valid_range[1])] = np.nan

    dsym = xr.Dataset(data_vars={var:(['lat','lon'],vals)},
                      coords={'lat':f.select('Latitude_Midpoint').get().astype(float)[0],
                              'lon':f.select('Longitude_Midpoint').get().astype(float)[0],
                              'time':datetime(y,m,15)},
                      attrs={'wavelength':'532nm'})
    return dsym




###################################################################################
#  ACROS CALIPSO (not doing ACROS MODIS because it doesn't include 2020)
###################################################################################

def read_acros(var,t0,tf):
    
    dpath = paths['acros']
    ymin,ymax = 2007,2020
    t0,tf,years,months = interpret_t0tf(t0,tf)

    fcal1 = 'GlobalClearSkyMeanDustClimatology_CALIOP20072019Monthly.nc'
    fcal2 = 'GlobalClearSkyMeanDustClimatology_CALIOP2020Monthly.nc'    
    rawds = xr.open_mfdataset([dpath+fcal1,dpath+fcal2])
    rawds = rawds.rename({'clr_hlc_aeraod':'aod','clr_hlc_dustaod':'dod'})

    time = [datetime(y,m,16) for y in rawds['year'] for m in rawds['month']]
    datlist = [rawds.sel(year=y,month=m)[var].values 
               for y in rawds['year'] for m in rawds['month']]

    ds = xr.Dataset(data_vars={var:(['time','lat','lon'],np.array(datlist))},
                    coords={'time':time,'lat':rawds.lat,'lon':rawds.lon})
    
    # don't use extend_dslist() b/c working with a single ds with time dim
    if years[0]<ymin:
        prelist = [build_single_dummy(ds.sel(time='2008-01'),y,m)
                   for y in np.arange(years[0],ymin+1) for m in months]
        ds = xr.concat(prelist,[ds],dim='time')
    if years[-1]>ymax:
        postlist = [build_single_dummy(ds.sel(time='2008-01'),y,m)
                    for y in np.arange(ymax,years[-1]+1) for m in months]
        ds = xr.concat([ds],postlist,dim='time')

    ds = standardize_grid(ds).sel(time=slice(t0,tf))
    
    return ds

    

###################################################################################
#  POLDER
###################################################################################

def read_polder(var,t0,tf):

    print('reading polder %s'%var)

    dpath = paths['polder']
    ymin,ymax = 2005,2013
    t0,tf,years,months = interpret_t0tf(t0,tf)

    if var=='aod': longvar = 'AOD565'
    elif var=='aaod': longvar = 'AAOD565'
    elif var=='ssa': longvar = 'SSA565'
    elif var=='ae': longvar = 'AExp'

    dslist = []
    for y in years:
        if (y<ymin) or (y>ymax): continue
        for m in months: 
            fname = 'GRASP_POLDER_L3_%d%02d.1degree.nc'%(y,m)
            dsraw = xr.open_dataset(dpath+fname)[[longvar]]
            # shape (x,y=lon) but  x = -lat? need to reflect.
            vals = dsraw[longvar].values[::-1,:] 
            dsym = xr.Dataset(data_vars={var:(['lat','lon'],vals)},
                              coords={'lat':dsraw.x.values,
                                      'lon':dsraw.y.values,
                                      'time':datetime(y,m,15)})
            dslist.extend([dsym])
            dsraw.close()

    dslist = extend_dslist(dslist,ymin,ymax,years,var)
    ds = xr.concat(dslist,dim='time').sel(time=slice(t0,tf))
    for dsi in dslist: dsi.close()
        
    ds = standardize_grid(ds)
    
    return ds

