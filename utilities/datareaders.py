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

config = yaml.safe_load(open('config_dirs_vars.yaml','r'))


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
        
    # if lat is on (90,90), reflect to (-90,90)
    if ds.lat[0]>0:
        ds.coords['lat'] = ds.coords['lat'][::-1]
        for var in ds.data_vars:
            ds[var].values = ds[var].values[:,::-1,:]

    return ds


#-------------------------------------------------------------------------
#  Standardize time coordinate: make into datetime if not
#-------------------------------------------------------------------------

def standardize_time(ds,years,months):
    
    # get a reference time value to assess
    d0 = ds.time.values[0]
    
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

    # space to add more exceptions if needed
    return ds

   
#-----------------------------------------------------------------------------
#  Read Model Ensemble Data (assuming dir structure /model/experiment/var/)
#-----------------------------------------------------------------------------

def read_ensemble(model,expt,var,t0,tf):

    print('reading %s-%s %s from %s to %s'%(model,expt,var,str(t0),str(tf)))
    print('progress bar = ensemble')
    
    try: 
        dpath = config['directory_path']['models']+'%s/%s/%s/'%(model,expt,var)
        if not os.path.isdir(dpath): raise FileNotFoundError()
    except FileNotFoundError:
        localvar = config['variable_name'][var]
        dpath = config['directory_path']['models']+'%s/%s/%s/'%(model,expt,localvar)
        if not os.path.isdir(dpath): 
            raise FileNotFoundError('path not found; tried %s and /%s/'%(dpath,var))
        
    flist = glob.glob(dpath+'*nc')
    t0,tf,years,months = interpret_t0tf(t0,tf)
    
    try: 
        dslist = [xr.open_dataset(f).sel(time=slice(t0,tf)) for f in tqdm(flist)]
    except:
        dslist = [standardize_time(xr.open_dataset(f),years,months) for f in tqdm(flist)]
        dslist = [ds.sel(time=slice(t0,tf)) for ds in dslist]
    
    if len(flist)>1:
        runs = np.array([ds.variant_label for ds in dslist])
        ds = xr.concat(dslist,pd.Index(runs,name='run'))
        for dsi in dslist: dsi.close()
    else: 
        ds = dslist[0]

    inv_keys = {v: k for k, v in config['variable_name'].items() if v in list(ds.variables)}
    ds = ds.rename(inv_keys)
    ds = standardize_grid(ds)

    return ds


#-----------------------------------------------------------------------------
#  Read Single Model File (path specified in config file)
#-----------------------------------------------------------------------------

def read_singlefile(model,t0,tf):
    
    dpath = config['single_files'][model]
    ds = xr.open_dataset(dpath)
    
    t0,tf,years,months = interpret_t0tf(t0,tf)
    ds = standardize_time(ds,years,months).sel(time=slice(t0,tf))

    inv_keys = {v: k for k, v in config['variable_name'].items() if v in list(ds.variables)}
    ds = ds.rename(inv_keys)
    ds = standardize_grid(ds)
    
    return ds


#-----------------------------------------------------------------------------
#  Read MISR
#-----------------------------------------------------------------------------

def read_misr(var,t0,tf):
    
    print('reading MISR %s from %s to %s'%(var.upper(),str(t0),str(tf)))
    print('progress bar = years')
    
    dpath = config['directory_path']['misr']
    t0,tf,years,months = interpret_t0tf(t0,tf)
    
    if var=='aod': longvar = 'Aerosol_Optical_Depth'
    elif var=='aaod': longvar = 'Absorbing_Optical_Depth'
    elif var=='ae': longvar = 'Angstrom_Exponent_550_860'

    dslist = []
    for y in tqdm(years):
        for m in months: 
            fname = 'MISR_AM1_CGAS_%s_%d_F15_0032.nc'%(calendar.month_abbr[m].upper(),y)
            dsraw = xr.open_dataset(dpath+fname,group='Aerosol_Parameter_Average')
            vals = dsraw.sel(Band='green_558nm',Optical_Depth_Range='all')[longvar].values
            dsym = xr.Dataset(data_vars={var:(['lat','lon'],vals)},
                              coords={'lat':dsraw.Latitude.values,
                                      'lon':dsraw.Longitude.values,
                                      'time':datetime(y,m,15)})
            dslist.extend([dsym])
            dsraw.close()

    ds = xr.concat(dslist,dim='time').sel(time=slice(t0,tf))
    for dsi in dslist: dsi.close()
        
    ds = standardize_grid(ds)
    
    return ds




#-----------------------------------------------------------------------------
#  Read MODIS (comes with a dedicated function for opening files)
#-----------------------------------------------------------------------------

def read_modis(sat,var,t0,tf):
    
    print('reading MODIS-%s %s from %s to %s'%(sat,var.upper(),str(t0),str(tf)))
    print('progress bar = years')
    
    t0,tf,years,months = interpret_t0tf(t0,tf)
    
    dslist = []
    dpath = config['directory_path']['modis-%s'%sat.lower()]

    for y in tqdm(years): 
        fpaths = sorted(glob.glob(dpath+'*M3.A%d*'%y))
        if len(fpaths)!=12: 
            print('found %d fpaths for %d'%(len(fpaths),y))
            print('searched for %s*M3.A%d*'%(dpath,y))
            print('paths found:\n',fpaths)
            raise FileNotFoundError
        for m,fpath in zip(months,fpaths):
            dslist.extend([open_modis_file(y,m,fpath,var)])

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
#  Read MIDAS 0.1x0.1 degree product
#-----------------------------------------------------------------------------

def read_midas(t0,tf):
    
    print('reading MIDAS DOD at 0.1x0.1deg from %s to %s'%(str(t0),str(tf)))
    print('progress bar = years')
    
    dpath = config['directory_path']['midas01x01']
    fbase = 'MODIS-AQUA-C061_AOD-and-DOD-V1-GRID_RESOLUTION_0.1-'
    t0,tf,years,months = interpret_t0tf(t0,tf)

    dslist = []
    for y in tqdm(years):

        ddir = dpath+'%d/GRID_RESOLUTION_0.1/'%y
        for m in months: 
            
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

    ds = xr.concat(dslist,dim='time')
    for dsi in dslist: dsi.close()
    if '-' in t0: ds = ds.sel(time=slice(t0,tf))
    ds = standardize_grid(ds)
        
    return ds





#-----------------------------------------------------------------------------
#  Read CALIOP (comes with a dedicated function for opening files)
#-----------------------------------------------------------------------------

def read_caliop(sky,daynight,var,t0,tf):
    
    print('reading CALIOP-%s-%s %s from %s to %s'%(sky,daynight,var,str(t0),str(tf)))
    print('progress bar = years')
    
    t0,tf,years,months = interpret_t0tf(t0,tf)
    dslist = []
    for y in tqdm(years): 
        for m in months:
            if (y==2016)&(m==2): dslist.append(caliop_dummy(dslist[-1],y,m,var))
            elif (y==2019)&(m==11)*(sky=='AllSky')*(daynight=='Day'):
                dslist.append(caliop_dummy(dslist[-1],y,m,var))
            else: dslist.append(open_caliop_file(sky,daynight,var,y,m))
                
    ds = xr.concat(dslist,dim='time').sel(time=slice(t0,tf))
    for dsi in dslist: dsi.close()
        
    ds = standardize_grid(ds)
    
    return ds


def open_caliop_file(sky,daynight,var,y,m):
    
    dpath = config['directory_path']['caliop']+'Aerosol_Tropos_%s_%s/'%(sky,daynight)
    prefix = 'CAL_LID_L3_Tropospheric_APro'
    fname = prefix+'_%s-Standard-V4-20.%d-%02d%s.hdf'%(sky,y,m,daynight[0])

    longvar = {'aod':'AOD_Mean','dod':'AOD_Mean_Dust'}[var]
    
    try: 
        f = SD(dpath + fname, SDC.READ)
    except: 
        print('no file found for %d %02d'%(y,m))
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


def caliop_dummy(ref,y,m,var):

    # create fill data for months that don't have any
    lats = ref.lat.values
    lons = ref.lon.values
    vals = np.ones((len(lats),len(lons)))*np.nan
    ds = xr.Dataset(data_vars={var:(['lat','lon'],vals)},
                    coords={'lat':lats,'lon':lons,
                            'time':datetime(y,m,15)})
    return ds



#-----------------------------------------------------------------------------
#  Read Voss dust product
#-----------------------------------------------------------------------------

def read_voss(t0,tf):
    
    print('reading Voss DOD from %s to %s'%(str(t0),str(tf)))
    
    dpath = config['directory_path']['voss']
    ds = xr.open_dataset(dpath+'Voss2020_DAOD_monthly_MODIS.nc')
    newtime = np.array([np.datetime64(datetime.fromordinal(datenum).strftime('%Y-%m-%d')) 
                        for datenum in ds.time.values.astype(int)])
    ds = ds.assign_coords({'time':newtime}).rename({'daod':'dod'})
    ds = ds.sel(time=slice(str(t0),str(tf)))
    ds = standardize_grid(ds)
    
    return ds












