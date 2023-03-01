#######################################################################
#
#  CONFIG FILE TO COORDINATE AEROSOL DIAGNOSTICS
#  
#  Use this file to specify which models you wish to evaluate
#  and what types of figures you wish to make. 
#  It then calls the various analysis functions.
#
#######################################################################

import utilities.datareaders as readers
from singlemodel_seasonalmaps import call_seasonal_maps
from comparemodels_globalproperties import call_compare_global

import numpy as np
import xarray as xr


#######################################################################
#  User defined setup
#######################################################################

figtag = 'bulkaero'

ens_models = ['CanESM5-1']
ens_expts = ['historical']

ind_models = ['CanAM-bulk']*2
ind_runids = ['v52rc0-amip','bulkaero2-aya']

y0,yf = 2003,2008
varlist = ['aod','aaod']

print()

#######################################################################
#  Read in obs to use in all of the analysis
#######################################################################

if 'aod' in varlist:
    obsAOD = {'MODIS Aqua':readers.read_modis('Aqua','aod',y0,yf),
              'MISR':readers.read_misr('aod',y0,yf),
              'CALIOP':readers.read_caliop('AllSky','Night','aod',y0,yf)}

if 'aaod' in varlist:
    obsAAOD = {'MISR':readers.read_misr('aaod',y0,yf)}

if 'ae' in varlist:
    obsAE = {'MODIS Aqua':readers.read_modis('Aqua','ae',y0,yf),
             'MISR':readers.read_misr('ae',y0,yf)}
             
obsdicts = {var: eval('obs%s'%var.upper()) for var in varlist}


#######################################################################
#  Call figures
#######################################################################

for model,runid in zip(ind_models,ind_runids):
    call_seasonal_maps(model,runid,varlist,obsdicts,y0,yf)

call_compare_global(figtag,ens_models,ens_expts,ind_models,ind_runids,
                    varlist,obsdicts,y0,yf)


#######################################################################
#  Clean up
#######################################################################

for var in varlist:
    obsdat = obsdicts[var]
    for sat,ds in obsdat.items(): ds.close()