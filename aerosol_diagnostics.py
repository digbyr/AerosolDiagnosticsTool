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
from analysis_scripts.singlemodel_seasonalmaps import call_seasonal_maps
from analysis_scripts.singlemodel_seasonalmaps_2obs import call_seasonal_maps_2obs
from analysis_scripts.comparemodels_globalproperties import call_compare_global
#from analysis_scripts.comparemodels_landvsocn import call_compare_landvsocn
from analysis_scripts.landvsocntest import call_compare_landvsocn

import numpy as np
import xarray as xr
import yaml


#######################################################################
#  User defined inputs
#######################################################################

config = yaml.safe_load(open('config.yaml'))

figtag = config['figtag']
y0,yf = config['y0'],config['yf'] 
varlist = config['varlist']

ens_models = config['ens_models']
ens_expts = config['ens_expts']

ind_models = config['ind_models']
ind_runids = config['ind_runids']

models_to_map = config['models_to_map']
modeltags_to_map = config['modeltags_to_map']

print()

#######################################################################
#  Read in obs to use in all of the analysis
#######################################################################

if 'aod' in varlist:
    obsAOD = {'MODIS Aqua':readers.read_modis('Aqua','aod',y0,yf),
              'MISR':readers.read_misr('aod',y0,yf),
              'CALIOP':readers.read_caliop('AllSky','Night','aod',y0,yf),
              #'POLDER-GRASP':readers.read_polder('aod',y0,yf),
              'POLDER-GRASP':readers.read_polder('aod',2006,yf)}

if 'aaod' in varlist:
    obsAAOD = {'MISR':readers.read_misr('aaod',y0,yf),
               #'POLDER-GRASP':readers.read_polder('aaod',y0,yf),
               'POLDER-GRASP':readers.read_polder('aaod',2006,yf)}

if 'dod' in varlist:
    obsDOD = {'ACROS':readers.read_acros('dod',y0,yf),
              'MIDAS':readers.read_midas01x01(y0,yf)}

if 'ae' in varlist:
    obsAE = {'MISR':readers.read_misr('ae',y0,yf),
             #'POLDER-GRASP':readers.read_polder('ae',y0,yf),
             'POLDER-GRASP':readers.read_polder('ae',2006,yf)}
             

obsdicts = {var: eval('obs%s'%var.upper()) for var in varlist}


#######################################################################
#  Call figures
#######################################################################


for model,tag in zip(models_to_map,modeltags_to_map):
    if (model in ens_models)*(tag in ens_expts): simtype='ens'
    elif (model in ind_models)*(tag in ind_runids): simtype='ind'
    call_seasonal_maps_2obs(model,tag,simtype,varlist,obsdicts,y0,yf)

call_compare_global(figtag,ens_models,ens_expts,ind_models,ind_runids,
                    varlist,obsdicts,y0,yf)

call_compare_landvsocn(figtag,ens_models,ens_expts,ind_models,ind_runids,
                       varlist,obsdicts,y0,yf)

#######################################################################
#  Clean up
#######################################################################

for var in varlist:
    obsdat = obsdicts[var]
    for sat,ds in obsdat.items(): ds.close()
