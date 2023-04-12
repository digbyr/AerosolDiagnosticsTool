# functions to regrid a dataset to a specified resolution
# regridder_global: (-180,180,-90,90) at default of 1x1
# will also do functions that accept specified lat/lon ranges


import numpy as np
import xarray as xr
import xesmf as xe

def regridder_global(ds,res=1):

    ds_target = xr.Dataset({'lon':(['lon'],np.arange(-180,180,res)),
                            'lat':(['lat'],np.arange(-90,91,res))})
    regridder = xe.Regridder(ds,ds_target,'bilinear',periodic=True)
    regridded = regridder(ds)

    return regridded




