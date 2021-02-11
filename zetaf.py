'''
This script calculates surface vertical vorticity for the TXLA model output in the nGOM from 1994-2017, subsampled 
from June-15 to July-25 each year to capture long-term variability for the sunrise cruise. 
'''
import numpy as np
import xgcm
from xgcm import Grid
from xhistogram.xarray import histogram
import xarray as xr
import xroms
from glob import glob 

#Need to subset the netcdf files because xroms will crash if I open them all. 

years = np.arange(2007, 2017)
for y in years:
    paths = glob('/d1/shared/TXLA_ROMS/output_20yr_obc/%i/ocean_his_00*.nc' % y)
    
    print('Loading Data')
    ds = xroms.open_mfnetcdf(paths, 
                             chunks = {'ocean_time':1})
    ds, grid = xroms.roms_dataset(ds, 
                                  Vtransform = None)
    print('Data Loaded')
    #This will yield a box in the SUNRISE area of interest
    etaslice = slice(30, 150)
    xislice = slice(270, 405)

    #Compute vertical relative vorticity
    rv = xroms.relative_vorticity(ds.u, ds.v, ds.u.attrs['grid'])
    #Interpolate to the rho points. 
    rv = grid.interp(rv, 'Z')
    rv = rv.isel(eta_v = etaslice, xi_u = xislice, s_rho = -1)

    #Note we're going to need to select the second to last vertical value since its the w-points

    fx = grid.interp(ds.f, 'X', boundary = 'extend')
    fxy = grid.interp(fx, 'Y', boundary = 'extend')
    f = fxy.isel(eta_v = etaslice, xi_u = xislice)

    RVn = rv/f

    zetabins = np.linspace(-2.5,2.5,100)

    RVn.name = 'relative_vorticity_n'
    RVn_slice = RVn.sel(ocean_time = slice(str(y)+'-06-15', str(y)+'-07-25'))
    #We need to remove the grid attributes from the variables. It will generate this annoying error code 
    #because of an invalid character.
    RVn_slice.attrs = ''
    
    print('Computing histogam')
    zetaf_hist = histogram(RVn_slice, bins = [zetabins], density = True)
    
    print('Saving histogram')
    zetaf_hist.to_netcdf('/home/dylan/Variability/histograms/rho/histograms_zetaf_%i.nc' % y)
