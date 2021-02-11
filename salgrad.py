'''
This script calculates surface horizontal salinity gradients for the TXLA model output 
in the nGOM from 1994-2017, subsampled from June-15 to July-25 each year to capture 
long-term variability for the sunrise cruise. 
'''
import numpy as np
import xgcm
from xgcm import Grid
from xhistogram.xarray import histogram
import xarray as xr
import xroms
import h5netcdf
from glob import glob 

#Need to subset the netcdf files because xroms will crash if I open them all. 

years = np.arange(1994, 2017)
for y in years:
    paths = glob('/d1/shared/TXLA_ROMS/output_20yr_obc/%i/ocean_his_00*.nc' % y)
    
    print('Loading Data')
    ds = xroms.open_mfnetcdf(paths, 
                             chunks = {'ocean_time':1})
    ds, grid = xroms.roms_dataset(ds, 
                                  Vtransform = None)
    print('Data Loaded')
    #This will yield a box in the SUNRISE area of interest: see histogram_demonstration for examples.
    etaslice = slice(30, 150)
    xislice = slice(270, 405)

#     dsdx = grid.derivative(ds.salt, 'X')
#     dsdy = grid.derivative(ds.salt, 'Y')
    dsdx = grid.interp(grid.diff(ds.salt, 'X'),'X', boundary = 'extend')/ds.dx
    dsdy = grid.interp(grid.diff(ds.salt, 'Y'),'Y', boundary = 'extend')/ds.dy
    
    dsdx.attrs = ''
    dsdy.attrs = ''
    
    dsdx_slice = dsdx.isel(eta_rho = etaslice, xi_rho = xislice, s_rho = -1).sel(ocean_time = slice(str(y)+'-06-15', str(y)+'-07-25'))
    dsdy_slice = dsdy.isel(eta_rho = etaslice, xi_rho = xislice, s_rho = -1).sel(ocean_time = slice(str(y)+'-06-15', str(y)+'-07-25'))

    sgradbins = np.linspace(-0.001,0.001,100)

    dsdx_slice.name = 'dsdx'
    dsdy_slice.name = 'dsdy'
    
#     #We need to remove the grid attributes from the variables. It will generate this annoying error code 
#     #because of an invalid character.
    dsdx_slice.attrs = ''
    dsdy_slice.attrs = ''
    
    print('Computing histogam')
    dsdx_hist = histogram(dsdx_slice, bins = [sgradbins], density = True)
    dsdy_hist = histogram(dsdy_slice, bins = [sgradbins], density = True)
    
    dsdx_hist.attrs = ''
    dsdy_hist.attrs = ''
    
    print('Saving histogram')
    dsdx_hist.to_netcdf('/home/dylan/Variability/histograms/dsdx/histograms_dsdx_%i.nc' % y)
    dsdy_hist.to_netcdf('/home/dylan/Variability/histograms/dsdy/histograms_dsdy_%i.nc' % y)