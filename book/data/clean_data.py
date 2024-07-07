import numpy as np
import xarray as xr
import xesmf as xe

### clean and downsample seawater data
d18o_sw = xr.open_dataset('/Users/starr/Downloads/D18O_Breitkreuz_et_al_2018.nc')
d18o_sw = d18o_sw[['depth_center','lat_1deg_center','lon_1deg_center','D18O_1deg','THETA_1deg','SALT_1deg']]

d18o_new = xr.Dataset(
    data_vars=dict(d18Osw=(["month","depth","lat","lon"],d18o_sw['D18O_1deg'].data),
        T=(["month","depth","lat","lon"],d18o_sw['THETA_1deg'].data),
    ),
    coords=dict(lon=d18o_sw['lon_1deg_center'][0,::].data,
                lat=d18o_sw['lat_1deg_center'][::,0].data,
                month=np.arange(1,13,1),
                depth=d18o_sw['depth_center'].data*-1,
    ),
    attrs=dict(description="gridded climatology from Breitkruez et al., 2018 (https://doi.org/10.1029/2018JC014300)"),
)

d18o_new['d18Osw'] = d18o_new['d18Osw'].assign_attrs(units="Permille VSMOW",desciption="Gridded seawater d18O (1950-1980 climatology)")
d18o_new['T'] = d18o_new['T'].assign_attrs(units="deg C",desciption="Gridded seawater temperature (1950-1980 climatology)")
d18o_new = d18o_new.sel(depth=slice(0,800))

# regrid to 2 degrees
ds_out = xr.Dataset(
    {
        "lat": (["lat"], np.arange(-40,80,2.5), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(-90,40,2.5), {"units": "degrees_east"}),
    }
)
ds_in = d18o_new
regridder = xe.Regridder(ds_in, ds_out, "bilinear",periodic=True)
d18o_new = regridder(ds_in,keep_attrs=True)



mas = '/Users/starr/Downloads/RECCAP2_region_masks_all_v20221025.nc'
mask = xr.open_dataset(mas)
ds_in = mask['atlantic']
ds_out = d18o_new
regridder = xe.Regridder(ds_in,ds_out,'nearest_s2d')
basin = regridder(ds_in,keep_attrs=True)
msk = basin>0
d18o_new = d18o_new.where(msk,drop=True)
d18o_new.to_netcdf('../data/gridded_seawater_data.nc')


##### clean and downsample Plafommodel output
warm = xr.open_dataset('/Users/starr/My Drive/Files/Data/Model/PLAFOM/PLAFOM2.0_GLOBAL_MONTHLY_CONC_warm-waterPlankForamSpecies.nc')
temp = xr.open_dataset('/Users/starr/My Drive/Files/Data/Model/PLAFOM/PLAFOM2.0_GLOBAL_MONTHLY_CONC_temperate-waterPlankForamSpecies.nc')
cold = xr.open_dataset('/Users/starr/My Drive/Files/Data/Model/PLAFOM/PLAFOM2.0_GLOBAL_MONTHLY_CONC_cold-waterPlankForamSpecies.nc')

lon = warm['longitude']
lon = lon[0,::].data
lat = warm['latitude']
lat = lat[::,0].data
month = ['Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan']
depth = warm['ndep'].data / 100

plafom = xr.Dataset(
    data_vars=dict(GRuberW=(["month","depth","lat","lon"],warm['GRuberW'].data),
                   GBulloides=(["month","depth","lat","lon"],temp['GBulloides'].data),
                   NPachyderma=(["month","depth","lat","lon"],cold['NPachyderma'].data),
    ),
    coords=dict(lon=lon,
                lat=lat,
                month=month,
                depth=depth,
    ),
    attrs=dict(description="Global monthly concentration of the cold-water planktonic foraminifera species from Kretschmer et al., 2018 ( https://doi.org/10.5194/bg-15-4405-2018)"),
)

plafom['GRuberW'] = plafom['GRuberW'].assign_attrs(description='Monthly concentration of G. ruber white',units='mmol C/m3')
plafom['NPachyderma'] = plafom['NPachyderma'].assign_attrs(description='Monthly concentration of N. pachyderma',units='mmol C/m3')
plafom['GBulloides'] = plafom['GBulloides'].assign_attrs(description='Monthly concentration of G. bulloides',units='mmol C/m3')

# regrid to 2 degrees
ds_out = xr.Dataset(
    {
        "lat": (["lat"], np.arange(-40,80,2.5), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(-90,40,2.5), {"units": "degrees_east"}),
    }
)
ds_in = plafom
regridder = xe.Regridder(ds_in, ds_out, "bilinear",periodic=True)
plafom_re = regridder(ds_in,keep_attrs=True)


mas = '/Users/starr/Downloads/RECCAP2_region_masks_all_v20221025.nc'
mask = xr.open_dataset(mas)
ds_in = mask['atlantic']
ds_out = plafom_re
regridder = xe.Regridder(ds_in,ds_out,'nearest_s2d',periodic=True)
basin = regridder(ds_in,keep_attrs=True)
msk = basin>0
plafom_re = plafom_re.where(msk,drop=True)
plafom_re = plafom_re.sel(depth=slice(0,600))
plafom.to_netcdf('../data/plafom_foram_abundance.nc')