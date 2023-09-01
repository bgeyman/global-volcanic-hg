import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import xesmf as xe
import time

print('STARTING SCRIPT: 1-Create_Emission_Pulse_Files.py')

# make plume altitude mask
altitude = pd.read_csv('../Data/GEOS-Chem_72_Level_Reference.csv')
area = xr.open_mfdataset('../Data/grid_area_1x1.nc')

dummy_array = xr.DataArray(np.ones((72,180,360)), 
             dims=['lev','lat', 'lon'],
             coords={'lev': np.arange(1,73),
                     'lat': area['lat'].values,
                     'lon': area['lon'].values},)

dummy_array_2D = xr.DataArray(np.ones((180,360)), 
             dims=['lat', 'lon'],
             coords={'lat': area['lat'].values,
                     'lon': area['lon'].values},)

def make_one_day_point_emission(df=dummy_array, lev=1, lat=46.5, lon=-122.5,
                                start_date="2010/01/01", end_date="2010/01/04"):
    test_file1 = df.copy()
    test_file1 = test_file1.where(((test_file1['lon']==lon)), 0)
    test_file1 = test_file1.where(((test_file1['lat']==lat)), 0)
    if lev != None:
        test_file1 = test_file1.where(((test_file1['lev']==lev)), 0)
    else:
        print('lev set to None, skipping')
    
    area = xr.open_mfdataset('../Data/grid_area_1x1.nc')
    test_file1['area'] = area['area']
    
    time_list = pd.date_range(start_date, end_date, freq='D')
    test_file1 = test_file1.expand_dims(time=time_list)
    
    test_file1['lat'] = test_file1['lat'].assign_attrs({'units':'degrees_north',
                                     'long_name':'latitude', 'axis':'Y'})
    test_file1['lon'] = test_file1['lon'].assign_attrs({'units':'degrees_east',
                                     'long_name':'longitude', 'axis':'X'})
    if lev != None:
        test_file1['lev'] = test_file1['lev'].assign_attrs({
                                        'units':'level',
                                        'long_name':'GEOS-Chem level',
                                        'positive':'up', 'axis':'Z'})
    
    test_file1 = test_file1.to_dataset(name='Hg')
    
    test_file1 = test_file1.where(test_file1['time'] == time_list[0], 0)
    
    ## scale emission magnitude to 100 Mg in units of kg m-2 s-1
    pixel_area = (test_file1['Hg']*test_file1['area']).sum().values
    seconds_per_day = 86400
    emis_mass = 100*1000 # 10,000 kg = 100 Mg
    
    test_file1['Hg'] = (emis_mass*test_file1['Hg'])/(pixel_area*seconds_per_day)
    
    test_file1['Hg'] = test_file1['Hg'].assign_attrs({'long_name':'mercury', 'units':'kg m-2 s-1'})
    import datetime as dt
    import time
    from netCDF4 import date2num,num2date
    
    def convert_numpy_datetime64_to_datetime_seconds(in_date):
        dt64 = in_date
        unix_epoch = np.datetime64(0,'s')
        one_second = np.timedelta64(1, 's')
        seconds_since_epoch = (dt64 - unix_epoch) / one_second
        out_date = dt.datetime.utcfromtimestamp(seconds_since_epoch)
        return out_date

    def convert_numpy_datetime64_to_datetime_days(in_date):
        dt64 = in_date
        unix_epoch = np.datetime64(0,'D')
        one_day = np.timedelta64(1, 'D')
        days_since_epoch = (dt64 - unix_epoch) / one_day
        out_date = dt.datetime.utcfromtimestamp(days_since_epoch)
        return out_date

    dates = []
    for in_date in test_file1['time'].values:
        out_date = convert_numpy_datetime64_to_datetime_seconds(in_date)
        dates.append(out_date)
        
    days_since = start_date.split('/')[0]+'-'+start_date.split('/')[1]+'-'+start_date.split('/')[2]
    t_units = f'days since {days_since} 00:00:00'
    begin_date = start_date.split('/')[0]+start_date.split('/')[1]+start_date.split('/')[2] # ex: '20100101'
    times = date2num(dates, t_units)
    #print(dates)
    
    test_file1['time'] = times
    test_file1['time'] = test_file1['time'].assign_attrs({
        'units':t_units,
        'long_name':'time',
        'calendar':'standard',
        'delta_t':'0000-00-01 00:00:00',
        'begin_date':begin_date,
        'begin_time':'000000',
        'axis':'T'
    })
    
    return test_file1

def compress_netcdf(path_root, src_fn, trg_fn, delete_src=False):

    import netCDF4 as nc

    src = nc.Dataset(path_root+src_fn) # uncompressed file
    trg = nc.Dataset(path_root+trg_fn, mode='w') # compressed file

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions, zlib=True)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Save the file
    trg.close()
    src.close()
    
    # delete srcfile (uncompressed)
    if delete_src != False:
        # delete src line using os.
        import os
        os.remove(path_root+src_fn) 
        
    return

#--------------------------------
# make vertical emission fields
path_root = '../Output/Pulse_Emission_Fields/'

for lat, lon in [(46.5,-122.5), (15.5,120.5)]:
    print('lat',lat, ' | lon',lon)
    if ((lat==46.5) & (lon==-122.5)):
        volcano_name = 'MSH'
    elif ((lat==15.5) & (lon==120.5)):
        volcano_name = 'Pinatubo'
        
    for lev in [5, 10, 20, 30, 35, 45]:
        for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
            out = make_one_day_point_emission(df=dummy_array, lev=lev, lat=lat, lon=lon,
                                              start_date=f"2010/{month}/01", end_date=f"2010/{month}/04")
            # check that emission is within tolerance of 100 Mg (tolerance is 10 kg)
            if np.abs((out['Hg']*out['area']).sum().values*86400 - 100000) > 10:
                print('wrong value - does not make 100 Mg emission')

            out.to_netcdf(f'{path_root}{volcano_name}_L{lev}_mo{month}.nc')
            # compress file and delete original
            compress_netcdf(path_root = path_root, 
                           src_fn = f'{volcano_name}_L{lev}_mo{month}.nc', 
                           trg_fn = f'{volcano_name}_L{lev}_mo{month}_c.nc', 
                           delete_src=True)

print('COMPLETED SCRIPT: 1-Create_Emission_Pulse_Files.py')
