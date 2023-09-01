import numpy as np
import pandas as pd
import xarray as xr

print('STARTING SCRIPT: 2-Merge_SO2_Datasets.py')

old_carn = False

def get_large_continuous_sources_v2(filepath='../Data/Catalogue_SO2_2022_complete.xls'):
    '''
    General Notes: this function loads data which supercedes that in `get_Carn_df()`
    '''
    df = pd.read_excel(filepath, sheet_name='Catalogue')

    df = df[(df['TYPE']=='Volcano')]
    
    # convert kt --> Mg (1kt = 1000 Mg)
    kt_to_Mg = 1e3

    out_dict = {}
    for site in df['NUMBER']:
        sel = df[df['NUMBER']==site]
        name = sel['NAME'].values.item()
        lat = sel['LATITUDE'].values.item()
        lon = sel['LONGITUDE'].values.item()
        country = sel['COUNTRY'].values.item()
        altitude = sel['ELEVATION'].values.item()
        AMF = sel['AMF'].values.item()
        out_dict[name] = {'Volcano':name, 'Country':country, 'Latitude':lat, 'Longitude':lon,
                          'Alt (m)':altitude,  'AMF':AMF,
                          'Year':[], 'SO2 emissions (Mg yr-1)':[], 'Emission uncertainty (1 sigma)':[]}

        for y in np.arange(2005, 2023):
            out_dict[name]['Year'].append(y)
            
            E     = sel[f'y{y}'].values.item()*kt_to_Mg # convert from kt to Mg SO2
            sigma = sel[f's{y}'].values.item()*kt_to_Mg # convert from kt to Mg SO2
            out_dict[name]['SO2 emissions (Mg yr-1)'].append(E)
            out_dict[name]['Emission uncertainty (1 sigma)'].append(sigma)

    df_out = pd.DataFrame()
    for site in out_dict.keys():
        sel = out_dict[site]
        n_yrs = len(sel['Year'])
        tmp_df = pd.DataFrame({'Volcano':np.repeat(sel['Volcano'], n_yrs),
                               'Country':np.repeat(sel['Country'], n_yrs),
                               'Lat':np.repeat(sel['Latitude'], n_yrs),
                               'Lon':np.repeat(sel['Longitude'], n_yrs),
                               'Alt (m)':np.repeat(sel['Alt (m)'], n_yrs),
                               'AMF':np.repeat(sel['AMF'], n_yrs),
                               'Year':sel['Year'],
                               'SO2 emissions (Mg yr-1)':sel['SO2 emissions (Mg yr-1)'],
                               'Emission uncertainty (1 sigma)':sel['Emission uncertainty (1 sigma)'],
                              })

        df_out = pd.concat((df_out, tmp_df))
    
    return df_out

def get_MSVOL_df(filepath='../Data/MSVOLSO2L4_v04-00-2022m0505.txt'):

    f=open(filepath,"r")
    lines=f.readlines()

    column_header = lines[47]
    data=[]
    for x in lines[48:]:
        data.append(x.split('\n')[0])
    f.close()

    separated=[]
    for x in data:
        # most data are tab separated
        row = x.split('\t')[0:12]
        if len(row)==12:
            separated.append(row)
        # catch entries which are whitespace separated
        else:
            row = x.split()[0:12]
            separated.append(row)

    columns=column_header.split('\t')[0:12]

    # create dataframe and assign type to each column
    df = pd.DataFrame(separated, columns=columns)

    numeric_cols = ['lat','lon','v_alt','yyyy', 'mm', 'dd','vei','p_alt_obs','p_alt_est','so2(kt)']
    str_cols = ['volcano','type']

    df[numeric_cols] = df[numeric_cols].astype(float)
    df[str_cols] = df[str_cols].astype(str)
    
    return df

def merge_datasets(all_eruptions_output_path='../Data/all_SO2.csv',):
    
    old_carn = False
    
    #print('getting pvf data from Fioletov et al. (2023)')
    PVF = get_large_continuous_sources_v2()
    # -----------------
    PVF['Volcano'] = PVF['Volcano'].str.replace('Miyake-jima','Miyakejima')
    PVF['Alt (m)'] = PVF['Alt (m)'].astype(float)
    PVF['P_alt'] = PVF['Alt (m)']*1e-3
    PVF['Type'] = 'pvf'

    # load eruptive volcanic flux (EVF)
    EVF = get_MSVOL_df(filepath='../Data/MSVOLSO2L4_v04-00-2022m0505.txt')
    EVF = EVF.rename(columns={'volcano':'Volcano', 'lat':'Lat', 'lon':'Lon', 'vei':'VEI', 'type':'Type'})
    # add a flag for whether plume altitude was observed (True if observed, False if estimated)
    EVF['f_alt_obs'] = EVF['p_alt_obs'] == -999
    # consolidate `p_alt_est` and `p_alt_obs` into single column
    EVF['p_alt'] = np.where((EVF.p_alt_obs == -999.), EVF.p_alt_est, EVF.p_alt_obs)
    EVF = EVF.drop(columns=['p_alt_est','p_alt_obs'])

    dates_obs = pd.date_range(start='1/1/2005', end='1/1/2022')
    dates_obs = pd.DataFrame({'Datetime':dates_obs})

    direct_obs = PVF.merge(dates_obs, left_on=PVF.Year, right_on=dates_obs.Datetime.dt.year)

    PVF_merge = direct_obs.copy()

    PVF_merge['V_alt'] = PVF_merge['P_alt']
    PVF_merge['f_alt_obs'] = True
    PVF_merge['VEI'] = np.nan
    PVF_merge['SO2 (Mg/day)'] = PVF_merge['SO2 emissions (Mg yr-1)']/365.25
    PVF_merge['1 sigma (Mg/day)'] = PVF_merge['Emission uncertainty (1 sigma)']/365.25

    PVF_merge = PVF_merge.drop(columns=['Alt (m)', 'SO2 emissions (Mg yr-1)', 'Emission uncertainty (1 sigma)','Year'])
 
    PVF_merge['Source'] = 'Fioletov et al. (2023)'

    EVF_merge = EVF.copy()
    ymd = EVF_merge['yyyy'].astype(int).astype(str)+(EVF_merge['mm'].astype(int).astype(str).str.zfill(2))+(EVF_merge['dd'].astype(int).astype(str).str.zfill(2))
    EVF_merge['Datetime'] = pd.to_datetime(ymd, format='%Y%m%d')

    EVF_merge = EVF_merge.rename(columns={'v_alt':'V_alt', 'p_alt':'P_alt'})
    EVF_merge['SO2 (Mg/day)'] = EVF_merge['so2(kt)']*1000
    EVF_merge['1 sigma (Mg/day)'] = np.nan
    EVF_merge['Source'] = 'Carn (2021)'
    EVF_merge['f_imputed'] = False
    EVF_merge['AMF'] = np.nan
    EVF_merge['Country'] = ''

    EVF_merge = EVF_merge.drop(columns=['yyyy','mm','dd','so2(kt)'])

    all_eruptions = pd.concat((PVF_merge, EVF_merge), axis=0)
    all_eruptions = all_eruptions.reset_index().drop(columns=['key_0','index'])
    for v in ['Volcano','Country','Type','Source']:
        all_eruptions[v] = all_eruptions[v].astype(str)

    col_order = ['Volcano','Country','Lat','Lon','Datetime','Type','V_alt', 'P_alt','SO2 (Mg/day)', '1 sigma (Mg/day)',
                 'VEI','AMF','f_imputed','f_alt_obs','Source']

    all_eruptions = all_eruptions[col_order]
    
    col_type = {'Volcano':str,
                'Country':str,
                'Lat':float,
                'Lon':float,
                'Type':str,
                'V_alt':float, 
                'P_alt':float,
                'SO2 (Mg/day)':float, 
                '1 sigma (Mg/day)':float,
                'f_imputed':bool,
                'f_alt_obs':bool,
                'Source':str}
    
    for c in col_type.keys():
        all_eruptions[c] = all_eruptions[c].astype(col_type[c])
    
    all_eruptions.to_csv(all_eruptions_output_path, index=False)
    return all_eruptions

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def make_netcdf(df, Type='exp', Name='Explosive', attrs={}):
    df = df[df['Type']==Type]
    data = {"lat": df['grid_lat'], 
            "lon": df['grid_lon'], 
            "time":df['Datetime'], 
            "SO2 (Mg/day)": df['SO2 (Mg/day)'], 
            "n sources":    df['n sources'],
            "P_alt_weighted": df['P_alt_weighted'], 
            "V_alt_weighted": df["V_alt_weighted"]}
    
    df_rows = pd.DataFrame(data).set_index(["time", "lon", "lat"])
    ds = xr.Dataset.from_dataframe(df_rows)
    ds[Name] = ds['SO2 (Mg/day)']
    ds = ds.drop(['SO2 (Mg/day)'])
    
    ds[Name] = ds[Name].assign_attrs(attrs)
    return ds

def update_netcdf_metadata(ds, name, lat, lon, nc_output_path=f'~/Downloads/exp.nc'):
    ds = ds.reindex({"lat": lat, "lon":lon, "time":pd.date_range(start='1/1/2005', end='1/1/2022')})
    
    # set plume altitude, volcano altitude, and # sources as non-dimensional coordinates
    ds.coords['plume_altitude'] = (['time','lon','lat'], ds['P_alt_weighted'].data)
    ds.coords['volcano_altitude'] = (['time','lon','lat'], ds['V_alt_weighted'].data)
    ds.coords['n_sources'] = (['time','lon','lat'], ds['n sources'].data)

    # add metadata to coordinates
    ds['time'] = ds['time'].assign_attrs({'long_name' :'Time', 'axis':'T', 'delta_t': '0000-00-01 00:00:00'})
    ds['lat'] = ds['lat'].assign_attrs({'long_name': 'Latitude', 'units': 'degrees_north', 'axis': 'Y'})
    ds['lon'] = ds['lon'].assign_attrs({'long_name': 'Longitude', 'units': 'degrees_east', 'axis': 'X'})
    ds['plume_altitude'] = ds['plume_altitude'].assign_attrs({'long_name': 'Plume emission height', 'units': 'km', 
                                            'Notes':'Where multiple sources contribute, plume height is an average weighted by SO2 emission magnitude'})
    ds['volcano_altitude'] = ds['volcano_altitude'].assign_attrs({'long_name': 'Volcano altitude', 'units': 'km', 
                                            'Notes':'Where multiple sources contribute, volcano height is an average weighted by SO2 emission magnitude'})
    ds['n_sources'] = ds['n_sources'].assign_attrs({'long_name': 'Number of volcanic sources in grid cell', 'units': ''})

    # drop extraneous variables
    ds = ds.drop(['n sources', 'P_alt_weighted', 'V_alt_weighted'])
    
    # set compression level with 'encoding' argument
    ds.to_netcdf(nc_output_path, engine="h5netcdf",
                 encoding=({name:{"zlib": True, "dtype": 'float32', "complevel": 9},
                            'lat':{"zlib": True, "dtype": 'float32', "complevel": 9},
                            'lon':{"zlib": True, "dtype": 'float32', "complevel": 9},
                            'plume_altitude':{"zlib": True, "dtype": 'float32', "complevel": 9},
                            'volcano_altitude':{"zlib": True, "dtype": 'float32', "complevel": 9},
                            'n_sources':{"zlib": True, "dtype": 'float32', "complevel": 9}}))
    return ds

def clean_SO2_names():
    # script to enforce consistent naming in SO2 dataset between
    # Carn et al. (2021) and Fioletov et al. (2023) datasets

    # requires that all_SO2.csv exists
    df = pd.read_csv('../Output/all_SO2.csv')

    # structure: {`existing name`:`new name`}
    renaming_dict = {'Barren Island':'Barren_Island',
                    'Asosan':'Aso',
                    'Mt. Etna':'Etna',
                    'Jebel at Tair':'Jebel_al_Tair',
                    'La Palma':'La_Palma',
                    'Lokon-Empung':'Lokon_Empung',
                    'San Cristobal':'San_Cristobal',
                    'Sangeang Api':'Sangeang_Api',
                    'Santa Ana':'Santa_Ana',
                    'Sarychev_Peak':'Sarychev',
                    'Shiveluch':'Sheveluch',
                    'Sierra Negra':'Sierra_Negra',
                    'Soufriere Hills':'Soufriere_Hills'}
    
    for key in renaming_dict.keys():
        df['Volcano'] = df['Volcano'].str.replace(key, renaming_dict[key], regex=True)

    df.to_csv('../Output/all_SO2.csv', index=False)

    return

def call_all(all_eruptions_csv_output_path='../Output/all_SO2.csv', 
             output_exp='~/Downloads/exp_SO2.nc', 
             output_eff='~/Downloads/eff_SO2.nc', 
             output_pvf='~/Downloads/pvf_SO2.nc',
             save_netcdf_output=False):

    df =  merge_datasets(all_eruptions_output_path=all_eruptions_csv_output_path)

    lon = np.linspace(-179.5,179.5,360)
    lon_res = lon[1]-lon[0]

    lat = np.linspace(-89.5,89.5,180)
    lat_res = lat[1]-lat[0]

    grid_reference = []
    for volcano in df['Volcano'].unique():

        # get nearest latitude center in grid
        grid_lat = find_nearest(lat, df[df['Volcano']==volcano]['Lat'].iloc[0])
        # get nearest latitude center in grid
        grid_lon = find_nearest(lon, df[df['Volcano']==volcano]['Lon'].iloc[0])    

        grid_reference.append([volcano, grid_lat, grid_lon])

    grid_reference = pd.DataFrame(grid_reference, columns=['Volcano','grid_lat','grid_lon'])

    df = df.merge(grid_reference, left_on='Volcano', right_on='Volcano')

    # Compute SO2-weighted average of P_alt for situations where multiple emission sources exist in one grid box
    df['P_alt_weight'] = df['P_alt']*df['SO2 (Mg/day)']
    df['V_alt_weight'] = df['V_alt']*df['SO2 (Mg/day)']
    df = df.groupby(by=['Datetime','grid_lat','grid_lon','Type'],as_index=False).agg({'SO2 (Mg/day)':'sum',
                                                                                      'P_alt_weight':'sum',
                                                                                      'V_alt_weight':'sum',
                                                                                      'Volcano':'count'})
    
    ## drop entries where SO2 (Mg/day) equal to or less than 0.
    df = df[df['SO2 (Mg/day)']>0.]
    
    df['P_alt_weighted'] = df['P_alt_weight']/df['SO2 (Mg/day)']
    df['V_alt_weighted'] = df['V_alt_weight']/df['SO2 (Mg/day)']
    # drop weighted plume and volcano altitude columns now that they're not useful
    df = df.drop(columns='P_alt_weight')
    df = df.drop(columns='V_alt_weight')

    df = df.rename(columns={'Volcano':'n sources'})

    if save_netcdf_output == True:
        exp = make_netcdf(df, Type='exp', Name='Explosive_SO2', 
                          attrs={'long_name': 'Explosive volcanic SO2 emission magnitude', 'units': 'Mg SO2 day-1'})
        eff = make_netcdf(df, Type='eff', Name='Effusive_SO2',
                          attrs={'long_name': 'Effusive volcanic SO2 emission magnitude', 'units': 'Mg SO2 day-1'})
        pvf = make_netcdf(df, Type='pvf', Name='Passive_SO2',
                          attrs={'long_name': 'Passive volcanic SO2 emission magnitude', 'units': 'Mg SO2 day-1'})
        
        exp = update_netcdf_metadata(ds=exp, name='Explosive_SO2', lat=lat, lon=lon, nc_output_path=output_exp)
        eff = update_netcdf_metadata(ds=eff, name='Effusive_SO2',  lat=lat, lon=lon, nc_output_path=output_eff)
        pvf = update_netcdf_metadata(ds=pvf, name='Passive_SO2',   lat=lat, lon=lon, nc_output_path=output_pvf)

    # update SO2 names to be consistent between Fioletov et al. (2023) and Carn et al. (2021)
    clean_SO2_names()

    return

## ---------------------
call_all()
## ---------------------

print('COMPLETED SCRIPT: 2-Merge_SO2_Datasets.py')
