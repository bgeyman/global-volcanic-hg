## Calculate Atmospheric Effects
import numpy as np
import xarray as xr

directory_path = '/Volumes/Expansion/Volcano_Hg/Publication_Data/'

def load_df(df_fn, blank_fn, remove_blank=True):
    ''' Loads concentration and deposition influence from volcanic emission scenario.
        Optional `remove_blank` argument subtracts background when TRUE'''
    
    df = xr.open_mfdataset(f'{directory_path}/{df_fn}')

    ## This section calculates the direct influence of primary anthropogenic emissions by subtracting 
    ## field with same background emissions and zero primary anthropogenic emissions
    if remove_blank==True:
        blank = xr.open_mfdataset(f'{directory_path}/{blank_fn}')
        # variables to calculate a difference for
        blank_vars = ['SpeciesConc_Hg0','SpeciesConcHg2','SpeciesConcHgP',
                      'Total_Hg_Wet_Dep','Total_Hg_Dry_Dep','Net_Ocean_Hg0_Uptake',
                      'FluxHg0fromOceanToAir','FluxHg0fromAirToOcean']    
    
        for v in blank_vars:
            tmp_attr = df[v].attrs
            tmp_attr['processing_notes'] = 'concentration and deposition from *PRIMARY* emissions; blank removed'
            df[v] = df[v]-blank[v]
            df[v] = df[v].assign_attrs(tmp_attr)
    
    # variables to merge with dataset, unchanged
    conserved_vars = ['AREA','Met_AD','Met_AIRDEN','Met_PRECTOT']
    df = xr.merge([df, blank[conserved_vars]])
            
    df['Total_Hg_Dep'] = df['Total_Hg_Wet_Dep']+df['Total_Hg_Dry_Dep']+df['Net_Ocean_Hg0_Uptake']
    df['Total_Hg_Dep'] = df['Total_Hg_Dep'].assign_attrs({'Deposition Type':'Wet + Dry + Net Ocean Uptake', 'units':'μg m-2 yr-1',})
        
    df['Total_Hg_Conc'] = (df['SpeciesConc_Hg0']+df['SpeciesConcHg2']+df['SpeciesConcHgP'])
    df['Total_Hg_Conc'] = df['Total_Hg_Conc'].assign_attrs({'Notes':'Sum of all species: Hg0 + Hg2 + HgP', 'units':'mol mol-1 dry',})
    
    df = df.assign_attrs({'general_processing_notes':'blank removed from concentration and emission values'})
    return df

def get_mixture(Hg0_endmember, HGCL2_endmember, HGCLP_endmember, 
                fraction_dict = {'f_Hg0':1.0, 'f_HGCL2':0.0, 'f_HGCLP':0.0}):
    ''' Constructs atmospheric Hg concentration and deposition resulting from a mixture of Hg emission species.
        Here, three emission endmembers (Hg0, HGCL2, HGCLP) are hardcoded as arguments. 
        Arguments:
            Hg0_endmember: xarray dataset
            HGCL2_endmember: xarray dataset
            HGCLP_endmember: xarray dataset
            fraction_dict: dictionary with 3 keys (f_Hg0, f_HGCL2, f_HGCLP). Each corresponding item is a float such that the three sum to 1. 
        Assertion statement confirms that endmember fractions in `fraction_dict ` sum to within 0.1% of 1. 
        '''
    # confirm that Hg species fractions sum to 1.
    f_Hg0   = fraction_dict['f_Hg0']
    f_HGCL2 = fraction_dict['f_HGCL2']
    f_HGCLP = fraction_dict['f_HGCLP']
    assert (np.abs(f_Hg0 + f_HGCL2 + f_HGCLP - 1.)<0.001)
    # -- 
    
    df = xr.Dataset()

    for v in ['Total_Hg_Dep', 'Total_Hg_Wet_Dep', 'Total_Hg_Dry_Dep','Net_Ocean_Hg0_Uptake',
              'SpeciesConc_Hg0', 'SpeciesConcHg2','SpeciesConcHgP', 'Total_Hg_Conc',
              'FluxHg0fromOceanToAir','FluxHg0fromAirToOcean']:
        df[v] = ( (Hg0_endmember[v]*f_Hg0)+(HGCL2_endmember[v]*f_HGCL2)+(HGCLP_endmember[v]*f_HGCLP) )

        attrs = Hg0_endmember[v].attrs
        df[v] = df[v].assign_attrs(attrs)
        
    for v in ['AREA','Met_AD','Met_PRECTOT','Met_AIRDEN']:
        if v in df.data_vars:
            df = df.drop(v)
            
    return df

def create_Von_Glasow_scenario(year_start:int, year_end:int, scaling_method:str,
                               fraction_dict = {'f_Hg0':0.69, 'f_HGCL2':0.23, 'f_HGCLP':0.08}):
    ''' A function to construct a volcanic Hg scenario from a collection of  '''
    blank_fn = f'BLANK_GC_aggregated_{year_start}-{year_end}.nc'
    Von_Glasow = xr.Dataset()

    for eruption_type in ['pvf','eff','exp']:
        Hg0_fn    = f'{eruption_type}_{scaling_method}_Hg0_GC_aggregated_{year_start}-{year_end}.nc'
        HGCL2_fn  = f'{eruption_type}_{scaling_method}_HGCL2_GC_aggregated_{year_start}-{year_end}.nc'
        HG2CLP_fn = f'{eruption_type}_{scaling_method}_HG2CLP_GC_aggregated_{year_start}-{year_end}.nc'

        Hg0_endmember    = load_df(Hg0_fn,   blank_fn, remove_blank=True)
        HGCL2_endmember  = load_df(HGCL2_fn, blank_fn, remove_blank=True)
        HG2CLP_endmember = load_df(HGCL2_fn, blank_fn, remove_blank=True)
        
        attrs_dict = {}
        for v in Hg0_endmember.data_vars:
            attrs_dict[v] = Hg0_endmember[v].attrs

        if list(Von_Glasow.data_vars) == []:
            Von_Glasow = get_mixture(Hg0_endmember, HGCL2_endmember, HG2CLP_endmember, fraction_dict)
        else:
            Von_Glasow += get_mixture(Hg0_endmember, HGCL2_endmember, HG2CLP_endmember, fraction_dict)

    for v in Von_Glasow.data_vars:
        Von_Glasow[v] = Von_Glasow[v].assign_attrs(attrs_dict[v])

    Von_Glasow = xr.merge((Von_Glasow, Hg0_endmember[['AREA','Met_AD','Met_PRECTOT','Met_AIRDEN']]))
    return Von_Glasow

def get_atmospheric_mass(df, conc_var='Total_Hg_Conc'):
    
    molar_mass_Hg      = 201.0 # g mol-1
    molar_mass_species = 201.0 # g mol-1
    cf = molar_mass_species/molar_mass_Hg # correction factor for historical issue in SpeciesConc diagnostics
    molar_mass_air     = 28.97 # g mol-1
    
    # sum over all species - values in mol/mol
    df['Hg_mass_atm'] = df[conc_var].copy()
    df['Hg_mass_atm'] = df['Hg_mass_atm']*(molar_mass_Hg/molar_mass_air)*cf # g Hg g-1 air
    # dry air mass (kg) --> dry air mass (g)
    df['g_air'] = 1e3*df['Met_AD']
    # Hg mass (grams) is [ (g Hg g-1 air) * (g air) ]
    df['Hg_mass_atm'] = df['Hg_mass_atm']*df['g_air']
    # convert to Mg Hg
    df['Hg_mass_atm'] = df['Hg_mass_atm']*1e-6 

    return df['Hg_mass_atm']

def convert_atm_concentrations(df, var_list:list, units_from:str, units_to:str):
    ''' 
        Description: converts atmospheric mercury concentrations for a given 
                     set of variables (with common units) in an xarray Dataset (`df`) 
        Arguments:
            - `df` : an xarray Dataset
            - `var_list` : (str or list) contains names of variables with common
                           units for which to perform unit conversion
            - `units_from` : (str) units to convert from
            - `units_to` : (str) units to convert to
        
        Returns: xarray Dataset passed in, with units converted for selected vars
        
        Notes: supported conversions are listed in `allowed_conversions` dictionary below

        Example usage: 
           df = convert_atm_concentrations(df, var_list=['SpeciesConc_Hg0', 'SpeciesConcHg2'], units_from='mol mol-1 dry', units_to='ng m-3')
    '''

    # currently supported conversions. build this out if helpful in the future.
    allowed_conversions = {'mol mol-1 dry': ['ng m-3', 'g Hg'],
                           'ng m-3': ['mol mol-1 dry'],
                           'g Hg':['mol mol-1']}

    # values useful for conversion
    molar_mass_Hg      = 201.0 # g mol-1
    molar_mass_species = 201.0 # g mol-1
    cf = molar_mass_species/molar_mass_Hg # correction factor for issue in SpeciesConc diagnostics
    molar_mass_air     = 28.97 # g mol-1

    # if `var` is passed as a string, convert it to a list with 
    # 1 item for consistency. This allows lists or strings to be 
    # passed for conversion
    if type(var_list)==str:
        var_list = [var_list]

    # record original attributes to prevent losing them
    original_attrs = {}
    ## check that DataArray has the correct input units
    for v in var_list:
        assert df[v].units == units_from
        original_attrs[v] = df[v].attrs

    if ( (units_from=='mol mol-1 dry') & (units_to=='ng m-3') ):
        for v in var_list:
            df[v] = df[v]*(molar_mass_Hg/molar_mass_air)*cf # g Hg g-1 air
            df[v] = df[v]*(1e3*df['Met_AIRDEN']) # (1000 g air kg-1 air)*(kg air m-3 air)
            df[v] = df[v]*1e9 # ng Hg m-3 air
            df[v].attrs['units'] = 'ng m-3'

    if ( (units_from=='ng m-3') & (units_to=='mol mol-1 dry') ):
        for v in var_list:
            df[v] = df[v]*1e-9 # ng Hg --> g Hg
            df[v] = df[v]/(1e3*df['Met_AIRDEN']) # m-3 air --> g-1 air
            df[v] = df[v]/((molar_mass_Hg/molar_mass_air)*cf) # g Hg g-1 air --> mol Hg mol-1 air (dry)
            df[v].attrs['units'] = 'mol mol-1 dry'

    if ((units_from=='mol mol-1 dry') & (units_to=='g Hg')):
        ''' converts from `mol Hg mol-1 dry air` to `g Hg` (per grid box)'''
        for v in var_list:
            df[v] = df[v]*(molar_mass_Hg/molar_mass_air)*cf # g Hg g-1 air
            df[v] = df[v]*(1e3*df['Met_AD']) # g Hg g-1 air --> g Hg (per grid cell)
            df[v].attrs['units'] = 'g Hg'

    if ((units_from=='g Hg') & (units_to=='mol mol-1 dry')):
        ''' converts `g Hg` (per grid box) back to `mol mol-1 dry` '''
        for v in var_list:
            df[v] = df[v]/(1e3*df['Met_AD'])
            df[v] = df[v]/((molar_mass_Hg/molar_mass_air)*cf)
            df[v].attrs['units'] = 'mol mol-1 dry'
    return df

def get_TGM(df):
    df['SpeciesConc_TGM'] = (df['SpeciesConc_Hg0']+
                             df['SpeciesConcHg2']+
                             df['SpeciesConcHgP'])

    # check that all species have the same units
    units = df['SpeciesConc_Hg0'].units
    assert (units == df['SpeciesConcHg2'].units == df['SpeciesConcHgP'].units)
    # get other attrs
    if df['SpeciesConc_Hg0'].long_name[0:28] == 'Dry mixing ratio of species ':
        long_name = 'Dry mixing ratio of species TGM'
    else:
        original_long_name = df['SpeciesConc_Hg0'].long_name
        print(f'long_name unknown -- original species name was {original_long_name}')
        long_name = 'unknown'
    averaging_method = df['SpeciesConc_Hg0'].averaging_method
    # now assign those units to TGM
    df['SpeciesConc_TGM'] = df['SpeciesConc_TGM'].assign_attrs({'units':units,
    'long_name':long_name, 'averaging_method': averaging_method})
    return df

def get_global_deposition_Mg(df):
    global_total_dep_Mg = (df['Total_Hg_Dep']*df['AREA']).sum().values.item()*1e-12
    return {'total_dep':global_total_dep_Mg, 'units':'Mg/y'}