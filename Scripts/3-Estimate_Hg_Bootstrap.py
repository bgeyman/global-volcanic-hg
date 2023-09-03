import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import json
from scipy.stats import linregress
from numpy.random import default_rng

# set seed for reproducible output
seed = 9
np.random.seed(seed)
rng = default_rng(seed=seed)

print('STARTING SCRIPT: 3-Estimate_Hg_Bootstrap.py')

# ------------------------------------------------------------------------------------------------------------------
n_iterations = 10000
output_path = '../Output/Data/Bootstrap_Output/'
F_exclude_pre_1990_ratios = True

# -- call flags
flag_call_bootstrap_annual        = True
# condition for whether to run individual bootstrap for all volcanoes in SO2 dataset.
# this is set to `False` by default because it takes hours to run. 
flag_call_bootstrap_volcano_list  = False 

# read unique volcano names (saved from all_SO2.csv)
with open('../Data/volcano_list.json', 'r') as f:
  volcano_list = json.load(f)

print(f"calculating bootstrap diagnostics with {n_iterations} iterations")

# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
# define functions used in bootstrap estimation procedure
# ------------------------------------------------------------------------------------------------------------------

def choose(array):
    ''' Randomly select a single value from an array '''
    choice = rng.choice(array)
    return choice

def select_SO2_uniform(central_estimate, uncertainty=0.50):
    ''' Randomly select an SO2 emission magnitude from a uniform distribution bounded as +/- fraction given by uncertainty.
        Used for EXP and EFF emissions.'''
    scale =  np.random.uniform(low=(1-uncertainty), high=(1+uncertainty), size=None)
    SO2_emission = central_estimate*scale
    return SO2_emission

def select_SO2_normal(mu, sigma):
    ''' Randomly select an SO2 emission magnitude from within a normal distribution given by mean and standard deviation.
        Used for PVF emissions.'''
    SO2_emission = np.random.normal(mu, sigma, size=None)
    return SO2_emission

def select_SO2(SO2, columns_to_keep=['Volcano','Type','year','Lat', 'Lon', 'VEI', 'SO2 (Mg/a)','1 sigma (Mg/a)', 'HgSO2_ratio']):
    '''
        Randomly gnerate SO2 emission estimates for all volcanic emissions in database.
        
        Arguments:
             - SO2 : the original dataset of SO2 emission estimates from Fioletov et al. (2022) and Carn et al. (2021)
        
        External Function Calls:
            - subset_df()
            - select_SO2_uniform()
            - select_SO2_normal()
    '''
    out = pd.DataFrame()
    
    # subset of columns to be used
    SO2 = SO2[columns_to_keep]
    
    # sample EXP SO2 emissions from uniform distribution
    # VEI >=4, uncertainty 20% -- n=31
    exp = subset_df(df=SO2, column='Type', match='exp')
    exp_large = exp[exp['VEI']>=4]
    new_SO2_exp = select_SO2_uniform(central_estimate=exp_large['SO2 (Mg/a)'].values, uncertainty=0.20)
    tmp = pd.DataFrame({'Volcano':   exp_large['Volcano'].values,
                        'Type':      exp_large['Type'].values,
                        'year':      exp_large['year'].values,
                        'Lat':       exp_large['Lat'].values,
                        'Lon':       exp_large['Lon'].values,
                        'SO2 (Mg/a)':new_SO2_exp,
                        'HgSO2_ratio': exp_large['HgSO2_ratio'].values})
    out = pd.concat((out, tmp))
    # else (smaller eruptions), uncertainty is 50%
    exp_small = exp[exp['VEI']<4]
    new_SO2_exp = select_SO2_uniform(central_estimate=exp_small['SO2 (Mg/a)'].values, uncertainty=0.50)
    tmp = pd.DataFrame({'Volcano':   exp_small['Volcano'].values,
                        'Type':      exp_small['Type'].values,
                        'year':      exp_small['year'].values,
                        'Lat':       exp_small['Lat'].values,
                        'Lon':       exp_small['Lon'].values,
                        'SO2 (Mg/a)':new_SO2_exp,
                        'HgSO2_ratio': exp_small['HgSO2_ratio'].values})
    out = pd.concat((out, tmp))
    
    # sample EFF SO2 emissions from uniform distribution
    eff = subset_df(df=SO2, column='Type', match='eff')
    new_SO2_eff = select_SO2_uniform(central_estimate=eff['SO2 (Mg/a)'].values, uncertainty=0.50)
    tmp = pd.DataFrame({'Volcano':   eff['Volcano'].values,
                        'Type':      eff['Type'].values,
                        'year':      eff['year'].values,
                        'Lat':       eff['Lat'].values,
                        'Lon':       eff['Lon'].values,
                        'SO2 (Mg/a)':new_SO2_eff,
                        'HgSO2_ratio': eff['HgSO2_ratio'].values})
    out = pd.concat((out, tmp))
    
    # sample PVF SO2 emissions from normal distribution
    pvf = subset_df(df=SO2, column='Type', match='pvf')
    
    new_SO2_pvf = select_SO2_normal(mu=pvf['SO2 (Mg/a)'].values, sigma=pvf['1 sigma (Mg/a)'].values)
    new_SO2_pvf = np.where(new_SO2_pvf<0, 0, new_SO2_pvf)# prevent negative values
    tmp = pd.DataFrame({'Volcano':   pvf['Volcano'].values,
                        'Type':      pvf['Type'].values,
                        'year':      pvf['year'].values,
                        'Lat':       pvf['Lat'].values,
                        'Lon':       pvf['Lon'].values,
                        'SO2 (Mg/a)':new_SO2_pvf,
                        'HgSO2_ratio': pvf['HgSO2_ratio'].values})
    out = pd.concat((out, tmp))
    
    return out

def subset_df(df, column=str, match=str, condition='=='):
    ''' Subset pandas dataframe based on column, match_value, and condition.
    
        This function simply allows me to preserve verbose naming of dataframes and 
        columns while maintaining *relatively* readable code
    '''
    if condition=='==':
        df = df[df[column]==match]
    else:
        print('no match found for argument `condition`; returning original df')
    return df

def update_Hg_ratio_from_range(Hg_ratios):
    # note - this line is necessary for case in which pre_1990 values
    # are excluded because index gets messed up from dropped values
    # -- this iterating through pandas rows works but is bad form.
    Hg_ratios = Hg_ratios.reset_index(drop=True)
    # --
    for i in range(len(Hg_ratios)):
        sel = Hg_ratios.iloc[i]
        if sel['Hg/SO2 range flag'] == True:
            lo = sel['Hg/SO2 lower bound']
            hi = sel['Hg/SO2 upper bound']
            # -- here, I draw from a log10-uniform distribution 
            new_HgSO2 =  np.random.uniform(low=np.log10(lo), high=np.log10(hi), size=None)
            new_HgSO2 = 10**new_HgSO2
            # --- 
            Hg_ratios.at[i, 'Hg/SO2'] = new_HgSO2
    
    return Hg_ratios

def sample_measured_ratios(Hg_ratios):
    ''' sample a ratio for each volcanic system in Hg_ratio '''
    measured_ratios = {}
    for v in Hg_ratios['Volcano'].unique():
        sel = subset_df(df=Hg_ratios, column='Volcano', match=v)
        choice = choose(sel['Hg/SO2'].values)
        measured_ratios[v] = choice
    return measured_ratios

def make_matches(measured_ratios={}):
    '''
    Description: assigns an Hg:SO2 measurement to all volcanoes in SO2 dataset. 
                 Accepts a dictionary containing a single Hg:SO2 ratio to assign
                 to every measured volcanoes. For each unmeasured volcano,
                 function randomly chooses an Hg:SO2 measurement from among values
                 in dictionary. 
    
    Arguments: `measured_ratios` is a dictionary with following structure:
                keys: volcano name (str), values: Hg:SO2 ratios (float)
    '''
    matches = {'SO2_key':[], 'HgSO2_key':[], 'HgSO2_ratio':[]}

    # make an array of all of the measured ratios drawn (and passed) in this iteration
    # use this to randomly sample a ratio for each unmeasured site
    arr_measured_ratios = [measured_ratios[key] for key in measured_ratios.keys()]
    
    for v in SO2['Volcano'].unique():
        sel = subset_df(conv_keys, column='SO2_key', match=v)
        if len(sel['SO2_key'].values) > 0:
            v_name_HgSO2 = sel['HgSO2_key'].values.item()

            if len(sel) == 0:
                print(f'{v} not in keys') # debugging -- change conv_keys to include all volcanoes in SO2
                R_HgSO2 =  choose(arr_measured_ratios) #unmatched_ratio
            elif v_name_HgSO2 in list(measured_ratios.keys()):
                measured_ratio = measured_ratios[v_name_HgSO2]
            else:
                R_HgSO2 = choose(arr_measured_ratios) #unmatched_ratio

            matches['SO2_key'].append(v)
            matches['HgSO2_key'].append(v_name_HgSO2)
            matches['HgSO2_ratio'].append(R_HgSO2)

    # convert matches to DataFrame to merge with SO2
    matches = pd.DataFrame(matches)

    return matches

def draw_sample(SO2, Hg_ratios):
    # bootstrap 0
    ## sample Hg:SO2 from within reported range
    Hg_ratios = update_Hg_ratio_from_range(Hg_ratios)
    
    # bootstrap 1
    ## sample for each matched volcano
    measured_ratios = sample_measured_ratios(Hg_ratios)
    
    # bootstrap 2
    ## sample ratio from among unmatched volcanoes
    all_measured_ratios = [measured_ratios[v] for v in measured_ratios.keys()]

    ## assign an Hg:SO2 ratio to all unique volcanoes in SO2 dataset
    matches = make_matches(measured_ratios)

    ## merge Hg:SO2 ratios with SO2 magnitudes
    df = pd.merge(left=SO2, right=matches, left_on='Volcano', right_on='SO2_key')

    # bootstrap 3
    # select SO2 emission magnitude from within uncertainty for each (volcano, type, year)
    df = select_SO2(df)
    
    # calculate Hg emission magnitude
    df['Hg (Mg/a)'] = df['SO2 (Mg/a)']*df['HgSO2_ratio']

    return df

# ------------------------------------------------------------------------------------------------------------------
# load all of the base data 
# ------------------------------------------------------------------------------------------------------------------
def load_base_data(F_exclude_pre_1990_ratios=bool):
    # row-format merge of daily passive, effusive, explosive emissions from Fioletov et al. (2023) and MSVOLSO4
    SO2 = pd.read_csv('../Output/all_SO2.csv')
    SO2['Datetime'] = pd.to_datetime(SO2['Datetime'])
    SO2['year'] = SO2.Datetime.dt.year
    SO2 = SO2.groupby(by=['Volcano', 'Type', 'year'], as_index=False).agg({'Volcano':'first',
                                                                         'Country':  'first',
                                                                         'Lat':      'first',
                                                                         'Lon':      'first',
                                                                         'Type':     'first',
                                                                         'V_alt':    'first',
                                                                         'P_alt':    'first',
                                                                         'SO2 (Mg/day)':'sum',
                                                                         '1 sigma (Mg/day)':'sum',
                                                                         'VEI':   'mean',
                                                                         'AMF':   'mean',
                                                                         'Source':'first',
                                                                         'year':  'mean'})
    SO2 = SO2[SO2['year']<2022] # remove 2023 values

    SO2 = SO2.rename(columns={'SO2 (Mg/day)':'SO2 (Mg/a)', '1 sigma (Mg/day)':'1 sigma (Mg/a)'})

    # Hg:SO2 ratios from literature compilation
    Hg_ratios = pd.read_csv('../Data/Hg_SO2_ratios.csv')
    cols = ['Volcano','Method','Hg/SO2','Hg/SO2 range flag', 'Hg/SO2 lower bound', 'Hg/SO2 upper bound','Measured_Pre_1990','Reference']
    Hg_ratios = Hg_ratios[cols]

    Hg_ratios = Hg_ratios.dropna(subset=['Hg/SO2'], how='all')

    if F_exclude_pre_1990_ratios == True:
        Hg_ratios = Hg_ratios[Hg_ratios['Measured_Pre_1990']==False]

    # Dataframe matching volcano names in df to those in Hg_ratios
    conv_keys = pd.read_csv('../Data/conversion_key_SO2_and_HgSO2.csv')
    
    return SO2, Hg_ratios, conv_keys

# ------------------------------------------------------------------------------------------------------------------
# Call sequence -- Estimate Annual Total Hg emission magnitude for BOOTSTRAP scaling method
# ------------------------------------------------------------------------------------------------------------------
def call_bootstrap_annual(n_draws: int, output_path:str, F_exclude_pre_1990_ratios=False):

    if F_exclude_pre_1990_ratios == True:
        tag_1990 = 'excl_pre_1990'
    else:
        tag_1990 = 'incl_pre_1990'

    # ---- load input data
    SO2, Hg_ratios, conv_keys = load_base_data(F_exclude_pre_1990_ratios)

    start_time = time.time()

    # ---- create compilation of annual totals
    compilation          = pd.DataFrame() # stores individual years e.g., [1978 - 2022]
    compilation_averages = pd.DataFrame() # stores average over all years for a given iteration

    for i in range(n_draws):
        sample = draw_sample(SO2, Hg_ratios)
        # get annual totals for compilation
        annual_totals = sample.groupby(by=['Type','year'], as_index=False).sum(numeric_only=True)[['Type','year','Hg (Mg/a)']]
        annual_average = annual_totals.groupby(by=['Type'], as_index=False).mean(numeric_only=True)[['Type','Hg (Mg/a)']]

        annual_totals['iteration'] = i
        annual_average['iteration'] = i

        compilation = pd.concat((compilation, annual_totals))
        compilation_averages = pd.concat((compilation_averages, annual_average))
    
    print("--- %s seconds ---" % np.round((time.time() - start_time),2))
    compilation.to_csv(f'{output_path}bootstrap_emission_compilation_annual_totals_{tag_1990}.csv', index=False)
    compilation_averages.to_csv(f'{output_path}bootstrap_emission_compilation_annual_averages_{tag_1990}.csv', index=False)
    
    print('COMPLETED TASK: annual total emission estimates')
    return

# ------------------------------------------------------------------------------------------------------------------
# Call sequence -- Estimate Annual Total Hg emission magnitude for BOOTSTRAP scaling method
# ------------------------------------------------------------------------------------------------------------------
def call_bootstrap_volcano_list(n_draws: int, output_path:str, volcano_list:list, F_exclude_pre_1990_ratios=False):

    if F_exclude_pre_1990_ratios == True:
        tag_1990 = 'excl_pre_1990'
    else:
        tag_1990 = 'incl_pre_1990'

    # ---- load input data
    SO2, Hg_ratios, conv_keys = load_base_data(F_exclude_pre_1990_ratios)
    
    ## now go through each individual source and generate n_draws number of emission estimates
    for v in volcano_list:
        start_time = time.time()
        compilation = pd.DataFrame() # stores individual years e.g., [1978 - 2022]

        for i in range(n_draws):
            sample = draw_sample(SO2, Hg_ratios)
            sample = sample[sample['Volcano']==v] # subset sample draw to 
            sample['iteration'] = i
            compilation = pd.concat((compilation, sample))

        compilation_pvf = compilation[compilation['Type']=='pvf']
        compilation_eff = compilation[compilation['Type']=='eff']
        compilation_exp = compilation[compilation['Type']=='exp']

        compilation_pvf.to_csv(f'{output_path}All_Volcanoes/bootstrap_emission_compilation_{v}_pvf_{tag_1990}.csv', index=False)
        compilation_eff.to_csv(f'{output_path}All_Volcanoes/bootstrap_emission_compilation_{v}_eff_{tag_1990}.csv', index=False)
        compilation_exp.to_csv(f'{output_path}All_Volcanoes/bootstrap_emission_compilation_{v}_exp_{tag_1990}.csv', index=False)
        print(f"completed {v} --- %s seconds ---" % np.round((time.time() - start_time),2))
    print('COMPLETED TASK: single volcano emission estimates')
    return

# ------------------------------------------------------------------------------------------------------------------
# CALL SEQUENCE
# ------------------------------------------------------------------------------------------------------------------
SO2, Hg_ratios, conv_keys = load_base_data(F_exclude_pre_1990_ratios)

if flag_call_bootstrap_annual == True:
    call_bootstrap_annual(n_draws=n_iterations, output_path=output_path, F_exclude_pre_1990_ratios=F_exclude_pre_1990_ratios)
if flag_call_bootstrap_volcano_list == True:
    call_bootstrap_volcano_list(n_draws=n_iterations, output_path=output_path, volcano_list=volcano_list, F_exclude_pre_1990_ratios=F_exclude_pre_1990_ratios)

print('COMPLETED SCRIPT: 3-Estimate_Hg_Bootstrap.py')
