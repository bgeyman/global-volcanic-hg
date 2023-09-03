import numpy as np
import pandas as pd

print('STARTING SCRIPT: 5-Make_Hg_Emission_Table.py')

# -------------------------------------------------------------------
# Establish which quantiles to use for bootstrap estimates
central_estimate_quantile = 0.50
lower_bound_quantile      = 0.25
upper_bound_quantile      = 0.75
bootstrap_columns = ['Type','Bootstrap Central (Mg/a)', 'Bootstrap Lower (Mg/a)', 'Bootstrap Upper (Mg/a)', 'Quantiles Used']

# -------------------------------------------------------------------
# TABLE 1
# make summary table of average Hg emission and uncertainty (average over all years in record)
# ------------------------------------------------------------------

# -- Uniform -- 
df = pd.read_csv('../Output/Data/Hg_emissions_row_format_uniform_only.csv')

event_count    = df.groupby(by=['Type','Year','Volcano'], as_index=False).sum(numeric_only=True).groupby(by=['Type','Year'], as_index=False).agg({'Volcano':'count'})
annual_sum     = df.groupby(by=['Type','Year'], as_index=False).sum(numeric_only=True)
annual_average = annual_sum.groupby(by=['Type'], as_index=False).mean(numeric_only=True)

output_table = pd.DataFrame({'Type': annual_average['Type'],
                             'Year Start': 0, # filler values -- replaced in loop below
                             'Year End':   0,
                             'N': event_count['Volcano'],
                             'SO2 mean (Tg/a)': (annual_average['SO2 (Mg/day)']*1e-6),
                             'SO2 std (Tg/a)':  (annual_average['1 sigma SO2 (Mg/day)']*1e-6),
                             'Uniform mean (Mg/a)': annual_average['Hg uniform (Mg/day)'],
                             'Uniform std (Mg/a)':  annual_average['1 sigma Hg uniform (Mg/day)'],
                             })

for type_sel in ['pvf', 'eff', 'exp']:
    len_tmp  = len(df[df['Type']==type_sel]['Volcano'].unique())
    yr_start = int(df[df['Type']==type_sel]['Year'].min())
    yr_end   = int(df[df['Type']==type_sel]['Year'].max())
    output_table.loc[output_table['Type']==type_sel, 'N'] = len_tmp
    output_table.loc[output_table['Type']==type_sel, 'Year Start'] = yr_start
    output_table.loc[output_table['Type']==type_sel, 'Year End'] = yr_end

# -- Bootstrap --
bootstrap_average = pd.read_csv('../Output/Data/Bootstrap_Output/bootstrap_emission_compilation_annual_averages_excl_pre_1990.csv')

output = pd.DataFrame(columns=bootstrap_columns)
for type_sel in ['pvf', 'eff', 'exp']:
    c1 = (bootstrap_average['Type']==type_sel)

    central_estimate = np.quantile(bootstrap_average[c1]['Hg (Mg/a)'], central_estimate_quantile)
    lower_bound      = np.quantile(bootstrap_average[c1]['Hg (Mg/a)'], lower_bound_quantile)
    upper_bound      = np.quantile(bootstrap_average[c1]['Hg (Mg/a)'], upper_bound_quantile)

    tmp = pd.DataFrame({'Type':[type_sel],
                        'Bootstrap Central (Mg/a)': [central_estimate],
                        'Bootstrap Lower (Mg/a)':[lower_bound], 
                        'Bootstrap Upper (Mg/a)':[upper_bound], 
                        'Quantiles Used':[f'{central_estimate_quantile} ({lower_bound_quantile} - {upper_bound_quantile})']})
    
    output = pd.concat((output, tmp))
    output[output.columns[1:-1]] = output[output.columns[1:-1]].astype(float)

output_table = pd.merge(left=output_table, right=output[bootstrap_columns], on='Type')
output_table = output_table[output_table['Type'].isin(['pvf','eff','exp'])] # drop additional rows

for c in output_table.columns[4:-1]:
    output_table[c] = np.round(output_table[c],3) # Increase precision for now. Round at end.

output_table.to_csv('../Output/Tables/Table1_Emission_Estimate_Summary.csv', index=False)

print('COMPLETED SCRIPT: 5-Make_Hg_Emission_Table.py')