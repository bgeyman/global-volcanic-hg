import numpy as np
import pandas as pd
from helpers import subset_df

print('STARTING SCRIPT: 6-Combine_RF_Hg_Emissions.py')

# ---------------------------------------------------------------------------
# Adds bootstrap Hg emission estimate to daily Hg emission data by scaling
# uniform estimate for each source according to a scale factor based on the ratio
# between the long-term (summary) bootstrap Hg emission estimate (median) 
# for that category and the uniform estimate
# ---------------------------------------------------------------------------

df = pd.read_csv('../Output/Data/Hg_emissions_row_format_uniform_only.csv')
summary_table = pd.read_csv('../Output/Tables/Table1_Emission_Estimate_Summary.csv')

df['Hg bootstrap (Mg/day)'] = np.nan

for type_sel in ['pvf','eff','exp']:
    tmp = subset_df(df=summary_table, kv_dict={'Type':type_sel})
    uniform_estimate   = tmp['Uniform mean (Mg/a)'].values.item()
    bootstrap_estimate = tmp['Bootstrap Central (Mg/a)'].values.item()
    scale_factor = bootstrap_estimate/uniform_estimate

    df['Hg bootstrap (Mg/day)'] = np.where(df['Type']==type_sel, 
                                           df['Hg uniform (Mg/day)']*scale_factor, 
                                           df['Hg bootstrap (Mg/day)'])
    
df.to_csv('../Output/Data/Hg_emissions_row_format_all.csv', index=False)

print('COMPLETED SCRIPT: 6-Combine_RF_Hg_Emissions.py')