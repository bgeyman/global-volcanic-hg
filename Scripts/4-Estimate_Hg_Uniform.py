import numpy as np
import pandas as pd

print('STARTING SCRIPT: 4-Estimate_Hg_Uniform.py')

# ----------------------------------------------------------------
# Estimates Hg emissions based on SO2 emission magnitude and uniform Hg:SO2 ratio
# writes daily, (volcano, type)-level emission files to `../Output/Hg_emissions_row_format_uniform_only.csv`
# ----------------------------------------------------------------

uniform_HgSO2_ratio = 7.8e-6

# row-format merge of daily passive, effusive, explosive emissions
df = pd.read_csv('../Output/all_SO2.csv')
df['Datetime'] = pd.to_datetime(df['Datetime'])
df['Year'] = df.Datetime.dt.year
df = df.drop_duplicates(subset = df.columns, keep = 'first').reset_index(drop = True)
df = df.rename(columns={'1 sigma (Mg/day)': '1 sigma SO2 (Mg/day)'})

## get relative uncertainty for EXP and EFF emissions
df['uniform_relative_uncertainty'] = np.nan
c1 = (df['Type'].isin(['exp','eff']))
c2 = (df['VEI'] < 4.0)
df.loc[(c1 & c2), 'uniform_relative_uncertainty'] = 0.5
c3 = (df['VEI'] >= 4.0)
df.loc[(c1 & c3), 'uniform_relative_uncertainty'] = 0.2

df['a'] = (1-df['uniform_relative_uncertainty'])*df['SO2 (Mg/day)']
df['b'] = (1+df['uniform_relative_uncertainty'])*df['SO2 (Mg/day)']

df.loc[c1, '1 sigma SO2 (Mg/day)'] = np.sqrt( ((df['b']-df['a'])**2)/12 )

## now calcualte sigmas
sigma_HgSO2 = (1.5e-6)*np.sqrt(13)

# now calculate Hg emission
df['Hg uniform (Mg/day)'] = df['SO2 (Mg/day)']*uniform_HgSO2_ratio


df['1 sigma Hg (Mg/day)'] = (df['Hg uniform (Mg/day)']*
                             np.sqrt( (df['1 sigma SO2 (Mg/day)']/df['SO2 (Mg/day)'])**2 + 
                                      (sigma_HgSO2/uniform_HgSO2_ratio)**2 )) 
df = df[['Volcano', 'Country', 'Lat', 'Lon', 'Datetime', 'Type', 'V_alt',
        'P_alt', 'SO2 (Mg/day)', '1 sigma SO2 (Mg/day)', 'VEI', 'AMF',
        'f_imputed', 'f_alt_obs', 'Source', 'Year',
        'Hg uniform (Mg/day)', '1 sigma Hg (Mg/day)']]


df = df.rename(columns={'1 sigma Hg (Mg/day)': '1 sigma Hg uniform (Mg/day)'})

df = df[df['Year']<2022]

#### now export Hg_uniform emission estimates
df.to_csv('../Output/Data/Hg_emissions_row_format_uniform_only.csv', index=False)

print('COMPLETED SCRIPT: 4-Estimate_Hg_Uniform.py')