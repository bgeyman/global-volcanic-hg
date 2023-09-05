import pandas as pd
import numpy as np

import os
import os.path

# removes intermediate files
os.remove('../Output/Data/Hg_emissions_row_format_uniform_only.csv')

# -- adds `total` to category emission estimates
df = pd.read_csv('../Output/Tables/Table1_Emission_Estimate_Summary.csv')

total = {}
for column in df.columns:
    if column in ['Type']:
        total[column] = 'All Categories'
    if column in ['Year Start','Year End','N']:
        total[column] = ''
    if column in ['SO2 mean (Tg/a)','SO2 std (Tg/a)', 'Uniform mean (Mg/a)', 'Uniform std (Mg/a)',
                  'Bootstrap Central (Mg/a)', 'Bootstrap Lower (Mg/a)', 'Bootstrap Upper (Mg/a)']:
        total[column] = df[column].sum()
    if column in ['Quantiles Used']:
        total[column] = df[column].iloc[0]

df = pd.concat((df, pd.DataFrame(total, index=[0])))

# -- round to 1 decimal place for publication
for c in df.columns[4:-1]:
    df[c] = np.round(df[c],1)

df.to_csv('../Output/Tables/Table1_Emission_Estimate_Summary.csv', index=False)
# -- 