
import numpy as np
import pandas as pd
import pylab as pl
import matplotlib.pyplot as plt
import json
import sys

sys.path.append('./flexbox/')
import flexbox
from sevenbox import init_model
from flexbox_helpers import import_metastring, print_steady

print('STARTING SCRIPT: 10-sevenbox_volcano_Hg.py')

# load flexbox model inputs from `flexbox_model_inputs_volcano_Hg.json`
# this dictionary contains both the model values required as 
# inputs to sevenbox and also some metadata such as description
# and reference 
with open('./flexbox/inputs/flexbox_model_inputs_volcano_Hg.json', 'r') as fn:
    inputs_meta = json.load(fn)

# strip the metadata from inputs_meta to create inputs dictionary
# to be passed as argument to function `init_model()` from sevenbox
inputs = {}
for key in inputs_meta.keys():
    inputs[key] = inputs_meta[key]['value']

model = init_model(input_dict=inputs)
model.build()

# weights dict not used for natural budget
weights_dict = {'f_tf':0.4, 'f_ts':0.3, 'f_ta':0.3, 'w_tfm':0., 'w_tsm':0.05,'w_tam':0.95, 'f_sequestered':0.5}
meta = import_metastring(weights_dict=weights_dict)

# ------------------------------------------------------------------------------
# Calculate standard pre-industrial budget
#--------------------------------------------------------------------------------
compartments = {'atm':0,'tf':1,'ts':2,'ta':3,'ocs':4,'oci':5,'ocd':6}

steady = print_steady(model, volcanic_air=230, volcanic_ocean=50, print_output=False)

# -----------------------------------
# Write reservoir magnitudes to file
# --
df_res = pd.DataFrame()  
for i in range(len(model.compartments)):
    tmp_out = pd.DataFrame({'Reservoir':[model.compartments[i]],
                            'Mass (Mg)': [int(np.round(steady[i],0))]})
    df_res = pd.concat((df_res, tmp_out))

df_res['Reservoir']   = pd.Categorical(df_res['Reservoir'],   ['atm','ocs','oci','ocd','tf','ts','ta', 'none'])

df_res = df_res.sort_values(by=['Reservoir'])

with open('../Output/Publication_Data/pre_industrial_reservoirs.txt', 'w') as f:
    f.write(df_res.to_string(header=False, index=False))

print('------------------------------------------------')
print('Pre-industrial reservoirs (Mg)')
print('-----')
print(df_res.to_string(index=False))

# -----------------------------------
# Write fluxes to file
# --
df_flux = pd.DataFrame()
for rate in model.rates:
    source_index = compartments[rate.compartment_from]
    source_res = steady[source_index]
    tmp = rate.k*source_res
    tmp_out = pd.DataFrame({'compartment_from':[rate.compartment_from], 
                            'compartment_to':[rate.compartment_to],
                            'Flux (Mg/a)':[str(np.round(tmp,1))], 
                            'notes':[rate.notes]})
    
    df_flux = pd.concat((df_flux, tmp_out))
df_flux = df_flux.reset_index(drop=True)
df_flux['compartment_from'] = pd.Categorical(df_flux['compartment_from'], ['atm','ocs','oci','ocd','tf','ts','ta', 'none'])
df_flux['compartment_to']   = pd.Categorical(df_flux['compartment_to'],   ['atm','ocs','oci','ocd','tf','ts','ta', 'none'])

df_flux = df_flux.sort_values(by=['compartment_from','compartment_to'])
# read labels mapped from notes
labels = pd.read_csv('../Data/budget_table_labels.csv')[['notes','label']]

df_flux = pd.merge(df_flux, labels, on='notes').fillna('none')
df_flux = df_flux.rename(columns={'label':'Description'})

df_flux['Coupling']    = df_flux['compartment_from'].astype(str).str.pad(4, side='left') + ' --> ' + df_flux['compartment_to'].astype(str).str.pad(4, side='right')
df_flux['Description'] = df_flux['Description'].astype(str).str.pad(59, side='right')
df_flux['Flux (Mg/a)'] = df_flux['Flux (Mg/a)'].astype(str).str.pad(6, side='right')
df_flux = df_flux[['Coupling','Description','Flux (Mg/a)']]

print('------------------------------------------------')
print('Pre-industrial fluxes (Mg/a)')
print('-----')
print(df_flux.to_string(index=False))
with open('../Output/Publication_Data/pre_industrial_fluxes.txt', 'w') as f:
    f.write(df_flux.to_string(header=False, index=False))


# ---------------------------------------------------------------------------------------------------- #
# Sensitivity test
# -
# - store select budget output under range of hydrothermal and subaerial flux estimates in DataFrame-- #
output = {'hydrothermal flux [Mg a-1]':[], 'subaerial volcanism flux [Mg a-1]':[], 'atmosphere [Mg]':[], 'surface ocean [Mg]':[], 'intermediate ocean [Mg]':[]}

for hydrothermal_flux in [20, 50, 80]:
    for volcano_air in [160, 230, 330]:
        steady = print_steady(model, volcanic_air=volcano_air, volcanic_ocean=hydrothermal_flux, print_output=False)
        output['hydrothermal flux [Mg a-1]'].append(hydrothermal_flux)
        output['subaerial volcanism flux [Mg a-1]'].append(volcano_air)
        output['atmosphere [Mg]'].append(np.round(steady[0],1))
        output['surface ocean [Mg]'].append(np.round(steady[4],1))
        output['intermediate ocean [Mg]'].append(np.round(steady[5],1))
output = pd.DataFrame(output)

def subset_df(df, kv_dict:dict):
    for key in kv_dict.keys():
        df = df[df[key]==kv_dict[key]]
    return df

# effect of hydrothermal inputs on surface ocean Hg
subaerial_volcanic_flux = 230    # Mg a-1
hydrothermal_flux = [80, 50, 20] # [high, mid, low] Mg a-1

hydrothermal_hi  = subset_df(df=output, kv_dict={'subaerial volcanism flux [Mg a-1]':subaerial_volcanic_flux, 'hydrothermal flux [Mg a-1]':hydrothermal_flux[0]})
hydrothermal_lo  = subset_df(df=output, kv_dict={'subaerial volcanism flux [Mg a-1]':subaerial_volcanic_flux, 'hydrothermal flux [Mg a-1]':hydrothermal_flux[2]})
hydrothermal_mid = subset_df(df=output, kv_dict={'subaerial volcanism flux [Mg a-1]':subaerial_volcanic_flux, 'hydrothermal flux [Mg a-1]':hydrothermal_flux[1]})
diff_surf     = (hydrothermal_hi['surface ocean [Mg]'].item() - hydrothermal_lo['surface ocean [Mg]'].item())/2
rel_diff_surf = diff_surf/hydrothermal_mid['surface ocean [Mg]'].item()
diff_int     = (hydrothermal_hi['intermediate ocean [Mg]'].item() - hydrothermal_lo['intermediate ocean [Mg]'].item())/2
rel_diff_int = diff_int/hydrothermal_mid['intermediate ocean [Mg]'].item()

# print comparisons described in manuscript
print('------------------------------------------------')
print('Sensitivity of surface ocean Hg to hydrothermal Hg flux:')
print('-----')
print('Difference in ocean Hg mass between high and low hydrothermal fluxes (80 vs. 20 Mg a-1): ')
print('-- Relative differences [%] --')
print('   surface ocean', np.round(rel_diff_surf,3)*1e2, '%')
print('   intermediate ocean', np.round(rel_diff_int,3)*1e2, '%')
print('-- Absolute differences [Mg] --')
print('   surface ocean', np.round(diff_surf,1), 'Mg')
print('   intermediate ocean', np.round(diff_int,1), 'Mg')
print('------------------------------------------------')
print('Atmospheric Hg mass under range of subaerial volcanism estimates:')
print('-----')
for volcano_air in [160, 230, 330]:
    steady = print_steady(model, volcanic_air=volcano_air, volcanic_ocean=50, print_output=False)
    print(f'subaerial volcanism {volcano_air} Mg a-1 --> atmospheric masss {np.round(steady[0],1)} Mg')
print('------------------------------------------------')

print('COMPLETED SCRIPT: 10-sevenbox_volcano_Hg.py')