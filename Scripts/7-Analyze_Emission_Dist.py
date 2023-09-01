import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

print('STARTING SCRIPT: 7-Analyze_Emission_Dist.py')

# ---------------------------------------------------------------
data_path = '../Output/Data/'
# ---------------------------------------------------------------

df = pd.read_csv(data_path+'Hg_emissions_row_format_uniform_only.csv')
df = df.groupby(by=['Year','Type','Volcano'], as_index=False).agg({'Lat':'first','SO2 (Mg/day)':'sum'})
df = df.rename(columns={'SO2 (Mg/day)':'SO2 (Mg/a)'})

lat_bounds = {'Antarctic':   [  -90, -66.5],
              'MLSH':        [-66.5, -23.5],
              'Tropical SH': [-23.5,     0],
              'Tropical NH': [    0,  23.5],
              'MLNH':        [ 23.5,  66.5],
              'Arctic':      [ 66.5,    90],
              }

regions = list(lat_bounds.keys())

lat_emissions = {'region':regions,
                 'pvf':[],
                 'pvf_lo':[],
                 'pvf_hi':[],
                 'eff':[],
                 'eff_lo':[],
                 'eff_hi':[],
                 'exp':[],
                 'exp_lo':[],
                 'exp_hi':[],
                 }

estimate_summary = pd.read_csv('../Output/Tables/Table1_Emission_Estimate_Summary.csv')

r_HgSO2 = {} # -- get Hg:SO2 ratio from emission totals table
n_years = {} # -- get n years from emission totals table
for type_sel in ['pvf','eff','exp']:
    tmp = estimate_summary[estimate_summary['Type']==type_sel]
    n_years[type_sel] = (tmp['Year End'] - tmp['Year Start']).values.item() + 1
    r_HgSO2[type_sel] = {}
    for val in ['Central', 'Lower', 'Upper']:
        r_HgSO2[type_sel][val] = (tmp[f'Bootstrap {val} (Mg/a)'].values/tmp['SO2 mean (Tg/a)'].values).item()*1e-6

for type_sel in ['pvf','eff','exp']:
    c1  = (df['Type']==type_sel)
    for region in regions:
        c2 = (df['Lat']> lat_bounds[region][0])
        c3 = (df['Lat']<=lat_bounds[region][1])
        tmp = df[c1&c2&c3]
        E   = tmp['SO2 (Mg/a)'].sum()/n_years[type_sel]
        lat_emissions[f'{type_sel}'].append(E*r_HgSO2[type_sel]['Central'])
        lat_emissions[f'{type_sel}_lo'].append(E*r_HgSO2[type_sel]['Lower'])
        lat_emissions[f'{type_sel}_hi'].append(E*r_HgSO2[type_sel]['Upper'])

# ----------- now get areas for each region
def do_grid (resolution=1):
    """Calculate the area of each grid cell for a user-provided
    grid cell resolution. Area is in square meters, but resolution
    is given in decimal degrees."""
    # Calculations needs to be in radians
    lats = np.deg2rad(np.arange(-90+0.5*resolution, 90+resolution, resolution))
    r_sq = 6371000**2
    n_lats = int(360./resolution) 
    area = r_sq*np.ones(n_lats)[:, None]*np.deg2rad(resolution)*(
                np.sin(lats[1:]) - np.sin(lats[:-1]))
    return area.T

res = 0.1
area = do_grid(resolution=res)
ds = xr.Dataset({'area': xr.DataArray(
                data   = area,
                dims   = ['lat','lon'],
                coords = {'lat': np.arange(-90+0.5*res, 90, res),
                          'lon': np.arange(-180, 180, res)},
                attrs  = {'units' : 'm2'}),
            })

# confirm area is within 0.1% of reference value
assert np.abs((ds['area'].sum()/5.10072e14)-1) < 1e-3

region_areas = {}
for region in regions:
    lat_lo = lat_bounds[region][0]
    lat_hi = lat_bounds[region][1]
    region_areas[region] = {}
    region_areas[region]['m2'] = ds.sel(lat=slice(lat_lo, lat_hi))['area'].sum().values.item()
    region_areas[region]['f_total'] = region_areas[region]['m2']/ds['area'].sum().values.item()

print('-----')
print('Printing diagnostic output for annotation of Fig. S1')
print('-----')
#print('region areas')
#print(region_areas)

# ----------------------------------------------------
# make plot

# -- set plotting parameters
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = False
mpl.rcParams['xtick.top'] = False
mpl.rcParams['xtick.bottom'] = False
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5

colors = [(21/255, 94/255, 171/255),
          (209/255, 72/255, 35/255),
          (233/255, 160/255, 29/255),]

# create the xticks beginning a index 0
x_labels = lat_emissions['region']
xticks = range(len(x_labels))

# plot
fig, ax = plt.subplots(figsize=(8, 4))
bottom = np.zeros(len(x_labels))
for (type_sel, c) in zip(['pvf','eff','exp'], colors):
    ax.bar(x=xticks, height=lat_emissions[type_sel], bottom=bottom, tick_label=x_labels, label=type_sel, facecolor=c, edgecolor='k', lw=1, zorder=2)
    bottom += lat_emissions[type_sel]

# --- Print out the values for the figure
for i in range(len(lat_emissions['region'])):
    absolute = bottom[i]
    relative = bottom[i]/np.sum(bottom)
    print(f'{x_labels[i]}'.ljust(19), f'| Absolute (Mg) : {np.round(absolute, 1)}'.ljust(23), f'| Relative : {np.round(relative, 2)}')

print('------')
SH_absolute = np.sum(bottom[:3])
SH_relative = SH_absolute/np.sum(bottom)
NH_absolute = np.sum(bottom[3:])
NH_relative = NH_absolute/np.sum(bottom)
print(f'Southern Hemisphere | Absolute (Mg) : {np.round(SH_absolute, 1)} | Relative : {np.round(SH_relative, 2)}')
print(f'Northern Hemisphere | Absolute (Mg) : {np.round(NH_absolute, 1)} | Relative : {np.round(NH_relative, 2)}')

plt.axhline(y=130, xmin=0.02, xmax=0.49, c='dimgrey', lw=3)
plt.axhline(y=130, xmin=0.51, xmax=0.98, c='dimgrey', lw=3)

plt.ylim(0,140)
ax.minorticks_on()
ax.yaxis.set_tick_params(which='major',right=True, labelright=True)

plt.savefig('../Output/Figures/Emissions_by_Latitude_Mg.png', format='png', dpi=500)

print('COMPLETED SCRIPT: 7-Analyze_Emission_Dist.py')