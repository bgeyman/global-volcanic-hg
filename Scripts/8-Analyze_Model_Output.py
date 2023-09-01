import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
from   Calculate_Atmospheric_Effects import *
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from   cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from   matplotlib.colors import ListedColormap, LinearSegmentedColormap
 
print('STARTING SCRIPT: 8-Analyze_Model_Output.py')

# ---------------------------------------------------------------
data_path = '../Output/Publication_Data/'
# ---------------------------------------------------------------

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

# -------------------------------
# Load data
ds = xr.open_mfdataset(data_path+'Model_Data.nc')
# now average over all years in record
ds_avg = ds.mean(dim='time')
# assign units
for v in ds_avg.data_vars:
    ds_avg[v].attrs = ds[v].attrs
# fully load into memory for faster calculations below
ds_avg = ds_avg.load()

# -------------------------------
# Calculate total atmospheric Hg content (Mg) by hemisphere
# -------------------------------
Hg_mass_total      = get_atmospheric_mass(ds_avg.sel(lat=slice(-90, 90)), conc_var='Total_Hg_Conc').sum().values.item()
Hg_mass_ETNH       = get_atmospheric_mass(ds_avg.sel(lat=slice(24,  90)), conc_var='Total_Hg_Conc').sum().values.item()
Hg_mass_tropics    = get_atmospheric_mass(ds_avg.sel(lat=slice(-24, 24)), conc_var='Total_Hg_Conc').sum().values.item()
Hg_mass_tropics_NH = get_atmospheric_mass(ds_avg.sel(lat=slice(0,   24)), conc_var='Total_Hg_Conc').sum().values.item()
Hg_mass_tropics_SH = get_atmospheric_mass(ds_avg.sel(lat=slice(-24, 0)),  conc_var='Total_Hg_Conc').sum().values.item()
Hg_mass_ETSH       = get_atmospheric_mass(ds_avg.sel(lat=slice(-90,-24)), conc_var='Total_Hg_Conc').sum().values.item()

print('----------------------------------')
print('Printing Model Output Diagnostics')
print('----------------------------------')

print('Total Hg mass in atmosphere: '.ljust(29), '{:.2f} Mg'.format(Hg_mass_total))
print('Hg mass in ETNH: '.ljust(29),             '{:.2f} Mg'.format(Hg_mass_ETNH))
print('Hg mass in tropics: '.ljust(29),          '{:.2f} Mg'.format(Hg_mass_tropics))
print('Hg mass in ETSH: '.ljust(29),             '{:.2f} Mg'.format(Hg_mass_ETSH))
print('--')
print('Hg mass is {:.2f}%'.format((Hg_mass_ETNH/Hg_mass_ETSH-1)*100) + ' greater in the extratropical NH than in SH')
print('Hg mass in tropics (NH): {:.2f} Mg'.format(Hg_mass_tropics_NH))
print('Hg mass in tropics (SH): {:.2f} Mg'.format(Hg_mass_tropics_SH))
print('----------------------------------')

# -------------------------------
# calculate surface Hg concentration enhancement for:
#   - extratropical NH
#   - extratropical SH
#   - tropics
# -------------------------------
def get_weighted_mean_surface_conc(ds, lat_bounds=[24, 90], conc_var='Total_Hg_Conc'):
    ds_surf = ds.sel(lat=slice(lat_bounds[0], lat_bounds[1])).isel(lev=-1)
    # weight concentrations by air density
    weighted_mean = (ds_surf[conc_var]*ds_surf['Met_AD']).sum()/ds_surf['Met_AD'].sum()
    return weighted_mean.values.item()

conc_surf_Global  = get_weighted_mean_surface_conc(ds_avg, lat_bounds=[-90, 90])*1e15
conc_surf_ETNH    = get_weighted_mean_surface_conc(ds_avg, lat_bounds=[24,  90])*1e15
conc_surf_tropics = get_weighted_mean_surface_conc(ds_avg, lat_bounds=[-24, 24])*1e15
conc_surf_ETSH    = get_weighted_mean_surface_conc(ds_avg, lat_bounds=[-90,-24])*1e15

rel_diff_NH_SH_pct = ((conc_surf_ETNH/conc_surf_ETSH)-1)*100
rel_diff_tropics_global_pct = ((conc_surf_tropics/conc_surf_Global)-1)*100
print('Surface Hg concentration Global: '.ljust(37),      '{:.2f} ppq'.format(conc_surf_Global))
print('Surface Hg concentration in ETNH: '.ljust(37),     '{:.2f} ppq'.format(conc_surf_ETNH))
print('Surface Hg concentration in tropics: '.ljust(37),  '{:.2f} ppq'.format(conc_surf_tropics))
print('Surface Hg concentration in ETSH: '.ljust(37),     '{:.2f} ppq'.format(conc_surf_ETSH))
print('--')
print('Extratropical NH surface Hg conc is {:.2f}%'.format(rel_diff_NH_SH_pct) + ' greater than in SH')
print('--')
print('Surface Hg conc is {:.2f}%'.format(rel_diff_tropics_global_pct) + ' greater in tropics than the global average')
print('----------------------------------')

# -------------------------------
# Calculate total atmospheric Hg deposition by hemisphere (and within tropics)
# -------------------------------
dep_global     = get_global_deposition_Mg(ds_avg)['total_dep']
dep_ETNH       = get_global_deposition_Mg(ds_avg.sel(lat=slice(24, 90)))['total_dep']
dep_tropics    = get_global_deposition_Mg(ds_avg.sel(lat=slice(-24, 24)))['total_dep']
dep_tropics_NH = get_global_deposition_Mg(ds_avg.sel(lat=slice(0, 24)))['total_dep']
dep_tropics_SH = get_global_deposition_Mg(ds_avg.sel(lat=slice(-24, 0)))['total_dep']
dep_ETSH       = get_global_deposition_Mg(ds_avg.sel(lat=slice(-90, -24)))['total_dep']

area_global  = ds_avg['AREA'].sum()
area_ETNH    = ds_avg.sel(lat=slice(24,90))['AREA'].sum()
area_tropics = ds_avg.sel(lat=slice(-24,24))['AREA'].sum()
area_ETSH    = ds_avg.sel(lat=slice(-90,-24))['AREA'].sum()

print('Total Hg deposition in ETNH: '.ljust(32),    '{:.2f} Mg'.format(dep_ETNH))
print('Total Hg deposition in tropics: '.ljust(32), '{:.2f} Mg'.format(dep_tropics))
print('Total Hg deposition in ETSH: '.ljust(32),    '{:.2f} Mg'.format(dep_ETSH))
print('--')
print('Avg. Hg deposition, Global (area weighted): '.ljust(47),    '{:.2f} ug m-2 a-1'.format((dep_global*1e12)/area_global))
print('Avg. Hg deposition in ETNH (area weighted): '.ljust(47),    '{:.2f} ug m-2 a-1'.format((dep_ETNH*1e12)/area_ETNH))
print('Avg. Hg deposition in tropics (area weighted): '.ljust(47), '{:.2f} ug m-2 a-1'.format((dep_tropics*1e12)/area_tropics))
print('Avg. Hg deposition in ETSH (area weighted): '.ljust(47),    '{:.2f} ug m-2 a-1'.format((dep_ETSH*1e12)/area_ETSH))

dep_rel_diff_tropics_global_pct = ( ((dep_tropics*1e12)/area_tropics)/((dep_global*1e12)/area_global) - 1)*100
print('Avg. Hg deposition is {:.2f}%'.format(dep_rel_diff_tropics_global_pct) + ' greater in tropics than the global average')

print('----------------------------------')

# ---------------------------------------------------------------------------------------------------
# Plot concentration and deposition figure
# ---------------------------------------------------------------------------------------------------

# load custom colormap from color_lists.json
with open('../Data/color_lists.json') as f:
    color_lists = json.load(f)
    f.close()

# set levels for concentration color bar
conc_levels = np.linspace(0,12,25)
conc_levels = np.append(conc_levels, 20)
conc_cmap = ListedColormap(color_lists['concentration'])
conc_cmap.set_over('red')

# set levels for deposition color bar
dep_levels = np.linspace(0,1.6,33)
dep_levels = np.append(dep_levels, 2.)
dep_cmap = ListedColormap(color_lists['deposition'])
dep_cmap.set_over('red')

# set scale factor to convert from mol/mol to ppq
scale_factor = 1e15

# set ticks for color bars
ticks_conc = [0, 3, 6, 9, 12]
ticks_dep  = [0,0.4,0.8,1.2,1.6]

central_longitude = 180
projection = ccrs.Robinson(central_longitude=central_longitude)

plt.figure(figsize=(12,5))
# -- plot concentration (Hg0)
ax = plt.subplot(121, projection=projection)
ax.coastlines(lw=0.5)
ds_volc = ds_avg['SpeciesConc_Hg0'].isel(lev=slice(-3,None)).mean(dim=['lev'])*scale_factor
ds = (ds_volc)
ds.plot(x='lon', y='lat', ax=ax, 
        transform=ccrs.PlateCarree(), 
        cmap=conc_cmap,
        levels=conc_levels,
        cbar_kwargs={'shrink':0.8, 'orientation':'horizontal', 'label':'$\Delta$ GEM [ppq]','ticks':ticks_conc})
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle=':')
ax.set_title('Concentration')

# plot deposition
ax = plt.subplot(122, projection=projection)
ax.coastlines(lw=0.5)
ds_volc = ds_avg['Total_Hg_Dep']
ds = (ds_volc)
ds.plot(x='lon', y='lat', ax=ax, 
        transform=ccrs.PlateCarree(), 
        cmap=dep_cmap, #new_cmaps['lajolla'], #['lajolla']
        levels=dep_levels,
        cbar_kwargs={'shrink':0.8, 'orientation':'horizontal','label':'$\Delta$ Deposition [$\mu g \ m^{-2} \ a^{-1}$]','ticks':ticks_dep})
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle=':')
ax.set_title('Deposition')
# save figure
plt.savefig(f'../Output/Figures/Fig1_Concentration_Deposition_Map.png', dpi=500)

print('COMPLETED SCRIPT: 8-Analyze_Model_Output.py')