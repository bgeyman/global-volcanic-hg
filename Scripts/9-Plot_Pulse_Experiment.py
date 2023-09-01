import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
import json

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

print('STARTING SCRIPT: 9-Plot_Pulse_Experiment.py')


# ---------------------------------------------------------------
data_path = '../Output/Publication_Data/'
# ---------------------------------------------------------------

# read colors from color_lists.json
with open('../Data/color_lists.json') as json_file:
    color_lists = json.load(json_file)

# read categorical colors for pulse levels
colors = color_lists['pulse_levels_categorical']

# -------------------------------------
df = pd.read_csv(data_path+'pulsed_emission_compilation.csv')
df = df[df['time']>=1]
df = df[df['level']!='2D']
df = df[df['level']!='D']
df['level'] = df['level'].astype('int')
df =df.sort_values(by=['location','level','species','month','time'])
# -------------------------------------

relative_y_axis = False # flag for whether y-axis is represented as a fraction of initial pulse mass
# initialize figure
plt.figure(figsize=(8,6))
# outer loop over species
for (plt_index, species) in zip([1,2,3,4], ['Hg0', 'HGCL2', 'HGBR2', 'HG2CLP']):
    plt_number = int(f'22{plt_index}')
    ax = plt.subplot(plt_number)
    ax.set_title(species)
    # subset dataframe to {species}
    df_species = df[df['species']==species]
    level_output = {}
    # inner loop over levels
    for (lev, c) in zip([5, 10, 20, 30, 35], colors):
        # subset dataframe to {species, level}
        df_species_level = df_species[df_species['level']==lev]
        # calculate time-series of Hg mass remaining
        tmp = df_species_level.groupby(by=['time'], as_index=False).agg({'Hg [Mg]':['mean', 'std', 'max', 'min']})
        
        # illustrate initial mass 
        ax.scatter(0, 100, facecolor='k', s=20, edgecolor='grey', zorder=3,)
        ax.plot([0, tmp['time'][0]], [100, tmp['Hg [Mg]']['mean'][0]], 
                markersize=5,
                c='lightgrey', ls=':', zorder=2)
        
        # plot time-series
        ax.plot(tmp['time'], tmp['Hg [Mg]']['mean'], c=c, lw=2, 
                markersize=5,
                markerfacecolor=c, markeredgecolor='k', markeredgewidth=0.5, label=lev, marker='o')
        # -- fill between (mean, std) range
        ax.fill_between(x=tmp['time'], 
                        y1=tmp['Hg [Mg]']['mean']-tmp['Hg [Mg]']['std'], 
                        y2=tmp['Hg [Mg]']['mean']+tmp['Hg [Mg]']['std'], 
                        facecolor=c, alpha=0.1)
        
    ax.set_xticks(np.arange(0,85,10))
    xlabels = list(np.arange(0,85,10).astype(int))
    for i in range(1, len(xlabels), 2):
        xlabels[i] = ""
    ax.set_xticklabels(labels=xlabels, fontsize=10)
    ax.set_xlim(-5, 90)
    ax.set_xlabel('Time After Eruption (days)', fontsize=10)
    ## rescale y-axis to fraction of mass remaining
    if relative_y_axis == True:
        ax.set_ylabel('Fraction of Hg Mass Remaining', fontsize=10)
    else:
        ax.set_ylim(0,105)
        ax.set_yticks(np.linspace(0,100,11))
        ylabels = list(np.linspace(0,100,11).astype(int))
        for i in range(1, len(ylabels), 2):
            ylabels[i] = ""
        ax.set_yticklabels(labels=ylabels, fontsize=10)
    
    ax.yaxis.set_tick_params(which='both', right='true')
    ax.xaxis.set_tick_params(which='both', top='true')
    ax.grid(which='major', axis='y', color='#DDDDDD', lw=0.2)
    ax.grid(which='minor', axis='y', color='#EEEEEE', ls=':', lw=0.1)
    ax.minorticks_on()
    if plt_index==4:
        ax.legend(fontsize=14, ncol=1, loc=(1.1,0))
plt.tight_layout()
plt.savefig('../Output/Figures/Fig2b_Pulse_Test.png', format='png', dpi=600)

levels = [5, 10, 20, 30, 35, 45]

def calculate_e_folding_time(f_mass, time):
    '''
    Inputs: `f_mass`, an array of length n with the fraction of the original mass remaining
    Example:
    f_mass = [1, 0.8, 0.7, 0.65, 0.62]
    calculate_e_folding_time(f_mass=f_mass)
    '''
    # take natural log of fraction of mass remaining
    y = np.log(list(f_mass))
    # plot number of elapsed months
    end_month = len(f_mass)-0.5
    x = time
    lm = stats.linregress(x,y)
    print_diags = False
    if print_diags==True:
        print(lm.rvalue**2)
        print(lm.intercept)
        print(lm.slope)

        plt.scatter(x, np.exp(y))
        y_mod = np.exp(lm.slope*x + lm.intercept)
        plt.plot(x, y_mod)
    
    halflife_days = np.log(2)/(-1*lm.slope) # halflife in months    
    return halflife_days

results = pd.DataFrame()

for location in ['MSH','Pinatubo']:
    for level in levels:
        for month in np.arange(1,4):
            for species in ['Hg0','HGCL2','HGBR2','HG2CLP']:

                sel = df.copy()
                time = 9.9
                c1 = (sel['level']==level)
                c2 = (sel['species']==species)
                c3 = (sel['location']==location)
                c4 = (sel['month']==month)
                c5 = (sel['time']>time)

                sel = sel[(c1&c2&c3&c4&c5)]
                ## intial decline - 10 days
                initial_decline = (100. - sel[sel['time']==10]['Hg [Mg]'])
                initial_decline = initial_decline.values.item()

                if sel['Hg [Mg]'].iloc[0] > 0:
                    f_mass = sel['Hg [Mg]']/100.
                    time = sel['time']
                    e_folding_time = calculate_e_folding_time(f_mass, time)
                else:
                    e_folding_time = np.nan

                tmp_output = pd.DataFrame({'location':[location], 'level':[level], 
                                           'month':[month], 'species':[species], 
                                           'initial decline':[initial_decline],
                                           'lifetime [days]':[e_folding_time]})
                results = pd.concat((results, tmp_output))

full_results = results.copy()
results = results.groupby(by=['level','species'], as_index=False).agg({'initial decline':['mean','std'], 'lifetime [days]':['mean','std']})

levels = [5, 10, 20, 30, 35]
data = {}
for species in ['Hg0','HGCL2', 'HGBR2', 'HG2CLP']:
    data[species] = {'initial_decline':[],
                     'initial_decline_std':[],
                     'fractional_initial_decline':[],
                     'fractional_initial_decline_std':[],
                     'lifetime':[],
                     'lifetime_std':[]}
    for lev in levels:
        sel = results[((results['species']==species)&(results['level']==lev))]
        data[species]['initial_decline'].append(sel['initial decline']['mean'].values.item())
        data[species]['initial_decline_std'].append(sel['initial decline']['std'].values.item())
        
        data[species]['fractional_initial_decline'].append(1-(sel['initial decline']['mean'].values.item())/100)
        data[species]['fractional_initial_decline_std'].append((sel['initial decline']['std'].values.item())/100)
        
        data[species]['lifetime'].append(sel['lifetime [days]']['mean'].values.item())
        data[species]['lifetime_std'].append(sel['lifetime [days]']['std'].values.item())
        
        initial_decline = sel['initial decline']['mean'].values.item()
        initial_decline_std = sel['initial decline']['std'].values.item()        

# --- now plot vertical version
markersize = 12
lw=3
X = np.arange(5)
fig = plt.figure(figsize=(6,6))
ax = plt.subplot(111)
x_offset = 0
for (species, color, label) in zip(['Hg0','HGCL2','HGBR2','HG2CLP'],
                                   [(0.8,0.8,0.8), (0.6,0.6,0.6), (0.4,0.4,0.4), (0.2,0.2,0.2)],
                                   ['Hg$^0$(g)','HgCl${_2}$(g)','HgBr$_2$(g)','HgCl${_2}$(p)']):
    ax.errorbar(data[species]['fractional_initial_decline'], 
                X+x_offset,  
                xerr=data[species]['fractional_initial_decline_std'], 
                color=color, markeredgecolor='k', marker='o', markersize=markersize, 
                linewidth=lw, ls='none', zorder=4, label=label)
    x_offset+=0.2

for i in np.arange(5):
    ax.axhspan(ymin=X[i]-0.2, ymax=X[i]+0.8, facecolor=colors[i], edgecolor='None', alpha=0.1)
    ax.axhspan(ymin=X[i]-0.2, ymax=X[i]+0.8, xmin=0, xmax=0.03, facecolor=colors[i], edgecolor='None', alpha=1, zorder=2)

ax.xaxis.set_tick_params(which='both', right='true')
ax.yaxis.set_tick_params(which='both', length=0)
ax.set_yticks(ticks=X+0.3, labels=levels)
ax.set_xticks(ticks=np.arange(0,1.2,0.2))
ax.set_ylim(-0.2, 4.8)
ax.set_ylabel('Grid Level')
ax.set_xlabel('fractional mass remaining (t=10 days)')
ax.grid(which='major', axis='x', color='#DDDDDD', lw=0.5)
ax.grid(which='minor', axis='x', color='#EEEEEE', ls=':', lw=0.4)
ax.minorticks_on()
ax.legend(loc=(0,1.05), ncols=4)

plt.savefig('../Output/Figures/Fig2a_Pulse_Test.png', format='png', dpi=600)
# -----------------------------------------------------
# now print lifetime post-10 days
cat = 'lifetime [days]'
levels = [5,10,20,30]
all_species = ['HGCL2','HGBR2','HG2CLP','Hg0']
c1 =(full_results['species'].isin(all_species))
c2 =(full_results['level'].isin(levels))
mu = np.nanmean(full_results[(c1&c2)][cat])
std = np.nanstd(full_results[(c1&c2)][cat])
print(f'{cat} of {all_species} emitted to L{levels}: {np.round(mu, 1)} +/- {np.round(std, 1)}')
# -----------------

print('COMPLETED SCRIPT: 9-Plot_Pulse_Experiment.py')
