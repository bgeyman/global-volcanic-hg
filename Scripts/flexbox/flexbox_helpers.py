import numpy as np
import pandas as pd
import flexbox

# -------------------------------------------------------------------------------------------
# --- Get steady state reservoir concentrations (and print)
# -------------------------------------------------------------------------------------------

def print_steady(model, volcanic_air=int, volcanic_ocean=int, print_output=True):
    # Geogenic emissions array
    Egeo = np.zeros(7)
    Egeo[0] = volcanic_air # Mg/y to air
    Egeo[-1] = volcanic_ocean # Mg/y to deep ocean
    steady = model.get_steady_state(Egeo)

    if print_output==True:
        print('---')
        print('Steady State Reservoir Masses')
        print('---')
        for (i, r) in enumerate(model.compartments):
            print(r, f'{steady[i]:.0f} Mg')
        print('---')
    return steady

# -------------------------------------------------------------------------------------------
# --- Import Metastring
# -------------------------------------------------------------------------------------------
def import_metastring(weights_dict={}):
    '''
        Arguments: 
            1. `weights_dict` : a dictionary with weights corresponding to the following keys 
                                [f_tf, f_ts, f_ta, 'w_tfm', 'w_tsm', 'w_tam','f_sequestered']
                                example: weights_dict = {'f_tf':0.2, 'f_ts':0.3, 'f_ta':0.5,
                                                         'w_tfm':0.1, 'w_tsm':0.2, 'w_tam':0.7,
                                                          'f_sequestered':0.4}
    '''

    # ---- 
    # initialize weights
    f_tf = weights_dict['f_tf']
    f_ts = weights_dict['f_ts']
    f_ta = weights_dict['f_ta']
    w_tfm = weights_dict['w_tfm']
    w_tsm = weights_dict['w_tsm']
    w_tam = weights_dict['w_tam']

    fraction_sequestered = weights_dict['f_sequestered'] # for land/water emissions
    frac = 1 - fraction_sequestered

    w_tf = f_tf*frac
    w_ts = f_ts*frac
    w_ta = f_ta*frac + fraction_sequestered 

    # check that weights add up to 1
    assert np.abs(w_tf+w_ts+w_ta-1.0) < 0.01
    assert np.abs(w_tfm+w_tsm+w_tam-1.0) < 0.01

    # end weight initialization

    metastring = f"""
    name: test run
    emissions_filename: 'inputs/AnthroPost1510.csv'
    fragments: 
        #tag:          ## DEFAULT VALUES CAN BE LEFT BLANK
        #    filename: 'path/to/file.csv' # DEFAULTS TO ABOVE FILENAME
        #    column_name: 'COLUMN1'
        #    time_column: 'Year' # DEFAULT
        #    interpolation: 'linear' # DEFAULT
        #    csvheader: 1 # DEFAULT
        #    compartment: 0 # DEFAULT
        #    weight: 1. # DEFAULT
        #
        # Emissions to air before this inventory:
        # Natural emissions (this study)
        Volcanic EXP:
            filename: './inputs/Volcanic_Hg.csv'
            column_name: 'EXP'
            csvheader: 0
            numboxes: 7
        Volcanic EFF:
            filename: './inputs/Volcanic_Hg.csv'
            column_name: 'EFF'
            csvheader: 0
            numboxes: 7
        Volcanic PVF:
            filename: './inputs/Volcanic_Hg.csv'
            column_name: 'PVF'
            csvheader: 0
            numboxes: 7
        # Pre-1510 emissions:
        Preindustrial:
            filename: './inputs/AnthroPre1510.csv'
            column_name: 'Air'
            csvheader: 0
            numboxes: 7
        # Total world emissions:
        World air:
            column_name: 'Air'
            numboxes: 7
        World tf:
            column_name: 'No Mining'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        World ts:
            column_name: 'No Mining'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        World ta:
            column_name: 'No Mining'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        Mining tf:
            column_name: 'Mining'
            weight: {w_tfm}
            compartment: 1
            numboxes: 7
        Mining ts:
            column_name: 'Mining'
            weight: {w_tsm}
            compartment: 2
            numboxes: 7
        Mining ta:
            column_name: 'Mining'
            weight: {w_tam}
            compartment: 3
            numboxes: 7
        Mining tf alt:
            column_name: 'Mining'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        Mining ts alt:
            column_name: 'Mining'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        Mining ta alt:
            column_name: 'Mining'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        Landfill tf:
            column_name: 'Landfill'
            time_column: 'LFYear'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        Landfill ts:
            column_name: 'Landfill'
            time_column: 'LFYear'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        Landfill ta:
            column_name: 'Landfill'
            time_column: 'LFYear'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        Landfill Seq:
            column_name: 'Landfill'
            time_column: 'LFYear'
            weight: 1.0
            compartment: 3
            numboxes: 7
            # These are the gold and silver emissions for North and South America
        NA Gold->air:
            column_name: 'NA Gold'
            compartment: 0
            numboxes: 7
        NAAg tf:
            column_name: 'North America'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        NAAg ts:
            column_name: 'North America'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        NAAg ta:
            column_name: 'North America'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        NAAg tf seq:
            column_name: 'North America'
            weight: {w_tfm}
            compartment: 1
            numboxes: 7
        NAAg ts seq:
            column_name: 'North America'
            weight: {w_tsm}
            compartment: 2
            numboxes: 7
        NAAg ta seq:
            column_name: 'North America'
            weight: {w_tam}
            compartment: 3
            numboxes: 7
        SAAg tf:
            column_name: 'South America'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        SAAg ts:
            column_name: 'South America'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        SAAg ta:
            column_name: 'South America'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        Pre1850 air:
            column_name: 'AIRpre1850'
            time_column: 'YearSplit'
            numboxes: 7
        Pre1850 tf:
            column_name: 'LWpre1850'
            time_column: 'YearSplit'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        Pre1850 ts:
            column_name: 'LWpre1850'
            time_column: 'YearSplit'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        Pre1850 ta:
            column_name: 'LWpre1850'
            time_column: 'YearSplit'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        1850-1900 air:
            column_name: 'AIR1850-1900'
            time_column: 'YearSplit'
            numboxes: 7
        1850-1900 tf:
            column_name: 'LW1850-1900'
            time_column: 'YearSplit'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        1850-1900 ts:
            column_name: 'LW1850-1900'
            time_column: 'YearSplit'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        1850-1900 ta:
            column_name: 'LW1850-1900'
            time_column: 'YearSplit'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        1900-1950 air:
            column_name: 'AIR1900-1950'
            time_column: 'YearSplit'
            numboxes: 7
        1900-1950 tf:
            column_name: 'LW1900-1950'
            time_column: 'YearSplit'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        1900-1950 ts:
            column_name: 'LW1900-1950'
            time_column: 'YearSplit'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        1900-1950 ta:
            column_name: 'LW1900-1950'
            time_column: 'YearSplit'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        1950-1990 air:
            column_name: 'AIR1950-1990'
            time_column: 'YearSplit'
            numboxes: 7
        1950-1990 tf:
            column_name: 'LW1950-1990'
            time_column: 'YearSplit'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        1950-1990 ts:
            column_name: 'LW1950-1990'
            time_column: 'YearSplit'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        1950-1990 ta:
            column_name: 'LW1950-1990'
            time_column: 'YearSplit'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        1990-2000 air:
            column_name: 'AIR1990-2000'
            time_column: 'YearSplit'
            numboxes: 7
        1990-2000 tf:
            column_name: 'LW1990-2000'
            time_column: 'YearSplit'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        1990-2000 ts:
            column_name: 'LW1990-2000'
            time_column: 'YearSplit'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        1990-2000 ta:
            column_name: 'LW1990-2000'
            time_column: 'YearSplit'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        2000-2010 air:
            column_name: 'AIR2000-2010'
            time_column: 'YearSplit'
            numboxes: 7
        2000-2010 tf:
            column_name: 'LW2000-2010'
            time_column: 'YearSplit'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        2000-2010 ts:
            column_name: 'LW2000-2010'
            time_column: 'YearSplit'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        2000-2010 ta:
            column_name: 'LW2000-2010'
            time_column: 'YearSplit'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        Tambora:
            column_name: 'Tambora'
            time_column: 'YearVolc'
            numboxes: 7
        Krakatoa:
            column_name: 'Krakatoa'
            time_column: 'YearVolc'
            numboxes: 7
        Pinatubo:
            column_name: 'Pinatubo'
            time_column: 'YearVolc'
            numboxes: 7
        NALF tf:
            column_name: 'LandfillNA'
            time_column:  'LandfillYear'
            compartment: 1
            weight: {w_tf}
            numboxes: 7
        NALF ts:
            column_name: 'LandfillNA'
            time_column:  'LandfillYear'
            compartment: 2
            weight: {w_ts}
            numboxes: 7
        NALF ta:
            column_name: 'LandfillNA'
            time_column:  'LandfillYear'
            compartment: 3
            weight: {w_ta}
            numboxes: 7
        NA air:
            column_name: 'NA Air'
            compartment: 0
            numboxes: 7
        NA tf:
            column_name: 'NA Land/Water'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        NA ts:
            column_name: 'NA Land/Water'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        NA ta:
            column_name: 'NA Land/Water'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        NA Gold->air:
            column_name: 'NA Gold'
            compartment: 0
            numboxes: 7
        NA Gold->tf:
            column_name: 'NA Gold'
            compartment: 1
            numboxes: 7
        NA Gold->ts:
            column_name: 'NA Gold'
            compartment: 2
            numboxes: 7
        NA Gold->ta:
            column_name: 'NA Gold'
            compartment: 3
            numboxes: 7
        Gold air:
            column_name: 'WorldGoldAir'
            time_column: 'GoldYear'
            numboxes: 7
        Gold tf:
            column_name: 'WorldGoldLW'
            time_column: 'GoldYear'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        Gold ts:
            column_name: 'WorldGoldLW'
            time_column: 'GoldYear'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        Gold ta:
            column_name: 'WorldGoldLW'
            time_column: 'GoldYear'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        Coal air:
            column_name: 'Coal Air'
            time_column: 'Coal Year'
            numboxes: 7
        Coal tf:
            column_name: 'Coal LW'
            time_column: 'Coal Year'
            weight: {w_tf}
            compartment: 1
            numboxes: 7
        Coal ts:
            column_name: 'Coal LW'
            time_column: 'Coal Year'
            weight: {w_ts}
            compartment: 2
            numboxes: 7
        Coal ta:
            column_name: 'Coal LW'
            time_column: 'Coal Year'
            weight: {w_ta}
            compartment: 3
            numboxes: 7
        Minamata air:
            column_name: 'Minamata'
            time_column: 'MinamataYear'
            compartment: 0
            numboxes: 7
        Minamata tf:
            column_name: 'MinamataLW'
            time_column: 'MinamataYear'
            compartment: 1
            weight: {w_tf}
            numboxes: 7
    """
    meta = flexbox.Meta(metastring)

    return meta