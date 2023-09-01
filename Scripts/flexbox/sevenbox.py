import numpy as np
from scipy.interpolate import interp1d
from flexbox import Model, Rate


def init_model(input_dict={}):

    compartments = ['atm','tf','ts','ta','ocs','oci','ocd']
    model = Model(compartments)
#==========================================================================
# OBJECTIVE
#   This script calculates rate constants as first order linear processes
#   from the present global Hg budget.  We assume that rate constants are
#   the same between the preindustrial and present-day period. Rates 
#   compiled by Helen Amos and Colin Thackray. 
#
# REFERENCES
#    Amos, H. M., et al. (2014), Global biogeochemical implications of
#      mercury discharges from rivers and sediment burial, Environ. Sci.
#      Technol., 48(16), 9514-9522.
#    Bagnato, E., et al. (2011), New clues on the contribution of earth's
#      volcanism to the global mercury cycle, Bull. Volcanol., 73(5), 497-
#      510.
#    Hararuk, O., et al. (2013), Modeling the sensitivity of soil mercury
#      storage to climate-induced changes in soil carbon pools,
#      Biogeosciences, 10, 2393?2407.
#    Holmes, C. D., et al. (2010), Global atmospheric model for mercury
#      including oxidation by bromine atoms, Atmos. Chem. Phys., 10, 12037
#      -12057.
#    Pirrone, N., et al. (2010), Global mercury emissions to the atmosphere
#      from anthropogenic and natural sources, Atmos. Chem. Phys., 10(13),
#      5951-5964.
#    Selin, N. E., et al. (2008), Global 3-d land-ocean-atmosphere model
#      for mercury: Present-day versus preindustrial cycles and
#      anthropogenic enrichment factors for deposition Global Biogeochem.
#      Cycles, 22(3), GB3099.
#    Smith-Downey, N. V., et al. (2010), Anthropogenic impacts on global
#      storage and emissions of mercury from terrestrial soils: Insights
#      from a new global model, J. Geophys. Res.-Biogeosci., 115, G03008.
#   Soerensen, A. L., et al. (2010), An improved global model for air-sea
#      exchange of mercury: High concentrations over the north atlantic,
#      Environ. Sci. Technol., 44(22), 8574-8580.
#
#==========================================================================

#--------------------------------------------------------------------------
# Reservoirs estimates for present day (Mg)
#--------------------------------------------------------------------------

# Atmosphere
    Ratm       = input_dict.get('Ratm',3900)               # Horowitz et al.(2017)

# Ocean
    Rocs       = input_dict.get('Rocs',2910)               # Soerensen et al. (2010)
    Roci       = input_dict.get('Roci',134000)             # Sunderland and Mason (2007)
    Rocd       = input_dict.get('Rocd',220649)             # Sunderland and Mason (2007)

#  Terrestrial
# instead of 240Gg -> up to 1000Gg
# shift fractions to slower pools?
    Rtf        = input_dict.get('Rtf',1.4*9620)               # Leaf, fast and intermediate pools from Smith-Downey et al (2010)
    Rts        = input_dict.get('Rts',1.4*34900)              # Slow pool from Smith-Downey et al. (2010)
    Rta        = input_dict.get('Rta',1.4*193600)             # Armored pool from Smith-Downey et al. (2010)

#--------------------------------------------------------------------------
# Fluxes for present day (Mg/year)
#--------------------------------------------------------------------------

# Atmosphere
    Dep_oHgII  = input_dict.get('Dep_oHgII',3700)#4600;            # Hg(II) deposition to the ocean (Horowitz et al., 2017)
    Dep_tHgII  = input_dict.get('Dep_tHgII',1900)#1000;            # Hg(II) deposition to land (Horowitz et al., 2017)
    Dep_tHg0   = input_dict.get('Dep_tHg0',1200)                 # Hg(0) deposition to land  (Horowitz et al., 2017)

# Hg0 air - sea exchange
    netEv_Hg0  = input_dict.get('netEv_Hg0',2900)            # net evasion from surface ocean to atmosphere (Horowitz et al., 2017)
    Ev_Hg0_ocs = input_dict.get('Ev_Hg0_ocs',4600)#4600;            # gross ocean evasion to the atmosphere (Horowitz et al., 2017)
    Upt_oHg0   = Ev_Hg0_ocs - netEv_Hg0; # gross uptake of Hg0 from the atmopshere (1700 Mg/yr from Horowitz et al., 2017)

# Surface ocean
    ps_ocs     = input_dict.get('ps_ocs',3320)            # particle settling
    vert_ocsi  = input_dict.get('vert_ocsi',5100)            # gross detrainment flux, surface to intermediate (ref: HMA, v9-01-02, 2005 avg)

# Intermediate ocean
    Kd_vert = input_dict.get('Kd_vert',1e-4)
    k_diffi     = Kd_vert*1e-4*4*3600*24*365/(1000**2) # 1cm2/s diffusivity, ~1000m vertical scale, unit change to per year
    ps_oci     = input_dict.get('ps_oci',480)             # particle settling
    vert_ocis  = input_dict.get('vert_ocis',7100)            # vertical seawater flow, intermediate to surface
    vert_ocid  = input_dict.get('vert_ocid',335)#0.003*Roci#335;             # vertical seawater flow, intermediate to deep
    vert_ocid  += k_diffi*Roci   # add diffusive flux to vertical flow

# Deep ocean
    k_diffd    = Kd_vert*1e-4*4*3600*24*365/(3000**2)
    ps_ocd     = input_dict.get('ps_ocd',210)             # particle settling, burial in deep sediments
    vert_ocdi  = input_dict.get('vert_ocdi',175)#0.001*Rocd#175;             # vertical sea water flow, deep to intermediate
    vert_ocdi += k_diffd*Rocd    # add diffusive flux to vertical flow

# rivers and biomass burning: assuming 75# vegetation (all fast) + 25#
# soils (fast,slow,armored) partitioned by C storage
    fveg       = input_dict.get('fveg',0.95)            # fraction to vegetation
    fsoil      = input_dict.get('fsoil',0.05)            # fraction to soil

# fraction of carbon in each pool (Smith-Downey et al., 2010)
    fCfast     = input_dict.get('fCfast',0.2185)          # fast reservoir
    fCslow     = input_dict.get('fCslow',0.5057)          # slow reservoir
    fCarmored  = input_dict.get('fCarmored',0.2759)          # armored reservoir

# E_Te_rf + E_Te_p + E_Te_rs + E_Te_ra should equal ~1200 (Horowitz et al. 2017)
# Fast terrestrial
    tot_resp = input_dict.get('tot_resp',735.)
    E_Te_rf    = tot_resp*460/735;#460;             # evasion due to respiration of organic carbon
    Te_exfs    = input_dict.get('Te_exfs',325)             # exchange among soil pools, fast pool to slow pool
    Te_exfa    = input_dict.get('Te_exfa',9)               # exchange among soil pools, fast pool to armored pool

# Slow terrestrial
    E_Te_rs    = tot_resp*250/735;#250             # evasion due to respiration of organic carbon
    Te_exsf    = input_dict.get('Te_exsf',205)             # exchange among soil pools, slow pool to fast pool
    Te_exsa    = input_dict.get('Te_exsa',0.5)             # exchange among soil pools, slow pool to armored pool

# Armored terrestrial
    E_Te_ra    = tot_resp*25/735;#25              # evasion due to respiration of organic carbon
    Te_exaf    = input_dict.get('Te_exaf',15)              # exchange among soil pools, armored pool to fast pool
    Te_exam    = input_dict.get('Te_exam',0)               # exchange from armored pool to mineral pool

    E_Te_p = Dep_tHgII + Dep_tHg0 + Dep_oHgII - netEv_Hg0 - E_Te_rf - E_Te_rs - E_Te_ra - 250 - 2270# Budget closing term representing everything else re-emis, ice, etc.

#--------------------------------------------------------------------------
# Atmospheric rates (1/year)
#--------------------------------------------------------------------------
    k_A_oHgII  = Dep_oHgII / Ratm;   
    model.add_rate(Rate('k_A_oHgII', k=k_A_oHgII, compartment_to='ocs', compartment_from='atm',
                        notes='HgII deposition to surface ocean; min = 0.7200 (CDH); max = 0.8913 (ESC)'
                        ))
    k_A_oHg0   = Upt_oHg0  / Ratm; 
    model.add_rate(Rate('k_A_oHg0', k=k_A_oHg0, compartment_to='ocs', compartment_from='atm',
                        notes='gross Hg0 uptake by the surface ocean'
                        ))
    k_A_tHgII  = Dep_tHgII / Ratm;   # 
    k_A_tHg0   = Dep_tHg0  / Ratm;   # Hg0  deposition to terrestrial surfaces;
                                 #   Assumed to represent terrestrial leaf uptake;
                                 #   min = 0.2927 (NES); max = 0.3478 (ESC)
    model.add_rate(Rate('k_A_tfHg0', k=k_A_tHg0, compartment_to='tf', compartment_from='atm',
                        notes='Hg0  deposition to terrestrial surfaces;Assumed to represent terrestrial leaf \
                    uptake; min = 0.2927 (NES); max = 0.3478 (ESC)'
                        ))
# fraction of atmopsheric deposition to...
    fdep_tf    = input_dict.get('fdep_tf',0.5027)             # the fast soil pool
    fdep_ts    = input_dict.get('fdep_ts',0.3213)             # the slow soil pool
    fdep_ta    = input_dict.get('fdep_ta',0.1760)             # the fast armored pool
    model.add_rate(Rate('k_A_tfHgII', k=k_A_tHgII*fdep_tf, compartment_to='tf', compartment_from='atm',
                        notes='HgII deposition to terrestrial fast pool;  min = 0.2927 (NES); max = 0.3043 (ESC)'
                        ))
    model.add_rate(Rate('k_A_tsHgII', k=k_A_tHgII*fdep_ts, compartment_to='ts', compartment_from='atm',
                        notes='HgII deposition to terrestrial slow pool;  min = 0.2927 (NES); max = 0.3043 (ESC)'
                        ))
    model.add_rate(Rate('k_A_taHgII', k=k_A_tHgII*fdep_ta, compartment_to='ta', compartment_from='atm',
                        notes='HgII deposition to terrestrial armored pool;  min = 0.2927 (NES); max = 0.3043 (ESC)'
                        ))

#--------------------------------------------------------------------------
# Surface ocean rates (1/year)
#--------------------------------------------------------------------------
    k_Oc_ev    = Ev_Hg0_ocs / Rocs;  
    model.add_rate(Rate('k_Oc_ev', k=k_Oc_ev, compartment_to='atm', compartment_from='ocs',
                        notes='evasion Hg0 (gross flux); min = 0.7400 (CDH); max = 1.8485 (ESC)'
                        ))
    k_Oc_sp1   = ps_ocs     / Rocs;  
    model.add_rate(Rate('k_Oc_sp1', k=k_Oc_sp1, compartment_to='oci', compartment_from='ocs',
                        notes='particle settling; surface to intermediate'
                        ))
    k_Oc_vsi   = vert_ocsi  / Rocs;  
    model.add_rate(Rate('k_Oc_vsi', k=k_Oc_vsi, compartment_to='oci', compartment_from='ocs',
                        notes='gross detrainment to intermediate ocean'
                        ))

#--------------------------------------------------------------------------
# Intermediate ocean rates (1/year)
#--------------------------------------------------------------------------
    k_Oc_sp2   = ps_oci    / Roci;   
    model.add_rate(Rate('k_Oc_sp2', k=k_Oc_sp2, compartment_to='ocd', compartment_from='oci',
                        notes='particle settling; intermediate to deep'
                        ))
    k_Oc_vis   = vert_ocis / Roci;   
    model.add_rate(Rate('k_Oc_vis', k=k_Oc_vis, compartment_to='ocs', compartment_from='oci',
                        notes='vertical seawater flow, intermediate to surface'
                        ))
    k_Oc_vid   = vert_ocid / Roci;   
    model.add_rate(Rate('k_Oc_vid', k=k_Oc_vid, compartment_to='ocd', compartment_from='oci',
                        notes='vertical seawater flow and diffusion, intermediate to deep'
                        ))

#--------------------------------------------------------------------------
# Deep ocean rates (1/year)
#--------------------------------------------------------------------------
    k_Oc_sp3   = ps_ocd    / Rocd;   
    model.add_rate(Rate('k_Oc_sp3', k=k_Oc_sp3, compartment_to=None, compartment_from='ocd',
                        notes='particle settling, deep ocean burial in deep sediments'
                        ))
    k_Oc_vdi   = vert_ocdi / Rocd;   
    model.add_rate(Rate('k_Oc_vdi', k=k_Oc_vdi, compartment_to='oci', compartment_from='ocd',
                        notes='vertical seawater flow and diffusion, deep to intermediate'
                        ))

#--------------------------------------------------------------------------
# Fast terrestrial reservoir rates (1/year)
#--------------------------------------------------------------------------

# Includes vegetation, ice, and fast+intermediate carbon pools from
# Smith-Downey et al. (2010)

    k_Te_rf    = E_Te_rf / Rtf;      
    model.add_rate(Rate('k_Te_rf', k=k_Te_rf, compartment_to='atm', compartment_from='tf',
                        notes='respiration; fast pool'
                        ))
    k_Te_p     = E_Te_p  / Rtf;      
    model.add_rate(Rate('k_Te_p', k=k_Te_p, compartment_to='atm', compartment_from='tf',
                        notes='photoreduction and re-release of deposited Hg0'
                        ))
    k_T_exfs   = Te_exfs / Rtf;      
    model.add_rate(Rate('k_T_exfs', k=k_T_exfs, compartment_to='ts', compartment_from='tf',
                        notes='exchange among soil pools, fast pool to slow pool'
                        ))
    k_T_exfa   = Te_exfa / Rtf;      
    model.add_rate(Rate('k_T_exfa', k=k_T_exfa, compartment_to='ta', compartment_from='tf',
                        notes='exchange among soil pools, fast pool to armored pool'
                        ))

#--------------------------------------------------------------------------
# Slow terrestrial reservoir rates (1/year)
#--------------------------------------------------------------------------

    k_Te_rs    = E_Te_rs / Rts;       
    model.add_rate(Rate('k_Te_rs', k=k_Te_rs, compartment_to='atm', compartment_from='ts',
                        notes='evasion due to respiration of organic carbon; slow pool'
                        ))
    k_T_exsf   = Te_exsf / Rts;       
    model.add_rate(Rate('k_T_exsf', k=k_T_exsf, compartment_to='tf', compartment_from='ts',
                        notes='exchange among soil pools, slow pool to fast pool'
                        ))
    k_T_exsa   = Te_exsa / Rts;       
    model.add_rate(Rate('k_T_exsa', k=k_T_exsa, compartment_to='ta', compartment_from='ts',
                        notes='exchange among soil pools, slow pool to armored pool'
                        ))

#--------------------------------------------------------------------------
# Armored terrestrial reservoir rates (1/year)
#--------------------------------------------------------------------------

    k_Te_ra    = E_Te_ra / Rta;       
    model.add_rate(Rate('k_Te_ra', k=k_Te_ra, compartment_to='atm', compartment_from='ta',
                        notes='evasion due to respiration of organic carbon; armored pool'
                        ))
    k_T_exaf   = Te_exaf / Rta;       
    model.add_rate(Rate('k_T_exaf', k=k_T_exaf, compartment_to='tf', compartment_from='ta',
                        notes='exchange among soil pools, armored pool to fast pool'
                        ))
    k_T_exam   = Te_exam / Rta;       #  (18 Dec 2011, hma)
    model.add_rate(Rate('k_T_exam', k=k_T_exam, compartment_from='ta',
                        notes='exchange from armored pool to mineral pool'
                        ))

#--------------------------------------------------------------------------
# Rivers
#--------------------------------------------------------------------------
    IHgD_pristine = input_dict.get('IHgD_pristine',78.)
    IHgP_pristine = input_dict.get('IHgP_pristine',659.)

# total discharged to ocean margins
    Te_riv_margin  = IHgD_pristine + IHgP_pristine;

# global fraction of riverine HgP reaching the open oceans (Table 2, Amos et al. 2014)
    Lriver_FHgP = 'WeightedWalsh'
    if Lriver_FHgP=='Walsh':
        f_HgPexport    = 0.30;
    elif Lriver_FHgP=='Chester':
        f_HgPexport    = 0.10;
    elif Lriver_FHgP=='WeightedWalsh':
        f_HgPexport   = 0.07
    f_HgPexport = input_dict.get('f_HgPexport',f_HgPexport)


# total reaching the open ocean
    Te_riv_ocean   = IHgD_pristine + f_HgPexport*IHgP_pristine;

# First-order rate coefficients (1/yr)

# Riverine discharge of terrestrial Hg to ocean margins
    k_T_riv_f      = (Te_riv_margin*fveg + (Te_riv_margin*fsoil*fCfast)) / Rtf;  # fast
    k_T_riv_s      = (Te_riv_margin*fsoil*fCslow) / Rts;                         # slow
    k_T_riv_a      = (Te_riv_margin*fsoil*fCarmored) / Rta;                      # armored

    k_O_riv_f      = (Te_riv_ocean*fveg + (Te_riv_ocean*fsoil*fCfast)) / Rtf;    
    model.add_rate(Rate('k_O_riv_f', k=k_O_riv_f, compartment_to='ocs', compartment_from='tf',
                        notes='riverine discharge of terrestrial Hg to open ocean; fast pool'
                        ))
    k_O_riv_s      = (Te_riv_ocean*fsoil*fCslow) / Rts;                          
    model.add_rate(Rate('k_O_riv_s', k=k_O_riv_s, compartment_to='ocs', compartment_from='ts',
                        notes='riverine discharge of terrestrial Hg to open ocean; slow pool'
                        ))
    k_O_riv_a      = (Te_riv_ocean*fsoil*fCarmored) / Rta;                       
    model.add_rate(Rate('k_O_riv_a', k=k_O_riv_a, compartment_to='ocs', compartment_from='ta',
                        notes='riverine discharge of terrestrial Hg to open ocean; armored pool'
                        ))

    k_L_riv_f      = k_T_riv_f - k_O_riv_f
    model.add_rate(Rate('k_L_riv_f', k=k_L_riv_f, compartment_from='tf',
                        notes='riverine discharge of terrestrial Hg to ocean margin sediment; fast pool'
                        ))
    k_L_riv_s      = k_T_riv_s - k_O_riv_s
    model.add_rate(Rate('k_L_riv_s', k=k_L_riv_s, compartment_from='ts',
                        notes='riverine discharge of terrestrial Hg to ocean margin sediment; slow pool'
                        ))
    k_L_riv_a      = k_T_riv_a - k_O_riv_a
    model.add_rate(Rate('k_L_riv_a', k=k_L_riv_a, compartment_from='ta',
                        notes='riverine discharge of terrestrial Hg to ocean margin sediment; armored pool'
                        ))

    return model
