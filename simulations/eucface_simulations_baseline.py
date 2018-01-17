#!/usr/bin/env python

""" EucFACE CO2 simulations

Full spin-up and simulations under amb and ele CO2 conditions for quasi-equil analysis
"""

import os
import shutil
import sys
import subprocess

USER = os.getlogin()
sys.path.append('/Users/%s/Documents/Research/Projects/eucface/Git/scripts' % (USER))
import adjust_gday_param_file as ad


__author__  = "Martin De Kauwe"
__version__ = "1.0 (14.12.2014)"
__email__   = "mdekauwe@gmail.com"

def main(experiment_id, site, SPIN_UP=True, POST_INDUST=True, ELE_INITIALIZATION=True, ELE_SPINUP=True, ELE_EQUILIB=True):

    GDAY_SPIN = "./gday -s -p "
    GDAY = "./gday -p "

    # dir names
    base_param_name = "base_start_with_P"
    base_param_dir = "/Users/%s/Documents/Research/Projects/eucface/Git/GDAY/params" % (USER)
    base_dir = os.path.dirname(os.getcwd())
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "met_data")
    run_dir = os.path.join(base_dir, "outputs")

    if SPIN_UP == True:

        # copy base files to make two new experiment files
        shutil.copy(os.path.join(base_param_dir, base_param_name + ".cfg"),
                    os.path.join(param_dir, "%s_%s_model_spinup.cfg" % \
                                                (experiment_id, site)))

        # Run model to equilibrium assuming forest, growing C pools from
        # effectively zero
        itag = "%s_%s_model_spinup" % (experiment_id, site)
        otag = "%s_%s_model_spunup" % (experiment_id, site)
        mtag = "%s_met_data_amb_var_co2.csv" % (site)
        out_fn = itag + "_equilib.out"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {
                        # files
                        "out_param_fname": "%s" % (out_param_fname),
                        "cfg_fname": "%s" % (cfg_fname),
                        "met_fname": "%s" % (met_fname),
                        "out_fname": "%s" % (out_fname),

                        # default C:N 25.
                        # Canopy height = 22 m average of 6 plots at UWS, site_description_stuff/EucFACE_Plot_Summary.doc
                        "activesoil": "0.001",
                        "activesoiln": "0.00004",
                        "activesoilp": "0.000002",
                        "age": "100.0",
                        "branch": "0.001",
                        "branchn": "0.00004",
                        "branchp": "0.000002",
                        "cstore": "0.0",
                        "nstore": "0.0",
                        "pstore": "0.0",
                        "inorgn": "0.0000",     # 0.00004
                        "inorglabp": "0.0000",  # 0.00004
                        "inorgsorbp": "0.0",
                        "inorgssorbp": "0.0",
                        "inorgoccp": "0.0",
                        "inorgparp": "0.054",
                        "metabsoil": "0.0",
                        "metabsoiln": "0.0",
                        "metabsoilp": "0.0",
                        "metabsurf": "0.0",
                        "metabsurfn": "0.0",
                        "metabsurfp": "0.0",
                        "passivesoil": "0.001",
                        "passivesoiln": "0.0004",
                        "passivesoilp": "0.000002",
                        "prev_sma": "1.0",
                        "root": "0.001",
                        "croot": "0.0",   # don't simulate coarse roots
                        "crootn": "0.0",  # don't simulate coarse roots
                        "crootp": "0.0",  # don't simulate coarse roots
                        "rootn": "0.00004",
                        "rootp": "0.000002",
                        "sapwood": "0.001",
                        "shoot": "0.001",
                        "shootn": "0.00004",
                        "shootp": "0.000002",
                        "slowsoil": "0.001",
                        "slowsoiln": "0.00004",
                        "slowsoilp": "0.000002",
                        "stem": "0.001",
                        "stemn": "0.00004",
                        "stemp": "0.000002",
                        "stemnimm": "0.00004",
                        "stempimm": "0.000002",
                        "stemnmob": "0.0",
                        "stempmob": "0.0",
                        "structsoil": "0.001",
                        "structsoiln": "0.00004",
                        "structsoilp": "0.000002",
                        "structsurf": "0.001",
                        "structsurfn": "0.00004",
                        "structsurfp": "0.0000024",

                        # parameters
                        "resp_coeff": "0.2",      
                        "alpha_j": "0.308",  # Taking the theoretical maximum (from Belinda) 0.385 x 0.8 (leaf absorptance) = 0.308
                        "intercep_frac": "0.15",
                        "max_intercep_lai": "3.0",
                        "latitude": "-33.61",
                        "albedo": "0.2",
                        "finesoil": "0.2",   # silt + clay fraction. Surface soil texture (upper 45 cm) for Clarenden sand: 80 +/- 8% sand, 9 +/- 5% silt, 11 +/- 3% clay
                        "slamax": "5.1",    # current unit: m2 kg-1; original unit: 43.7 +/- 1.5 cm2 g 1 dry mass
                        "sla": "5.1",       # current unit: m-2 kg-1; original unit: 43.7 +/-  1.5 cm2 g 1 dry mass
                        "slazero": "5.1",   # current unit: m-2 kg-1; original unit: 43.7+/-  1.5 cm2 g 1 dry mass
                        "lai_closed": "0.5",  # I am effectively turning this feature off by setting it so low
                        "c_alloc_fmax": "0.45",  # 0.35
                        "c_alloc_fmin": "0.05",  # 0.15
                        "c_alloc_rmax": "0.45",  # 0.35
                        "c_alloc_rmin": "0.05",  # 0.05
                        "c_alloc_bmax": "0.1",   # 0.1
                        "c_alloc_bmin": "0.1",   # 0.1
                        "c_alloc_cmax": "0.0", # turn off coarse roots!
                        "biochemical_p_constant": "150.0",
                        "fretrans": "0.5",
                        "fretransp": "0.5",
                        "rretrans": "0.0",
                        "bretrans": "0.0",
                        "wretrans": "0.7",
                        "cretrans": "0.0",
                        "crit_n_cost_of_p": "15.0",
                        "ncwnewz": "0.003",          #New stem ring N:C at zero leaf N:C (mobile)
                        "ncwnew": "0.003",           #New stem ring N:C at critical leaf N:C (mob)
                        "ncwimmz": "0.003",          #Immobile stem N C at zero leaf N C
                        "ncwimm": "0.003",           #Immobile stem N C at critical leaf N C
                        "ncbnewz": "0.003",          #new branch N C at zero leaf N C
                        "ncbnew": "0.003",           #new branch N C at critical leaf N C
                        "nccnewz": "0.003",          #new coarse root N C at zero leaf N C
                        "nccnew": "0.003",           #new coarse root N C at critical leaf N C
                        "ncrfac": "0.8",
                        "ncmaxfyoung": "0.04",
                        "ncmaxfold": "0.04",
                        "ncmaxr": "0.03",
                        "retransmob": "0.0",
                        "fdecay": "0.6",    # 18 mth turnover * 1/30
                        "fdecaydry": "0.6", # 18 mth turnover * 1/30
                        "max_p_biochemical": "0.001",
                        "rdecay": "0.6",
                        "rdecaydry": "0.6",
                        "crdecay": "0.00",           # turn off coarse roots!
                        "bdecay": "0.1",            # no idea, assuming 50 years
                        "wdecay": "0.1",            # no idea, assuming 50 years
                        "watdecaydry": "0.0",
                        "watdecaywet": "0.1",
                        "ligshoot": "0.18",          # Based on white et al. 2000 #"0.145",   # assuming leaf and root same as DE word document
                        "ligroot": "0.22",           # Based on white et al. 2000    # assuming leaf and root same as DE word document
                        "rateuptake": "1.8",
                        "rateloss": "0.05",           # was 0.1
                        "topsoil_depth": "450.0",    # Not needed as I have supplied the root zone water and topsoil water available
                        "rooting_depth": "2500.0",   # Not needed as I have supplied the root zone water and topsoil water available
                        "wcapac_root": "300.0",      # [mm] (FC-WP)*rooting_depth. But using 2.0 m, site_description_stuff/EucFACE_Plot_Summary.doc
                        "wcapac_topsoil": "67.5",    # [mm] (FC-WP)*rooting_depth. But using 0.45 m, site_description_stuff/EucFACE_Plot_Summary.doc
                        "ctheta_topsoil": "0.65",     # Derive based on soil type loamy_sand
                        "ntheta_topsoil": "8.0",     # Derive based on soil type loamy_sand
                        "ctheta_root": "0.525",      # Derive based on soil type sandy_clay_loam
                        "ntheta_root": "5.5",        # Derive based on soil type sandy_clay_loam
                        "topsoil_type": "loamy_sand",
                        "rootsoil_type": "sandy_clay_loam",
                        "soil_order": "andisol",
                        "ks": "0.5",
                        "kp": "0.3",
                        "krp": "0.00001",
                        #"dz0v_dh": "0.1",
                        #"z0h_z0m": "1.0",
                        #"displace_ratio": "0.67",

                        "dz0v_dh": "0.05",         # Using Value from JULES for TREE PFTs as I don't know what is best. However I have used value from Jarvis, quoted in Jones 1992, pg. 67. Produces a value within the bounds of 3.5-1.1 mol m-2 s-1 Drake, 2010, GCB for canht=17
                        "displace_ratio": "0.75",  # From Jones, pg 67, following Jarvis et al. 1976
                        "z0h_z0m": "1.0",

                        "g1": "3.8667",          # 3.8667 Fit by Me to Teresa's data 7th Nov 2013; or 2.78 from stomatal model
                        #"jmaxna": "14.891",      # 
                        #"jmaxpa": "291.4305",    # 
                        #"jmaxnb": "99.497",      # 
                        #"jmaxpb": "99.949",      # 88.56  
                        #"vcmaxna": "10.453",     # 6.426
                        #"vcmaxpa": "153.1748",    
                        #"vcmaxnb": "74.522",     # 60.526
                        #"vcmaxpb": "57.242",     # 27.66
                        "jmaxna": "49.930",      # forcing intercept to zero; if use all species df, 49.743
                        "jmaxpa": "933.90",      # forcing intercept to zero; if use all species df, 842.46 
                        "jmaxnb": "0.0",         # forcing intercept to zero
                        "jmaxpb": "0.0",         # forcing intercept to zero
                        "vcmaxna": "27.707",     # forcing intercept to zero; if use all species df, 27.627
                        "vcmaxpa": "516.83",     # forcing intercept to zero; if use all species df, 468.76
                        "vcmaxnb": "0.0",        # forcing intercept to zero
                        "vcmaxpb": "0.0",        # forcing intercept to zero
                        "measurement_temp": "25.0", # parameters obtained at 22 not 25 degrees
                        "heighto": "4.826",
                        "htpower": "0.35",
                        "height0": "5.0",
                        "height1": "25.0",
                        "leafsap0": "4000.0",     # "4000.0",
                        "leafsap1": "2700.0",     # 2700
                        "branch0": "5.61",
                        "branch1": "0.346",
                        "croot0": "0.34",
                        "croot1": "0.84",
                        "targ_sens": "0.5",
                        "density": "800.0",       # 480
                        "nf_min": "0.005", 
                        "nf_crit": "0.015",
                        "sapturnover": "0.1",
                        "p_atm_deposition": "0.000086",   # 1/4 of value from Table 4, Olander et al. 2005; Earth Interactions.
                        "p_rate_par_weather": "0.0001", # Calcualted so that weathering rate = atm deposition;
                        "passpcmin": "0.005",
                        "passpcmax": "0.05",
                        "psecmnp": "0.000022",
                        "pcbnew": "0.0003",
                        "pcbnewz": "0.0003",
                        "pccnew": "0.0003",
                        "pccnewz": "0.0003",
                        "pcmaxfold": "0.002",    # 0.0015 Table 3, Olander et al. 2005, Earth Interactions.
                        "pcmaxfyoung": "0.002",
                        "pcmaxr": "0.0006",
                        "pcrfac": "0.8",
                        "pcwimm": "0.00014",
                        "pcwimmz": "0.00014",
                        "pcwnew": "0.00014",
                        "pcwnewz": "0.00014",
                        "pf_crit": "0.002",
                        "pf_min": "0.0002",
                        "phmax": "7.6",
                        "phmin": "5.0",
                        "phtextmin": "0.000008",
                        "phtextmax": "0.00015",
                        "phtextslope": "0.00004",
                        "pmax": "0.002",
                        "pmin": "0.01",
                        "pmin0": "0.0",
                        "pmincrit": "2.0",
                        "prateloss": "0.05",
                        "prateuptake": "3.6",    # Fitted value to obtain balance between uptake N:P ratio and reasonable P labile pool
                        "slowpcmin": "0.005",
                        "slowpcmax": "0.011111",
                        "soilph": "4.5",          # Olander et al., 2005, Earth Interactions.
                        "sorpmx": "5.0",
                        "sorpaf": "1.0",
                        "structcp": "5500.0",
                        "structratp": "0.0",

                        # control
                        "adjust_rtslow": "false",  # priming, off
                        "alloc_model": "allometric",
                        "assim_model": "mate",
                        "calc_sw_params": "true",   #false=use fwp values, true=derive them
                        "deciduous_model": "false",
                        "disturbance": "false",
                        "exudation": "false",
                        "fixed_stem_nc": "true",
                        "fixed_stem_pc": "true",
                        "fixleafnc": "false",
                        "fixleafpc": "false",
                        "grazing": "false",
                        "gs_model": "medlyn",
                        "aci_relationship": "baseline",
                        "model_optroot": "false",
                        "modeljm": "1",
                        "ncycle": "true",
                        "pcycle": "false",
                        "nuptake_model": "1",
                        "puptake_model": "1",
                        "triose_p": "false",
                        "output_ascii": "true",
                        "passiveconst": "false",
                        "print_options": "end",
                        "ps_pathway": "c3",
                        "respiration_model": "fixed",
                        "strfloat": "0",
                        "strpfloat": "0",
                        "sw_stress_model": "1",  # Sands and Landsberg
                        "use_eff_nc": "0",
                        "text_effect_p": "1",
                        "water_stress": "true",

        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY_SPIN + cfg_fname)
        
        
    if POST_INDUST == True:

        # copy spunup base files to make two new experiment files
        shutil.copy(os.path.join(param_dir, "%s_%s_model_spunup.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_spunup_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_spunup_adj" % (experiment_id, site)
        otag = "%s_%s_model_indust" % (experiment_id, site)
        mtag = "%s_met_data_amb_var_co2.csv" % (site)
        out_fn = "%s_amb_equilib.csv" % (site)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
     
                         # control
                         "print_options": "end",
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)

    if POST_INDUST == True:

        # copy spunup base files to make two new experiment files
        shutil.copy(os.path.join(param_dir, "%s_%s_model_spunup.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_spunup_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_spunup_adj" % (experiment_id, site)
        otag = "%s_%s_model_indust" % (experiment_id, site)
        mtag = "%s_met_data_amb_var_co2.csv" % (site)
        out_fn = "%s_amb_equilib.csv" % (site)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
     
                         # control
                         "print_options": "daily",
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
    
    # elevated co2 initialization to store output 
    if ELE_INITIALIZATION == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_indust.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_indust_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_indust_adj" % (experiment_id, site)
        otag = "%s_%s_model_ele_initial" % (experiment_id, site)
        mtag = "%s_met_data_%s_var_co2.csv" % (site, treatment)
        out_fn = "%s_ele_initial.csv" % (site)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
    
    # elevated co2 initialization to store cfg
    if ELE_INITIALIZATION == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_indust.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_indust_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_indust_adj" % (experiment_id, site)
        otag = "%s_%s_model_ele_initial" % (experiment_id, site)
        mtag = "%s_met_data_%s_var_co2.csv" % (site, treatment)
        out_fn = "%s_ele_initial.csv" % (site)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
    
    # elevated CO2 during spin up to store cfg file
    if ELE_SPINUP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_ele_initial.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_ele_spinup.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_ele_spinup" % (experiment_id, site)
        otag = "%s_%s_model_ele_spunup" % (experiment_id, site)
        mtag = "%s_met_data_%s_var_co2.csv" % (site, treatment)
        out_fn = "FACE_EUC_ele_spunup_%s%s.csv" % (site, treatment.upper())
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY_SPIN + cfg_fname)
    
    # elevated co2 final equilibrium simulation
    if ELE_EQUILIB == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_ele_spunup.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_ele_equil.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_ele_equil" % (experiment_id, site)
        otag = "%s_%s_model_ele_final" % (experiment_id, site)
        mtag = "%s_met_data_%s_var_co2.csv" % (site, treatment)
        out_fn = "%s_ele_final_equilib.csv" % (site)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)

if __name__ == "__main__":

    experiment_id = "FACE"
    site = "EUC"
    treatment = "ele"
    main(experiment_id, site, SPIN_UP=True, POST_INDUST=True, ELE_INITIALIZATION=True, ELE_SPINUP=True, ELE_EQUILIB=True)
