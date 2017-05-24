#!/usr/bin/env python
# coding: utf-8
""" Translate GDAY output file

Match the NCEAS format and while we are at it carry out unit conversion so that
we matched required standard. Data should be comma-delimited
"""
import shutil
import os
import numpy as np
import csv
import sys
#import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
#from io import StringIO
from io import BytesIO

__author__  = "Martin De Kauwe"
__version__ = "1.0 (12.05.2014)"
__email__   = "mdekauwe@gmail.com"

def date_converter(*args):
    return dt.datetime.strptime(str(int(float(args[0]))) + " " +\
                                str(int(float(args[1]))), '%Y %j')

def translate_output(infname, met_fname):
    outdir = "../outputs"
    UNDEF = -9999.
    units = setup_units()
    variable, variable_names = setup_varnames()

    # load met stuff, i.e. the stuff needed for NCEAS output that G'day
    # does not output
    envir = load_met_input_data(met_fname)

    # load the rest of the g'day output
    (gday, git_ver) = load_gday_output(infname)

    # merge dictionaries to ease output
    data_dict = dict(envir, **gday)

    ofname = os.path.join(outdir, "temp.nceas")
    f = open(ofname, "w")
    f.write("%s" % (git_ver))

    # write output in csv format
    writer = csv.writer(f, dialect=csv.excel, lineterminator="\n")
    writer.writerow(variable)
    writer.writerow(units)
    writer.writerow(variable_names)
    for i in range(len(gday['DOY'])):
        writer.writerow([("%.8f" % (float(data_dict[k][i])) \
                         if k in data_dict else UNDEF)
                         for k in variable_names])

    # Need to replace the temp file with the infname which is actually
    # the filename we want to use
    shutil.move(ofname, infname)

def remove_comments_from_header(fname):
    """ I have made files with comments which means the headings can't be
    parsed to get dictionary headers for pandas! Solution is to remove these
    comments first """
    #s = StringIO()
    s = BytesIO()
    with open(fname) as f:
        for line in f:
            if '#' in line:
                line = line.replace("#", "").lstrip(' ')
            s.write(line)
    s.seek(0) # "rewind" to the beginning of the StringIO object

    return s

def remove_comments_from_header_and_get_git_rev(fname):
    """ I have made files with comments which means the headings can't be
    parsed to get dictionary headers for pandas! Solution is to remove these
    comments first """
    #s = StringIO()
    s = BytesIO()
    with open(fname) as f:
        line_counter = 0
        for line in f:
            if line_counter == 0:
                git_ver = line.rstrip(' ')
            if '#' in line:
                line = line.replace("#", "").lstrip(' ')
            s.write(line)
            line_counter += 1
    s.seek(0) # "rewind" to the beginning of the StringIO object

    return s, git_ver

def load_met_input_data(fname):
    MJ_TO_MOL = 4.6
    SW_TO_PAR = 0.48
    DAYS_TO_HRS = 24.0
    UMOL_TO_MOL = 1E-6
    tonnes_per_ha_to_g_m2 = 100.0

    s = remove_comments_from_header(fname)
    met_data = pd.read_csv(s, parse_dates=[[0,1]], skiprows=4, index_col=0,
                           sep=",", keep_date_col=True,
                           date_parser=date_converter)

    precip = met_data["rain"]
    par = (met_data["par_am"] + met_data["par_pm"]) * MJ_TO_MOL
    air_temp = met_data["tair"]
    soil_temp = met_data["tsoil"]
    vpd = (met_data["vpd_am"] + met_data["vpd_pm"]) / 2.0
    co2 = met_data["co2"]
    ndep = met_data["ndep"] * tonnes_per_ha_to_g_m2

    return {'CO2': co2, 'PREC':precip, 'PAR':par, 'TAIR':air_temp, 'TSOIL':soil_temp,
            'VPD':vpd, 'NDEP':ndep}

def load_gday_output(fname):
    SW_RAD_TO_PAR = 2.3
    UNDEF = -9999.
    tonnes_per_ha_to_g_m2 = 100
    yr_to_day = 365.25

    (s, git_ver) = remove_comments_from_header_and_get_git_rev(fname)
    out = pd.read_csv(s, parse_dates=[[0,1]], skiprows=1, index_col=0,
                      sep=",", keep_date_col=True, date_parser=date_converter)

    year = out["year"]
    doy = out["doy"]

    # state outputs
    pawater_root = out["pawater_root"]
    shoot = out["shoot"] * tonnes_per_ha_to_g_m2
    stem = out["stem"] * tonnes_per_ha_to_g_m2
    branch = out["branch"] * tonnes_per_ha_to_g_m2
    fine_root = out["root"] * tonnes_per_ha_to_g_m2
    coarse_root = out["croot"] * tonnes_per_ha_to_g_m2
    coarse_rootn = out["crootn"] * tonnes_per_ha_to_g_m2
    litterc = out["litterc"] * tonnes_per_ha_to_g_m2
    littercag = out["littercag"] * tonnes_per_ha_to_g_m2
    littercbg = out["littercbg"] * tonnes_per_ha_to_g_m2
    soilc = out["soilc"] * tonnes_per_ha_to_g_m2
    lai = out["lai"]
    shootn = out["shootn"] * tonnes_per_ha_to_g_m2
    stemn = out["stemn"] * tonnes_per_ha_to_g_m2
    branchn = out["branchn"] * tonnes_per_ha_to_g_m2
    rootn = out["rootn"] * tonnes_per_ha_to_g_m2
    crootn = out["crootn"] * tonnes_per_ha_to_g_m2
    litternag = out["litternag"] * tonnes_per_ha_to_g_m2
    litternbg = out["litternbg"] * tonnes_per_ha_to_g_m2
    nsoil = out["soiln"] * tonnes_per_ha_to_g_m2
    inorgn = out["inorgn"] * tonnes_per_ha_to_g_m2
    tnc = out["cstore"] * tonnes_per_ha_to_g_m2
    nstorage = out["nstore"] * tonnes_per_ha_to_g_m2
    pstorage = out["pstore"] * tonnes_per_ha_to_g_m2
    activesoiln = out["activesoiln"] * tonnes_per_ha_to_g_m2
    slowsoiln = out["slowsoiln"] * tonnes_per_ha_to_g_m2
    passivesoiln = out["passivesoiln"] * tonnes_per_ha_to_g_m2
    npoolo = activesoiln + slowsoiln + passivesoiln
    shootp = out["shootp"] * tonnes_per_ha_to_g_m2
    stemp = out["stemp"] * tonnes_per_ha_to_g_m2
    branchp = out["branchp"] * tonnes_per_ha_to_g_m2
    rootp = out["rootp"] * tonnes_per_ha_to_g_m2
    crootp = out["crootp"] * tonnes_per_ha_to_g_m2
    litterpag = out["litterpag"] * tonnes_per_ha_to_g_m2
    litterpbg = out["litterpbg"] * tonnes_per_ha_to_g_m2
    psoil = out["soilp"] * tonnes_per_ha_to_g_m2
    inorgp = out["inorgp"] * tonnes_per_ha_to_g_m2
    inorglabp = out["inorglabp"] * tonnes_per_ha_to_g_m2
    inorgsorbp = out["inorgsorbp"] * tonnes_per_ha_to_g_m2
    inorgavlp = out["inorgavlp"] * tonnes_per_ha_to_g_m2
    inorgssorbp = out["inorgssorbp"] * tonnes_per_ha_to_g_m2
    inorgoccp = out["inorgoccp"] * tonnes_per_ha_to_g_m2
    inorgparp = out["inorgparp"] * tonnes_per_ha_to_g_m2
    activesoilp = out["activesoilp"] * tonnes_per_ha_to_g_m2
    slowsoilp = out["slowsoilp"] * tonnes_per_ha_to_g_m2
    passivesoilp = out["passivesoilp"] * tonnes_per_ha_to_g_m2
    ppoolo = activesoilp + slowsoilp + passivesoilp


    # fluxes outputs
    beta = out["wtfac_root"]
    nep = out["nep"] * tonnes_per_ha_to_g_m2
    gpp = out["gpp"] * tonnes_per_ha_to_g_m2
    npp = out["npp"] * tonnes_per_ha_to_g_m2
    rh = out["hetero_resp"] * tonnes_per_ha_to_g_m2
    ra = out["auto_resp"] * tonnes_per_ha_to_g_m2
    et = out["et"] # mm of water' are same value as kg/m2
    trans = out["transpiration"] # mm of water' are same value as kg/m2
    soil_evap = out["soil_evap"] # mm of water' are same value as kg/m2
    can_evap = out["canopy_evap"] # mm of water' are same value as kg/m2
    runoff = out["runoff"] # mm of water' are same value as kg/m2
    gl = out["cpleaf"] * tonnes_per_ha_to_g_m2
    # gw summed from cpstem and cpbranch below
    cpstem = out["cpstem"] * tonnes_per_ha_to_g_m2
    cpbranch = out["cpbranch"] * tonnes_per_ha_to_g_m2
    gr = out["cproot"] * tonnes_per_ha_to_g_m2
    gcr = out["cpcroot"] * tonnes_per_ha_to_g_m2
    deadleaves = out["deadleaves"] * tonnes_per_ha_to_g_m2
    deadroots = out["deadroots"] * tonnes_per_ha_to_g_m2
    deadcroots = out["deadcroots"] * tonnes_per_ha_to_g_m2
    deadbranch = out["deadbranch"] * tonnes_per_ha_to_g_m2
    deadstems = out["deadstems"] * tonnes_per_ha_to_g_m2
    deadleafn = out["deadleafn"] * tonnes_per_ha_to_g_m2
    deadbranchn = out["deadbranchn"] * tonnes_per_ha_to_g_m2
    deadstemn = out["deadstemn"] * tonnes_per_ha_to_g_m2
    deadrootn = out["deadrootn"] * tonnes_per_ha_to_g_m2
    deadcrootn = out["deadcrootn"] * tonnes_per_ha_to_g_m2
    nup = out["nuptake"] * tonnes_per_ha_to_g_m2
    ngross = out["ngross"] * tonnes_per_ha_to_g_m2
    nmin = out["nmineralisation"] * tonnes_per_ha_to_g_m2
    npleaf = out["npleaf"] * tonnes_per_ha_to_g_m2
    nproot = out["nproot"] * tonnes_per_ha_to_g_m2
    npcroot = out["npcroot"] * tonnes_per_ha_to_g_m2
    npstemimm = out["npstemimm"] * tonnes_per_ha_to_g_m2
    npstemmob = out["npstemmob"] * tonnes_per_ha_to_g_m2
    npbranch = out["npbranch"] * tonnes_per_ha_to_g_m2
    apar = out["apar"] / SW_RAD_TO_PAR
    gcd = out["gs_mol_m2_sec"]
    ga = out["ga_mol_m2_sec"]
    nleach = out["nloss"] * tonnes_per_ha_to_g_m2
    activesoil = out["activesoil"] * tonnes_per_ha_to_g_m2
    slowsoil = out["slowsoil"] * tonnes_per_ha_to_g_m2
    passivesoil = out["passivesoil"] * tonnes_per_ha_to_g_m2
    cfretransn = out["leafretransn"] * tonnes_per_ha_to_g_m2
    deadleafp = out["deadleafp"] * tonnes_per_ha_to_g_m2
    deadbranchp = out["deadbranchp"] * tonnes_per_ha_to_g_m2
    deadstemp = out["deadstemp"] * tonnes_per_ha_to_g_m2
    deadrootp = out["deadrootp"] * tonnes_per_ha_to_g_m2
    deadcrootp = out["deadcrootp"] * tonnes_per_ha_to_g_m2
    pup = out["puptake"] * tonnes_per_ha_to_g_m2
    pgross = out["pgross"] * tonnes_per_ha_to_g_m2
    pmin = out["pmineralisation"] * tonnes_per_ha_to_g_m2
    ppleaf = out["ppleaf"] * tonnes_per_ha_to_g_m2
    pproot = out["pproot"] * tonnes_per_ha_to_g_m2
    ppcroot = out["ppcroot"] * tonnes_per_ha_to_g_m2
    ppstemimm = out["ppstemimm"] * tonnes_per_ha_to_g_m2
    ppstemmob = out["ppstemmob"] * tonnes_per_ha_to_g_m2
    ppbranch = out["ppbranch"] * tonnes_per_ha_to_g_m2
    pleach = out["ploss"] * tonnes_per_ha_to_g_m2
    cfretransp = out["leafretransp"] * tonnes_per_ha_to_g_m2

    # extra traceability stuff
    tfac_soil_decomp = out["tfac_soil_decomp"]
    c_into_active = out["c_into_active"] * tonnes_per_ha_to_g_m2
    c_into_slow = out["c_into_slow"] * tonnes_per_ha_to_g_m2
    c_into_passive = out["c_into_passive"] * tonnes_per_ha_to_g_m2
    active_to_slow = out["active_to_slow"] * tonnes_per_ha_to_g_m2
    active_to_passive = out["active_to_passive"] * tonnes_per_ha_to_g_m2
    slow_to_active = out["slow_to_active"] * tonnes_per_ha_to_g_m2
    slow_to_passive = out["slow_to_passive"] * tonnes_per_ha_to_g_m2
    passive_to_active = out["passive_to_active"] * tonnes_per_ha_to_g_m2
    co2_rel_from_surf_struct_litter = out["co2_rel_from_surf_struct_litter"] * tonnes_per_ha_to_g_m2
    co2_rel_from_soil_struct_litter = out["co2_rel_from_soil_struct_litter"] * tonnes_per_ha_to_g_m2
    co2_rel_from_surf_metab_litter = out["co2_rel_from_surf_metab_litter"] * tonnes_per_ha_to_g_m2
    co2_rel_from_soil_metab_litter = out["co2_rel_from_soil_metab_litter"] * tonnes_per_ha_to_g_m2
    co2_rel_from_active_pool = out["co2_rel_from_active_pool"] * tonnes_per_ha_to_g_m2
    co2_rel_from_slow_pool = out["co2_rel_from_slow_pool"] * tonnes_per_ha_to_g_m2
    co2_rel_from_passive_pool = out["co2_rel_from_passive_pool"] * tonnes_per_ha_to_g_m2

    # extra priming stuff
    rexc = [UNDEF] * len(doy)
    rexn = [UNDEF] * len(doy)
    co2x = [UNDEF] * len(doy)
    factive = [UNDEF] * len(doy)
    rtslow = [UNDEF] * len(doy)
    rexcue = [UNDEF] * len(doy)
    cslo = out["slowsoil"] * tonnes_per_ha_to_g_m2
    nslo = out["slowsoiln"] * tonnes_per_ha_to_g_m2
    cact = out["activesoil"] * tonnes_per_ha_to_g_m2
    nact = out["activesoiln"] * tonnes_per_ha_to_g_m2


    # Misc stuff we don't output
    drainage = [UNDEF] * len(doy)
    rleaf = [UNDEF] * len(doy)
    rwood = [UNDEF] * len(doy)
    rcr = [UNDEF] * len(doy)
    rfr = [UNDEF] * len(doy)
    rgrow = [UNDEF] * len(doy)
    rsoil = [UNDEF] * len(doy)
    cex = [UNDEF] * len(doy)
    cvoc = [UNDEF] * len(doy)
    lh = [UNDEF] * len(doy)
    sh = [UNDEF] * len(doy)
    ccoarse_lit = [UNDEF] * len(doy)
    ndw = [UNDEF] * len(doy)
    pclitb = [UNDEF] * len(doy)
    nvol = [UNDEF] * len(doy)
    gb = [UNDEF] * len(doy)
    grepr = [UNDEF] * len(doy)
    cwretransn = [UNDEF] * len(doy)
    ccrretransn = [UNDEF] * len(doy)
    cfrretransn = [UNDEF] * len(doy)
    plretr = [UNDEF] * len(doy)
    pwretr = [UNDEF] * len(doy)
    pcrretr = [UNDEF] * len(doy)
    pfrretr = [UNDEF] * len(doy)

    # Misc calcs from fluxes/state
    lma = shoot / lai
    ncon = shootn / shoot
    nflit = litternag + litternbg
    pflit = litterpag + litterpbg
    pcon = shootp / shoot
    recosys = rh + ra
    secp = inorgsorbp + inorgssorbp
    cw = stem + branch
    cwp = stemp + branchp
    gw = cpstem + cpbranch
    cwn = stemn + branchn
    cwin = deadstems + deadbranch
    ccrlin = deadcroots
    cfrlin = deadroots
    ndeadwood = deadbranchn + deadstemn
    pdeadwood = deadbranchp + deadstemp
    nwood_growth = npstemimm + npstemmob + npbranch
    pwood_growth = ppstemimm + ppstemmob + ppbranch

    return {'YEAR':year, 'DOY':doy, 'SW':pawater_root, 'SWPA':pawater_root,
            'NEP':nep, 'GPP':gpp, 'NPP':npp, 'CEX':cex, 'CVOC':cvoc,
            'RECO':recosys, 'RAU':ra, 'RL':rleaf, 'RW':rwood,
            'RCR':rcr, 'RFR':rfr,
            'RGR':rgrow, 'RHET':rh, 'RSOIL':rsoil, 'ET':et, 'T':trans,
            'ES':soil_evap, 'EC':can_evap, 'RO':runoff, 'DRAIN':drainage,
            'LE':lh, 'SH':sh, 'CL':shoot, 'CW':cw, 'CCR':coarse_root,
            'CFR':fine_root, 'CSTOR':tnc, 'CFLIT':litterc, 'CFLITA':littercag,
            'CFLITB':littercbg, 'CCLITB':ccoarse_lit, 'CSOIL':soilc,
            'CGL':gl, 'CGW':gw, 'CGCR':gcr, 'CGFR':gr, 'CREPR':grepr, 'CLITIN':deadleaves,
            'CCRLIN':ccrlin, 'CFRLIN':cfrlin, 'CWLIN':cwin, 'LAI':lai, 'LMA':lma, 'NCON':ncon,
            'NL':shootn, 'NW':cwn, 'NCR':coarse_rootn, 'NFR':rootn,
            'NSTOR':nstorage, 'NFLIT': nflit, 'NFLITA':litternag, 'NFLITB':litternbg, 'NCLITB':ndw,
            'NSOIL':nsoil, 'NPMIN':inorgn, 'NPORG':npoolo, 
            'NGL':npleaf, 'NGW':nwood_growth, 'NGCR':npcroot, 'NGFR':nproot,
            'NLITIN':deadleafn, 'NCRLIN':deadcrootn,
            'NFRLIN':deadrootn, 'NWLIN':ndeadwood, 'NUP':nup,
            'NGMIN':ngross, 'NMIN':nmin, 'NVOL': nvol, 'NLEACH':nleach,
            'NLRETR':cfretransn, 'NWRETR':cwretransn,
            'NCRRETR':ccrretransn, 'NFRRETR':cfrretransn,
            'APARd':apar, 'GCd':gcd, 'GAd':ga, 'Gbd':gb, 'Betad':beta,
            'PL':shootp, 'PW':cwp,
            'PCR':crootp, 'PFR':rootp,
            'PSTOR':pstorage, 'PFLIT':pflit,
            'PFLITA':litterpag, 'PFLITB':litterpbg, 'PCLITB':pclitb,
            'PSOIL':psoil, 'PLAB':inorglabp,
            'PSEC':secp, 'POCC':inorgoccp,
            'PPAR':inorgparp,
            'PPMIN':inorgp, 'PPORG':ppoolo,
            'PLITIN':deadleafp, 'PCRLIN':deadcrootp,
            'PFRLIN':deadrootp, 'PWLIN':pdeadwood, 'PUP':pup,
            'PGMIN':pgross, 'PMIN':pmin, 'PLEACH':pleach,
            'PGL':ppleaf, 'PGW':pwood_growth, 'PGCR':ppcroot, 'PGFR':pproot,
            'PLRETR':cfretransp, 'PWRETR':pwretr, 'PFRRETR':pcrretr, 'PFRRETR':pfrretr, 
            'CTOACTIVE':c_into_active, 'CTOSLOW':c_into_slow,
            'CTOPASSIVE':c_into_passive, 'CACTIVETOSLOW':active_to_slow,
            'CACTIVETOPASSIVE':active_to_passive, 'CSLOWTOACTIVE':slow_to_active,
            'CSLOWTOPASSIVE':slow_to_passive, 'CPASSIVETOACTIVE':passive_to_active,
            'CACTIVE':activesoil, 'CSLOW':slowsoil, 'CPASSIVE':passivesoil,
            'CO2SLITSURF':co2_rel_from_surf_struct_litter,
            'CO2SLITSOIL':co2_rel_from_soil_struct_litter,
            'CO2MLITSURF':co2_rel_from_surf_metab_litter,
            'CO2MLITSOIL':co2_rel_from_soil_metab_litter,
            'CO2FSOM':co2_rel_from_active_pool,
            'CO2SSOM':co2_rel_from_slow_pool,
            'CO2PSOM':co2_rel_from_passive_pool,
            'TFACSOM':tfac_soil_decomp,
            'REXC':rexc,
            'REXN':rexn,
            'CO2X':co2x,
            'FACTIVE':factive,
            'RTSLOW':rtslow,
            'REXCUE':rexcue,
            'CSLO':cslo,
            'NSLO':nslo,
            'CACT':cact,
            'NACT':nact}, git_ver


def setup_units():
    units = ['--','--','Mean ppm', 'mm d-1', 'mol m-2', 'Mean DegC', 'Mean DegC',
                'kPa h', 'mm', 'mm', 'gN m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'kgH2O m-2 d-1',
                'kgH2O m-2 d-1', 'kgH2O m-2 d-1', 'kgH2O m-2 d-1',
                'kgH2O m-2 d-1', 'kgH2O m-2 d-1', 'MJ m-2', 'MJ m-2',
                'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2',
                'gC m-2', 'gC m-2', 'gC m-2', 'gC m-2 0 to 30 cm',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'm2 m-2', 'gC m-2',
                'gN gd.m.-1', 'gN m-2', 'gN m-2', 'gN m-2', 'gN m-2', 'gN m-2',
                'gN m-2', 'gN m-2', 'gN m-2', 'gN m-2', 'gN m-2 0 to 30 cm',
                'gN m-2 0 to 30 cm', 'gN m-2 0 to 30 cm', 'gN m-2 d-1',
                'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1',
                'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1',
                'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1',
                'gN m-2 d-1', 'gN m-2 d-1', 'gN m-2 d-1',
                'gN m-2 d-1', 'gN m-2 d-1',
                'MJ m-2 d-1', 'mol H2O m-2 s-1', 'mol H2O m-2 s-1',
                'mol H2O m-2 s-1', 'frac',
                'gP m-2', 'gP m-2', 'gP m-2',
                'gP m-2', 'gP m-2', 'gP m-2',
                'gP m-2', 'gP m-2', 'gP m-2',
                'gP m-2', 'gP m-2',
                'gP m-2', 'gP m-2',
                'gP m-2', 'gP m-2',
                'gP m-2','gP m-2 d-1', 'gP m-2 d-1',
                'gP m-2 d-1','gP m-2 d-1', 'gP m-2 d-1',
                'gP m-2 d-1', 'gP m-2 d-1', 'gP m-2 d-1',
                'gP m-2 d-1', 'gP m-2 d-1', 'gP m-2 d-1', 'gP m-2 d-1',
                'gP m-2 d-1', 'gP m-2 d-1', 'gP m-2 d-1', 'gP m-2 d-1',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2', 'gC m-2', 'gC m-2',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'gC m-2 d-1', 'gC m-2 d-1',
                'frac', 'gC m-2 d-1', 'gN m-2 d-1', 'gC m-2 d-1',
                'gC m-2 d-1', 'years', 'frac', 'gC m-2 d-1', 'gN m-2 d-1',
                'gC m-2 d-1', 'gN m-2 d-1']


    return units

def setup_varnames():
    variable = ['Year', 'Day of the year', 'CO2', 'Precipitation', 'PAR',
                'Air temp canopy', 'Soil temp 10 cm', 'Vapour Pres Def',
                'Total soil water content', 'Plant available soil water content',
                'N deposition', 'Net Eco Prod',
                'Gross Prim Prod', 'Net Prim Prod', 'C exudation',
                'C VOC Flux', 'Resp ecosystem', 'Resp autotrophic',
                'Resp leaves (maint)', 'Resp Wood (maint)',
                'Resp coarse root (maint)',
                'Resp Fine Root (maint)', 'Resp growth',
                'Resp heterotrophic',
                'Evapotranspiration', 'Transpiration', 'Soil Evaporation',
                'Canopy evaporation', 'Runoff', 'Drainage', 'Latent Energy',
                'Sensible Heat', 'C Leaf Mass', 'C Wood Mass',
                'C Coarse Root mass', 'C Fine Root mass',
                'C Storage as TNC', 'C Fine Litter Total',
                'C Fine Litter above', 'C Fine Litter below',
                'C Coarse Litter', 'C Soil', 'C Leaf growth',
                'C Wood growth', 'C Coarse Root growth',
                'C Fine Root growth', 'C reproduction growth',
                'C Leaf Litterfall',
                'C Coarse Root litter inputs', 'C Fine Root litter inputs',
                'C Wood/branch inputs',
                'LAI projected', 'Leaf gC/leaf area', 'N Conc Leaves',
                'N Mass Leaves', 'N Mass Wood', 'N Mass Coarse Roots',
                'N Mass Fine Roots', 'N storage', 'N fine litter total', 'N litter aboveground',
                'N litter belowground', 'N Dead wood', 'N Soil Total',
                'N in Mineral form', 'N in Organic form', 'N fixation',
                'N Leaf growth', 'N Wood growth', 'N CR growth', 'N Fine Root growth',
                'N Leaf Litterfall',
                'N Coarse Root litter input', 'N Fine Root litter input',
                'N Wood/brch litterfall', 'N Biomass Uptake',
                'N Gross Mineralization', 'N Net mineralization',
                'N Volatilization', 'N Leaching',
                'Foliage retranslocation',
                'Wood/Branch retranslocation', 'Coarse Root retranslocation',
                'Fine Root retranslocation',
                'Aborbed PAR', 'Average daytime canopy conductance',
                'Average daytime aerodynamic conductance',
                'Average daytime leaf boundary conductance',
                'Soil moisture stress',
                'P Mass Leaves',
                'P Mass Wood', 'P Mass Coarse Roots', 'P Mass Fine Roots',
                'P storage', 'P litter total', 'P litter aboveground', 'P litter belowground',
                'P coarse litter',
                'P Soil Total', 'P in labile form',
                'P in secondary form',
                'P in occluded form', 'P parent pool',
                'P Inorganic pool',
                'P in Organic form','P Leaf Litterfall', 
                'P Coarse Root litter input','P Fine Root litter input', 'P Wood/brch litterfall',
                'P Biomass Uptake',
                'P Gross Mineralisation', 'P Net mineralisation', 'P Leaching',
                'P Leaf growth', 'P Wood growth', 'P CR growth', 'P Fine Root growth',
                'P Foliage retranslocation',
                'P Wood/Branch retranslocation', 'P Coarse Root retranslocation',
                'P Fine Root retranslocation',
                'C fluxes from litter & slow/passive to active soil pool',
                'C fluxes from litter & active soil pool to slow pool',
                'C fluxes from active & slow soil pool to passive pool',
                'C flux from active soil pool to slow soil pool',
                'C flux from active soil pool to passive soil pool',
                'C flux from slow soil pool to active soil pool',
                'C flux from slow pool to passive soil pool',
                'C flux from passive pool to active pool',
                'C Active SOM pool',
                'C Slow SOM pool',
                'C Passive SOM pool',
                'CO2 efflux from surf structural litter',
                'CO2 efflux from soil structural litter',
                'CO2 efflux from surf metabolic litter',
                'CO2 efflux from soil metabolic litter',
                'CO2 efflux from fast SOM pool',
                'CO2 efflux from slow SOM pool',
                'CO2 efflux from passive SOM pool',
                'Temperature scalar on C efflux from SOM pools',
                'Root Exudation of C',
                'Root Exudation of N',
                'CO2 released from exudation',
                'Total C flux from the active pool',
                'Residence time of slow pool',
                'REXC carbon use efficiency',
                'Total C in the slow pool',
                'Total N in the slow pool',
                'Total C in the active pool',
                'Total N in the active pool']




    variable_names = ['YEAR', 'DOY', 'CO2', 'PREC', 'PAR', 'TAIR', 'TSOIL', 'VPD',
                      'SW', 'SWPA', 'NDEP', 'NEP', 'GPP', 'NPP', 'CEX', 'CVOC',
                      'RECO', 'RAU', 'RL', 'RW', 'RCR', 'RFR', 'RGR',
                      'RHET',
                      'ET', 'T', 'ES', 'EC', 'RO', 'DRAIN', 'LE', 'SH',
                      'CL', 'CW', 'CCR', 'CFR', 'CSTOR', 'CFLIT', 'CFLITA',
                      'CFLITB', 'CCLITB', 'CSOIL', 'CGL', 'CGW', 'CGCR', 'CGFR',
                      'CREPR','CLITIN', 'CCRLIN', 'CFRLIN','CWLIN', 'LAI',
                      'LMA', 'NCON', 'NL', 'NW', 'NCR', 'NFR', 'NSTOR',
                      'NFLIT', 'NFLITA','NFLITB', 'NCLITB', 'NSOIL', 'NPMIN', 'NPORG', 'NFIX',
                      'NGL', 'NGW', 'NGCR', 'NGFR',
                      'NLITIN', 'NCRLIN', 'NFRLIN','NWLIN', 'NUP', 'NGMIN', 'NMIN',
                      'NVOL', 'NLEACH', 'NLRETR', 'NWRETR',
                      'NCRRETR', 'NFRRETR', 'APARd',
                      'GCd', 'GAd', 'GBd', 'Betad',
                      'PL', 'PW',
                      'PCR', 'PFR','PSTOR',
                      'PFLIT', 'PFLITA', 'PFLITB', 'PCLITB',
                      'PSOIL', 'PLAB',
                      'PSEC', 'POCC',
                      'PPAR',
                      'PPMIN', 'PPORG',
                      'PLITIN', 'PCRLIN',
                      'PFRLIN', 'PWLIN', 'PUP',
                      'PGMIN', 'PMIN', 'PLEACH',
                      'PGL', 'PGW', 'PGCR', 'PGFR',
                      'PLRETR', 'PWRETR', 'PCRRETR', 'PFRRETR',
                      'CTOACTIVE', 'CTOSLOW', 'CTOPASSIVE', 'CACTIVETOSLOW',
                      'CACTIVETOPASSIVE', 'CSLOWTOACTIVE', 'CSLOWTOPASSIVE',
                      'CPASSIVETOACTIVE', 'CACTIVE', 'CSLOW', 'CPASSIVE',
                      'CO2SLITSURF', 'CO2SLITSOIL', 'CO2MLITSURF',
                      'CO2MLITSOIL', 'CO2FSOM', 'CO2SSOM', 'CO2PSOM',
                      'TFACSOM','REXC','REXN','CO2X','FACTIVE','RTSLOW','REXCUE',
                      'CSLO','NSLO','CACT','NACT']



    return variable, variable_names

