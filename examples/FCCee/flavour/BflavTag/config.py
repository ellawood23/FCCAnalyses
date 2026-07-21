# config.py
# This file contains all of the important configuration options used throughtout the analysis
# This file, along with all other analysis scripts, are expected to be in the `FCCAnalysesPath` directory
# Also contains branching fractions, cut efficiencies, etc which are manually filled in for now
import os
import matplotlib.pyplot as plt
import numpy as np

# MANDATORY ----> replace the default string with the path to the working directory in the FCCAnalyses repo
FCCAnalysesPath = "/usera/ejnw2/PhD/FCC_FT/FCCAnalyses/examples/FCCee/flavour/BflavTag/"
FCCAnalysesPath = os.path.abspath(FCCAnalysesPath)
FT_outputDir = "/r02/lhcb/ejnw2/FCC_FT/FCC_FT_outputs_full_July2026"


# RUNNING MODE
run_mode_choices = ['BflavTag_full'] #when add run mode, now need to add to processList, fccana_opts AND PROCESS_TUPLES!!!


run_mode = 'BflavTag_full'
if run_mode not in run_mode_choices:
    raise RuntimeError(f'{run_mode} is not a valid run mode')


##############################
## CONFIG DICTS
##############################
# processList to pass to `fccanalysis run`
processList = {
    # Size of winter2023 samples in /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/:
    # p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu == 13G (2,000,000 events)
    # p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu ~= 13G (2,200,000 events)
  

    "BflavTag_full":{"p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu":{"fraction": 1, "chunks": 2},
                "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu":{"fraction": 1, "chunks": 2},},


}

# Default options to pass to `fccanalysis run`
fccana_opts = {
    "prodTag":   "FCCee/winter2023/IDEA",
    "outputDir": {
        "BflavTag_full": os.path.abspath(FT_outputDir),
    },

    "testFile": {
        "Bs": "root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu/events_026683563.root",
        "Bd": "root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu/events_004838962.root",
    },
    "analysisName":   "BflavTag",
    "nCPUS":          8,
    "runBatch":       True,
    "batchQueue":     "workday",
    "compGroup":      "group_u_FCC.local_gen",
    "yamlPath":       os.path.join(FCCAnalysesPath, "BflavTag_branches.yaml"),  # Path to the YAML file containing feature names
    "outputBranches": {
        "BflavTag_full":"fcc-FT-vars",
    },
}



#BSC taken from https://github.com/HEP-FCC/FCCeePhysicsPerformance/blob/master/General/README.md#generating-events-under-realistic-fcc-ee-environment-conditions and agreement checked with MC samples
#nb. if spring2021 values used for winter2023, BSC is too tight --> error and slow fitting: `VertexFit::RegInv: null determinant for N = 2`
BSC_opts = {
    "winter2023": [5.96,23.8e-3,0.397e3], # vertex sigma [x,y,z] in micrometers
    "spring2021": [4.5,20e-3,0.3e3],
}

quark_dictionary = {"b": 5,
                    "c":4,
                    "s":3,
                    "d":1,
                    "u":2,}



##############################
## SAMPLE OPTIONS
##############################

samples = ["p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu",
    "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu"]



sample_allocations = {
    "Bssignal":     ["p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu"],
    "Bdsignal":   ["p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu"],
    "combined_signal": ["p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu", "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu"],#want the signal last here so that it's plotted on top in binning plots
}


sample_colors = {
    "combined_signal": plt.cm.Blues( np.linspace(0, 1, 6)[3:-1] ),
    "Bssignal":  plt.cm.Blues( np.linspace(0, 1, 6)[3] ) ,
    "Bdsignal":  plt.cm.Blues( np.linspace(0, 1, 6)[-2] ) ,
}

sample_total = {
    "combined_signal": None,
}

sample_hatches = {
    "combined_signal": [r'////', r'\\\\'],
    "Bdsignal": [r'////'],
    "Bssignal": [r'\\\\'],
}

sample_shorthand = {
    "p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu": "Bs2NuNu",
    "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu": "Bd2NuNu",
}

titles = {
    "p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu": r"$B_s^0 \to \nu \bar{\nu}$",
    "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu": r"$B^0 \to \nu \bar{\nu}$",
    "p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu_invis": r"$B_s^0 \to$ invisible",
    "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu_invis": r"$B^0 \to$ invisible",
    "signal": r"$B_{(s)}^0 \to \nu \bar{\nu}$",
}

##############################
## NUMERICAL DATA
##############################
# from PDG -> B production fractions
LEP_prod_fracs = { # taken from https://hflav-eos.web.cern.ch/hflav-eos/osc/PDG_2021/#FRAC 
    "Bu":  (0.408, 0.007),# 0.43 - old #s from https://indico.in2p3.fr/event/23012/contributions/89940/attachments/61988/84706/Hill-fcc-france.pdf
    "Bd":  (0.408, 0.007),#0.43,
    "Bs": (0.100, 0.008),#0.096,
    "Bc": (0.0004,0.0004),#Not measured by LEP or at Z pole to date therefore add 100% error, https://arxiv.org/pdf/hep-ph/9707248
    #"Lb": #0.037, 
    }

prod_frac = {
    # Actually, use the sample name to make integration with Bd2NuNu easier
    "p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu": LEP_prod_fracs["Bs"],
    "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu": LEP_prod_fracs["Bd"],
}


# from PDG (value, error)
# Z branching fractions
# Z->ss = (Z->dd+ss+bb)/3 * 3 - Z->bb / 2
# Z->uu/cc = 2 * (11.6 +/- 0.6) = 23.2 +/- 1.2
# Z->dd/ss/bb = 3 * (15.6 +/- 0.4) = 46.8 +/- 1.2
#Z->τ+τ− = (3.3696±0.0083) %
#Z->mu+mu− = (3.3662±0.0066) %
#Z->e+e− = (3.3632±0.0042) %

branching_fractions = {
    "p8_ee_Zbb_ecm91_EvtGen_Bs2NuNu": (1, 0),  # a dummy value
    "p8_ee_Zbb_ecm91_EvtGen_Bd2NuNu": (1, 0),  # a dummy value
    "p8_ee_Zbb_ecm91": (0.1512, 0.0005),
}

mass_Z = 91.188  # Ecm used in the winter2023 samples

N_z = 6e12 # total number of Nz expected across all experiments during tera-Z run (from https://arxiv.org/pdf/2309.11353 Matt/Aidan paper)

B_z = 2.0 #magnetic field strength from delphes card: https://raw.githubusercontent.com/HEP-FCC/FCC-config/winter2023/FCCee/Delphes/card_IDEA.tcl

