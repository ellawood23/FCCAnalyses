########## INSTRUCTIONS ############
#
# Use this file to process B -> inv tuples for flavour tagging
# Runs configurations
#
#   -BflavTag_full:
#       Processes full FT B2inv samples, with only cut to ensure PV present
#
####################################


import os
import sys

# Config and yaml file must be in this directory by default
# Absolute path must be supplied for the script to work in batch mode
configPath = '/usera/ejnw2/PhD/FCC_FT/FCCAnalyses/examples/FCCee/flavour/BflavTag/'
sys.path.append(os.path.abspath(configPath))


import ROOT
from yaml import safe_load
from math import sqrt

import config as cfg

#Mandatory: List of processes
processList = cfg.processList[cfg.run_mode]

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag = cfg.fccana_opts['prodTag']

#Optional: output directory, default is local running directory
outputDir = cfg.fccana_opts['outputDir'][cfg.run_mode]
#outputDirEos = cfg.fccana_opts['outputDirEos'][cfg.run_mode]
#eosType = "eospublic"

#Optional: analysisName, default is ""
analysisName = cfg.fccana_opts['analysisName']

#Optional: ncpus, default is 4
nCPUS = cfg.fccana_opts['nCPUS']

#Optional running on HTCondor, default is False
runBatch = cfg.fccana_opts['runBatch']

#Optional test file
testFile = cfg.fccana_opts['testFile']['Bd']

print("----> INFO: Using config.py file from:")
print(f"{15*' '}{os.path.abspath(configPath)}")
print("----> INFO: Using branch names from:")
print(f"{15*' '}{cfg.fccana_opts['yamlPath']}")


class RDFanalysis():

    #__________________________________________________________
    def analysers(df):
        
        #BSC for vertexing
        #bsc = [ 6, 25e-3, 400 ]
        bsc = cfg.BSC_opts['winter2023'] # list of sigmax,sigmay,sigmaz
        
        """## name of collections in EDM root files
        collections = {
            "GenParticles": "Particle",
            "PFParticles": "ReconstructedParticles",
            "PFTracks": "EFlowTrack",
            "PFPhotons": "EFlowPhoton",
            "PFNeutralHadrons": "EFlowNeutralHadron",
            "TrackState": "EFlowTrack_1",
            "TrackerHits": "TrackerHits",
            "CalorimeterHits": "CalorimeterHits",
            "dNdx": "EFlowTrack_2",
            "PathLength": "EFlowTrack_L",
            "Bz": "magFieldBz",
        }"""

        df2 = (
            df
            #############################################
            ##          Aliases for # in python        ##
            #############################################
            .Alias("MCRecoAssociationsRec", "MCRecoAssociations#0.index")  # points to ReconstructedParticles
            .Alias("MCRecoAssociationsGen", "MCRecoAssociations#1.index")  # points to Particle
            .Alias("ParticleParents",       "Particle#0.index")            # gen particle parents
            .Alias("ParticleChildren",      "Particle#1.index")            # gen particle children

            #All MC particles needed for is to make event displays
            #############################################
            ##             MC IDs and Status           ##
            #############################################
            .Define("MC_n",           "MCParticle::get_n(Particle)")
            .Define("MC_genStatus",   "MCParticle::get_genStatus(Particle)")
            .Define("MC_PDG",         "MCParticle::get_pdg(Particle)")  # this true ID from MC
            .Define("MC_M1",          "myUtils::getMC_parent(0,Particle,ParticleParents)")
            .Define("MC_M2",          "myUtils::getMC_parent(1,Particle,ParticleParents)")
            .Define("MC_D1",          "myUtils::getMC_daughter(0,Particle,ParticleChildren)")
            .Define("MC_D2",          "myUtils::getMC_daughter(1,Particle,ParticleChildren)")
            .Define("MC_D3",          "myUtils::getMC_daughter(2,Particle,ParticleChildren)")
            .Define("MC_D4",          "myUtils::getMC_daughter(3,Particle,ParticleChildren)")


            #############################################
            ##               MC Particles              ##
            #############################################
            .Define("MC_e",           "MCParticle::get_e(Particle)")
            .Define("MC_m",           "MCParticle::get_mass(Particle)")
            .Define("MC_q",           "MCParticle::get_charge(Particle)")
            .Define("MC_p",           "MCParticle::get_p(Particle)")
            .Define("MC_pt",          "MCParticle::get_pt(Particle)")
            .Define("MC_px",          "MCParticle::get_px(Particle)")
            .Define("MC_py",          "MCParticle::get_py(Particle)")
            .Define("MC_pz",          "MCParticle::get_pz(Particle)")
            .Define("MC_eta",         "MCParticle::get_eta(Particle)")
            .Define("MC_phi",         "MCParticle::get_phi(Particle)")
            .Define("MC_orivtx_x",     "MCParticle::get_vertex_x(Particle)")
            .Define("MC_orivtx_y",     "MCParticle::get_vertex_y(Particle)")
            .Define("MC_orivtx_z",     "MCParticle::get_vertex_z(Particle)")
            .Define("MC_endPoint_x",    "MCParticle::get_endPoint_x(Particle, ParticleChildren)") #warning: if a neutral B has oscillated into a Bbar, this returns the production vertex of the Bbar
            .Define("MC_endPoint_y",    "MCParticle::get_endPoint_y(Particle, ParticleChildren)")  #warning: if a neutral B has oscillated into a Bbar, this returns the production vertex of the Bbar
            .Define("MC_endPoint_z",    "MCParticle::get_endPoint_z(Particle, ParticleChildren)")   #warning: if a neutral B has oscillated into a Bbar, this returns the production vertex of the Bbar


            ##################################################
            ## MC variables to help with understanding event## - nb. these cause issues for taus
            ##################################################
            # Pythia8 generatorStatus
            # 21 - incoming particles of hardest process (e+ e- beams)
            # 22 - intermediate particles of hardest process (Z)
            # 23 - outgoing particles of hardest process (quark pair produced from Z)
            #  1 - final-state particles
            .Define("MC_ee",          "MCParticle::sel_genStatus(21)(Particle)")   # INTERMEDIATE
            .Define("MC_Z",           "MCParticle::sel_genStatus(22)(Particle)")   # INTERMEDIATE
            .Define("MC_qq",          "MCParticle::sel_genStatus(23)(Particle)")   # INTERMEDIATE
            .Define("MC_FS",          "MCParticle::sel_genStatus(1)(Particle)") # INTERMEDIATE
        
            
            # --------------------------------------- #
            #    MC final-state particle variables    #
            # --------------------------------------- #
            .Define("MCfinal_PDG",       "MCParticle::get_pdg(MC_FS)")
            .Define("MCfinal_e",         "MCParticle::get_e(MC_FS)")
            .Define("MCfinal_m",         "MCParticle::get_mass(MC_FS)")
            .Define("MCfinal_q",         "MCParticle::get_charge(MC_FS)")
            .Define("MCfinal_p",         "MCParticle::get_p(MC_FS)")
            .Define("MCfinal_pt",        "MCParticle::get_pt(MC_FS)")
            .Define("MCfinal_px",        "MCParticle::get_px(MC_FS)")
            .Define("MCfinal_py",        "MCParticle::get_py(MC_FS)")
            .Define("MCfinal_pz",        "MCParticle::get_pz(MC_FS)")
            .Define("MCfinal_eta",       "MCParticle::get_eta(MC_FS)")
            .Define("MCfinal_phi",       "MCParticle::get_phi(MC_FS)")
            .Define("MCfinal_orivtx_x",  "MCParticle::get_vertex_x(MC_FS)")
            .Define("MCfinal_orivtx_y",  "MCParticle::get_vertex_y(MC_FS)")
            .Define("MCfinal_orivtx_z",  "MCParticle::get_vertex_z(MC_FS)")


            #############################################
            ##            MC PrimaryVertex             ##
            #############################################

            .Define("MC_PrimaryVertex",  "MCParticle::get_EventPrimaryVertex(21)(Particle)")

            .Define("MC_PV_x",  "MC_PrimaryVertex.X()") 
            .Define("MC_PV_y",  "MC_PrimaryVertex.Y()") 
            .Define("MC_PV_z",  "MC_PrimaryVertex.Z()")


            #############################################
            ##           Find MC Vertices              ##
            #############################################

            .Define("MC_VertexObject",   "myUtils::get_MCVertexObject(Particle, ParticleParents)")
            
            .Define("MC_vtx_n",        "int(MC_VertexObject.size())")
            .Define("MC_vtx_ntracks",  "myUtils::get_NTracksMCVertex(MC_VertexObject)")
            .Define("MC_vtx_indMC",    "myUtils::get_MCindMCVertex(MC_VertexObject)")
            .Define("MC_vtx_x",        "myUtils::get_MCVertex_x(MC_VertexObject)")
            .Define("MC_vtx_y",        "myUtils::get_MCVertex_y(MC_VertexObject)")
            .Define("MC_vtx_z",        "myUtils::get_MCVertex_z(MC_VertexObject)")
            
            # MCParticle variable that needed MC_VertexObject
            .Define("MC_orivtx_ind",  "myUtils::get_MCVertex_fromMC(Particle, MC_VertexObject)")
            


            ###############################################
            ##      Perform vertex fitting - seeded from MC (initial attempt)         
            ###############################################

            # get all MC vertices
            #.Define("MC_VertexObject",          "myUtils::get_MCVertexObject(Particle, ParticleParents)") #defined earlier
            # use this to seed the Rec vertexing
            .Define("Rec_VertexObject",        f"myUtils::get_VertexObject(MC_VertexObject, ReconstructedParticles, EFlowTrack_1, MCRecoAssociationsRec, MCRecoAssociationsGen, {bsc[0]}, {bsc[1]}, {bsc[2]})")
                        
            ## Define seeded vertex variables 
            # Filter events with no seeded PV. Otherwise Snapshot would throw a strop
            .Define("EVT_hasPV",                "myUtils::hasPV(Rec_VertexObject)")
            .Filter("EVT_hasPV==1")

            # add the PID hypothesis info to the RecParticles (based on MC truth - ie assume perfect PID here)
            .Define("RecoParticlesPID",          "myUtils::PID(ReconstructedParticles, MCRecoAssociationsRec, MCRecoAssociationsGen, Particle)")
            # now update reco momentum based on the rec vertex position
            .Define("RecoParticlesPIDAtVertex",  "myUtils::get_RP_atVertex(RecoParticlesPID, Rec_VertexObject)")


            ##############################################
            # Now also make variables where momentum updated from vtx fit
            ##############################################
            
            #-----------------------------------------------------------------------------------
            #Rec particle variables with momentum updated from seeded vertex fit
            #-----------------------------------------------------------------------------------

            .Define("Rec_n",         "ReconstructedParticle::get_n(RecoParticlesPIDAtVertex)")
            .Define("Rec_type",      "ReconstructedParticle::get_type(RecoParticlesPIDAtVertex)") # this PID from Delphes, ie. straight from particle flow. Therefore for charged particles usually same as MC, but for neutral hadrons, all KL
            .Define("Rec_m",         "ReconstructedParticle::get_mass(RecoParticlesPIDAtVertex)")
            .Define("Rec_q",         "ReconstructedParticle::get_charge(RecoParticlesPIDAtVertex)")
            .Define("Rec_eta",       "ReconstructedParticle::get_eta(RecoParticlesPIDAtVertex)")
            .Define("Rec_phi",       "ReconstructedParticle::get_phi(RecoParticlesPIDAtVertex)")
           
            .Define("Rec_indvtx",    "myUtils::get_Vertex_fromRP(RecoParticlesPIDAtVertex, Rec_VertexObject)")
            .Define("Rec_e",         "ReconstructedParticle::get_e(RecoParticlesPIDAtVertex)")
            .Define("Rec_p",         "ReconstructedParticle::get_p(RecoParticlesPIDAtVertex)")
            .Define("Rec_pt",        "ReconstructedParticle::get_pt(RecoParticlesPIDAtVertex)")
            .Define("Rec_px",        "ReconstructedParticle::get_px(RecoParticlesPIDAtVertex)")
            .Define("Rec_py",        "ReconstructedParticle::get_py(RecoParticlesPIDAtVertex)")
            .Define("Rec_pz",        "ReconstructedParticle::get_pz(RecoParticlesPIDAtVertex)")
           

            #--------------------------------------------------------
            # All Rec Vertices (seeded)
            #--------------------------------------------------------

            .Define("Rec_vtx_n",               "float(Rec_VertexObject.size())")
            .Define("Rec_vtx_indRP",           "myUtils::get_Vertex_ind(Rec_VertexObject)")
            .Define("Rec_vtx_chi2",            "myUtils::get_Vertex_chi2(Rec_VertexObject)")
            .Define("Rec_vtx_isPV",            "myUtils::get_Vertex_isPV(Rec_VertexObject)")
            .Define("Rec_vtx_ntracks",         "myUtils::get_Vertex_ntracks(Rec_VertexObject)")
            .Define("Rec_vtx_m",               "myUtils::get_Vertex_mass(Rec_VertexObject, RecoParticlesPIDAtVertex)")
            .Define("Rec_vtx_x",               "myUtils::get_Vertex_x(Rec_VertexObject)")
            .Define("Rec_vtx_y",               "myUtils::get_Vertex_y(Rec_VertexObject)")
            .Define("Rec_vtx_z",               "myUtils::get_Vertex_z(Rec_VertexObject)")
            .Define("Rec_vtx_xerr",            "myUtils::get_Vertex_xErr(Rec_VertexObject)")
            .Define("Rec_vtx_yerr",            "myUtils::get_Vertex_yErr(Rec_VertexObject)")
            .Define("Rec_vtx_zerr",            "myUtils::get_Vertex_zErr(Rec_VertexObject)")

            .Define("Rec_vtx_d2PV",            "myUtils::get_Vertex_d2PV(Rec_VertexObject,-1)")   # INTERMEDIATE
            .Define("Rec_vtx_d2PV_x",          "myUtils::get_Vertex_d2PV(Rec_VertexObject, 0)")
            .Define("Rec_vtx_d2PV_y",          "myUtils::get_Vertex_d2PV(Rec_VertexObject, 1)")
            .Define("Rec_vtx_d2PV_z",          "myUtils::get_Vertex_d2PV(Rec_VertexObject, 2)")
            .Define("Rec_vtx_d2PV_mag",        "sqrt(Rec_vtx_d2PV_x*Rec_vtx_d2PV_x+Rec_vtx_d2PV_y*Rec_vtx_d2PV_y+Rec_vtx_d2PV_z*Rec_vtx_d2PV_z)") 
            .Define("Rec_vtx_d2PV_err",        "myUtils::get_Vertex_d2PVError(Rec_VertexObject,-1)")
            .Define("Rec_vtx_d2PV_xerr",       "myUtils::get_Vertex_d2PVError(Rec_VertexObject, 0)")
            .Define("Rec_vtx_d2PV_yerr",       "myUtils::get_Vertex_d2PVError(Rec_VertexObject, 1)")
            .Define("Rec_vtx_d2PV_zerr",       "myUtils::get_Vertex_d2PVError(Rec_VertexObject, 2)")
            .Define("Rec_vtx_normd2PV",        "Rec_vtx_d2PV / Rec_vtx_d2PV_err")   # INTERMEDIATE
            .Define("Rec_vtx_normd2PV_x",      "Rec_vtx_d2PV_x / Rec_vtx_d2PV_xerr")
            .Define("Rec_vtx_normd2PV_y",      "Rec_vtx_d2PV_y / Rec_vtx_d2PV_yerr")
            .Define("Rec_vtx_normd2PV_z",      "Rec_vtx_d2PV_z / Rec_vtx_d2PV_zerr")


            #--------------------------------------------------------
            # EVT variables from seeded vertices
            #--------------------------------------------------------

            ## Construct the Thrust Axis       
            .Define("EVT_ThrustInfoNoPointing",     'Algorithms::minimize_thrust("Minuit2","Migrad")(Rec_px, Rec_py, Rec_pz)') 
            .Define("EVT_ThrustCosThetaNoPointing", "Algorithms::getAxisCosTheta(EVT_ThrustInfoNoPointing, Rec_px, Rec_py, Rec_pz)")
            .Define("EVT_ThrustInfo",               "Algorithms::getThrustPointing(1.)(EVT_ThrustCosThetaNoPointing, Rec_e, EVT_ThrustInfoNoPointing)")
            .Define("Rec_thrustCosTheta",           "Algorithms::getAxisCosTheta(EVT_ThrustInfo, Rec_px, Rec_py, Rec_pz)")
            .Define("Rec_in_hemisEmin",             "myUtils::get_RP_inHemis(1)(Rec_thrustCosTheta)")
            .Define("Rec_in_hemisEmax",             "myUtils::get_RP_inHemis(0)(Rec_thrustCosTheta)")

            .Define("EVT_Thrust_mag",          "EVT_ThrustInfo.at(0)")
            .Define("EVT_unitThrust_x",            "myUtils::norm_RVec_x(EVT_ThrustInfo.at(1),EVT_ThrustInfo.at(3),EVT_ThrustInfo.at(5))")
            .Define("EVT_unitThrust_y",            "myUtils::norm_RVec_x(EVT_ThrustInfo.at(3),EVT_ThrustInfo.at(1),EVT_ThrustInfo.at(5))")
            .Define("EVT_unitThrust_z",            "myUtils::norm_RVec_x(EVT_ThrustInfo.at(5),EVT_ThrustInfo.at(3),EVT_ThrustInfo.at(1))")

            #Sum all particle momenta 
            .Define("EVT_sum_Rec_px",  "ROOT::VecOps::Sum(Rec_px)")
            .Define("EVT_sum_Rec_py",  "ROOT::VecOps::Sum(Rec_py)")
            .Define("EVT_sum_Rec_pz",  "ROOT::VecOps::Sum(Rec_pz)")

            .Define("EVT_p", "sqrt(EVT_sum_Rec_px*EVT_sum_Rec_px+EVT_sum_Rec_py*EVT_sum_Rec_py+EVT_sum_Rec_pz*EVT_sum_Rec_pz)") 
            .Define("EVT_e", "ROOT::VecOps::Sum(Rec_e)") 



            ###############################################
            ##    Perform vertex fitting - NOT SEEDED from MC        
            ###############################################

            # Get collection of tracks consistent with a PV (i.e. not downstream Ks, Lb etc. tracks)
            # using the get_PrimaryTracks() method with a beam spot constraint under the following parameters
            # bsc_sigma(x,y,z) = (6, 25e-3, 400)
            # bsc_(x,y,z) = (0,0,0)
            
            # First the PV - select tracks reconstructed as primaries
            .Define("Rec_PrimaryTracks",       f"VertexFitterSimple::get_PrimaryTracks( EFlowTrack_1, true, {bsc[0]}, {bsc[1]}, {bsc[2]}, 0., 0., 0.)")
            .Define("Rec_n_primary_tracks",     "ReconstructedParticle2Track::getTK_n( Rec_PrimaryTracks )")
            # Then fit the PV using these tracks
            .Define("Rec_PrimaryVertexObject", f"VertexFitterSimple::VertexFitter_Tk( 1, Rec_PrimaryTracks, true, {bsc[0]}, {bsc[1]}, {bsc[2]} )")
            .Define("Rec_PrimaryVertex",        "Rec_PrimaryVertexObject.vertex")
            
            # Get secondary tracks
            .Define("Rec_SecondaryTracks",      "VertexFitterSimple::get_NonPrimaryTracks( EFlowTrack_1, Rec_PrimaryTracks )")
            .Define("Rec_n_secondary_tracks",   "ReconstructedParticle2Track::getTK_n( Rec_SecondaryTracks )")
            
            #Now add SV fit too using LCFIPlus (without jet clustering)
            #letting V0_rej, chi2_cut, invM_cut, chi2Tr_cut all resort to default values for now: V0_rej = True, chi2_cut= 9, invM_cut=10, chi2Tr_cut=5
            .Define("Rec_SecondaryVertexObject",  "VertexFinderLCFIPlus::get_SV_event(Rec_SecondaryTracks, EFlowTrack_1,  Rec_PrimaryVertexObject)")
            
            # Add back in V0 vertex objects (removed with V0_rej = True in LCFIPlus) - tight true to match what removed from SVs, chi2_cut left as default ie. 9 
            # V0 vertices essentially KS), Lambda0 and photon conversion
            .Define("Rec_V0VertexV0Object",  "VertexFinderLCFIPlus::get_V0s(Rec_SecondaryTracks, Rec_PrimaryVertexObject, true)") #FCCAnalysesV0 type
            .Define("Rec_V0VertexObject",  "Rec_V0VertexV0Object.vtx") #FCCAnalysesVertex type


            #---------------------------------
            #Additional PV fits without BSC for comparion
            #------------------------------------
            # First the PV - select tracks reconstructed as primaries
            #.Define("Rec_PrimaryTracks_noBSC",       f"VertexFitterSimple::get_PrimaryTracks( EFlowTrack_1, false)")
            #.Define("Rec_n_primary_tracks_noBSC",     "ReconstructedParticle2Track::getTK_n( Rec_PrimaryTracks_noBSC )")
            # Then fit the PV using these tracks
            #.Define("Rec_PrimaryVertexObject_noBSC", f"VertexFitterSimple::VertexFitter_Tk( 1, Rec_PrimaryTracks_noBSC, false)")
            #.Define("Rec_PrimaryVertex_noBSC",        "Rec_PrimaryVertexObject_noBSC.vertex")

            #.Define("Rec_PV_noBSC_ntracks",  "float(Rec_PrimaryTracks_noBSC.size())")
            #.Define("Rec_PV_noBSC_x",        "Rec_PrimaryVertex_noBSC.position.x")
            #.Define("Rec_PV_noBSC_y",        "Rec_PrimaryVertex_noBSC.position.y")
            #.Define("Rec_PV_noBSC_z",        "Rec_PrimaryVertex_noBSC.position.z")


            # BSC fit with tracks found without BSC
            #.Define("Rec_PrimaryVertexObject_noBSCtrkfind", f"VertexFitterSimple::VertexFitter_Tk( 1, Rec_PrimaryTracks_noBSC, true, {bsc[0]}, {bsc[1]}, {bsc[2]} )")
            #.Define("Rec_PrimaryVertex_noBSCtrkfind",        "Rec_PrimaryVertexObject_noBSCtrkfind.vertex")

            #.Define("Rec_PV_noBSCtrkfind_ntracks",  "float(Rec_PrimaryTracks_noBSC.size())")
            #.Define("Rec_PV_noBSCtrkfind_x",        "Rec_PrimaryVertex_noBSCtrkfind.position.x")
            #.Define("Rec_PV_noBSCtrkfind_y",        "Rec_PrimaryVertex_noBSCtrkfind.position.y")
            #.Define("Rec_PV_noBSCtrkfind_z",        "Rec_PrimaryVertex_noBSCtrkfind.position.z")

            #fit tracks found with BSC without BSC
            #.Define("Rec_PrimaryVertexObject_noBSCfit", f"VertexFitterSimple::VertexFitter_Tk( 1, Rec_PrimaryTracks, false)")
            #.Define("Rec_PrimaryVertex_noBSCfit",        "Rec_PrimaryVertexObject_noBSCfit.vertex")

            #.Define("Rec_PV_noBSCfit_ntracks",  "float(Rec_PrimaryTracks.size())")
            #.Define("Rec_PV_noBSCfit_x",        "Rec_PrimaryVertex_noBSCfit.position.x")
            #.Define("Rec_PV_noBSCfit_y",        "Rec_PrimaryVertex_noBSCfit.position.y")
            #.Define("Rec_PV_noBSCfit_z",        "Rec_PrimaryVertex_noBSCfit.position.z")
           
           
           
           
            #############################################
            ##       Define non-seeded vertex fit variables  
            #############################################

            # PV (non-seeded)
            .Define("Rec_PV_ntracks",  "float(Rec_PrimaryTracks.size())")
            .Define("Rec_PV_x",        "Rec_PrimaryVertex.position.x")
            .Define("Rec_PV_y",        "Rec_PrimaryVertex.position.y")
            .Define("Rec_PV_z",        "Rec_PrimaryVertex.position.z")
            .Define("Rec_PV_xerr",        "sqrt(Rec_PrimaryVertex.covMatrix[0])")
            .Define("Rec_PV_yerr",         "sqrt(Rec_PrimaryVertex.covMatrix[2])")
            .Define("Rec_PV_zerr",         "sqrt(Rec_PrimaryVertex.covMatrix[5])")
            .Define("Rec_PV_chi2",     "Rec_PrimaryVertex.chi2") #chi2 of PV fit - used to check for if PV actually fitter

            ## SV info 
            .Define("EVT_nSV",               "VertexingUtils::get_n_SV(Rec_SecondaryVertexObject)")
            .Define("Rec_SV_ntracks",               "VertexingUtils::get_VertexNtrk(Rec_SecondaryVertexObject)")
            .Define("Rec_SV_chi2",            "VertexingUtils::get_chi2_SV(Rec_SecondaryVertexObject)") # SV chi2 (unnormalised) can also have normalised
            .Define("Rec_SV_nDOF",            "VertexingUtils::get_nDOF_SV(Rec_SecondaryVertexObject)") 
            .Define("Rec_SV_m",               "VertexingUtils::get_invM(Rec_SecondaryVertexObject)") # invariant mass of a vertex (assuming all tracks to be pions)
            .Define("Rec_SV_p",               "VertexingUtils::get_pMag_SV(Rec_SecondaryVertexObject)")
            .Define("Rec_SV_position",               "VertexingUtils::get_position_SV(Rec_SecondaryVertexObject)") # now need to turn into x,y,z
            .Define("Rec_SV_x",              "ROOT::VecOps::Map(Rec_SV_position, std::mem_fn(&TVector3::X))")
            .Define("Rec_SV_y",              "ROOT::VecOps::Map(Rec_SV_position, std::mem_fn(&TVector3::Y))")
            .Define("Rec_SV_z",              "ROOT::VecOps::Map(Rec_SV_position, std::mem_fn(&TVector3::Z))")
            .Define("Rec_SV_xerr",            "myUtils::get_Vertex_xErr(Rec_SecondaryVertexObject)")
            .Define("Rec_SV_yerr",            "myUtils::get_Vertex_yErr(Rec_SecondaryVertexObject)")
            .Define("Rec_SV_zerr",            "myUtils::get_Vertex_zErr(Rec_SecondaryVertexObject)")
            .Define("Rec_SV_d2PV",            "VertexingUtils::get_d3d_SV(Rec_SecondaryVertexObject, Rec_PrimaryVertexObject)") # magnitude of vector of distances of all reconstructed SV from PV (in mm in 3D)
            .Define("Rec_SV_d2PV_xy",          "VertexingUtils::get_dxy_SV(Rec_SecondaryVertexObject, Rec_PrimaryVertexObject)") #perpendicular projection of vector of distances of all reconstructed SV from PV (in mm in xy plane)

            ##V0 info
            .Define("EVT_nV0",               "VertexingUtils::get_n_SV(Rec_V0VertexObject)")
            .Define("Rec_V0_chi2",            "VertexingUtils::get_chi2_SV(Rec_V0VertexObject)") # SV chi2 (unnormalised) can also have normalised
            .Define("Rec_V0_nDOF",            "VertexingUtils::get_nDOF_SV(Rec_V0VertexObject)") 
            .Define("Rec_V0_p",               "VertexingUtils::get_pMag_SV(Rec_V0VertexObject)")
            .Define("Rec_V0_position",          "VertexingUtils::get_position_SV(Rec_V0VertexObject)") # now need to turn into x,y,z
            .Define("Rec_V0_x",               "ROOT::VecOps::Map(Rec_V0_position, std::mem_fn(&TVector3::X))")
            .Define("Rec_V0_y",               "ROOT::VecOps::Map(Rec_V0_position, std::mem_fn(&TVector3::Y))")
            .Define("Rec_V0_z",               "ROOT::VecOps::Map(Rec_V0_position, std::mem_fn(&TVector3::Z))")
            .Define("Rec_V0_xerr",            "myUtils::get_Vertex_xErr(Rec_V0VertexObject)")
            .Define("Rec_V0_yerr",            "myUtils::get_Vertex_yErr(Rec_V0VertexObject)")
            .Define("Rec_V0_zerr",            "myUtils::get_Vertex_zErr(Rec_V0VertexObject)")
            .Define("Rec_V0_d2PV",            "VertexingUtils::get_d3d_SV(Rec_V0VertexObject, Rec_PrimaryVertexObject)") 
            .Define("Rec_V0_d2PV_xy",          "VertexingUtils::get_dxy_SV(Rec_V0VertexObject, Rec_PrimaryVertexObject)") 
            .Define("Rec_V0_type",                 "VertexingUtils::get_pdg_V0(Rec_V0VertexV0Object)") #vector of V0 ID from reconstruction (ie. KS, lambda0 or photon conversion)
            .Define("Rec_V0_m",                 "VertexingUtils::get_invM_V0(Rec_V0VertexV0Object)") # vector of invariant masses of all reconstructed V0

            .Define("EVT_nVtx",                 "1+EVT_nSV+EVT_nV0")

            #combine for all SV d2PV
            .Define("Rec_SVV0_d2PV",               "ROOT::VecOps::Concatenate(Rec_SV_d2PV, Rec_V0_d2PV)")
            .Define("Rec_SVV0_x",               "ROOT::VecOps::Concatenate(Rec_SV_x, Rec_V0_x)")
            .Define("Rec_SVV0_y",               "ROOT::VecOps::Concatenate(Rec_SV_y, Rec_V0_y)")
            .Define("Rec_SVV0_z",               "ROOT::VecOps::Concatenate(Rec_SV_z, Rec_V0_z)")
            .Define("Rec_SVV0_p",               "ROOT::VecOps::Concatenate(Rec_SV_z, Rec_V0_p)")
            .Define("Rec_PVSV_ntracks",           "ROOT::VecOps::Concatenate(ROOT::VecOps::RVec<float>{Rec_PV_ntracks}, Rec_SV_ntracks)")



            ##############################################
            # variables now updating momentum from non-seeded fit
            #############################################

            # now update reco momentum based on the rec vertex position
            .Define("RecoParticlesPIDAtPV",  f"myUtils::get_RP_atVertex(RecoParticlesPID, ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>{{Rec_PrimaryVertexObject}})")
            .Define("RecoParticlesPIDAtPVSV",  "myUtils::get_RP_atVertex(RecoParticlesPIDAtPV, Rec_SecondaryVertexObject)")
            .Define("RecoParticlesPIDAtPVSVV0",  "myUtils::get_RP_atVertex(RecoParticlesPIDAtPVSV, Rec_V0VertexObject)")

            .Define("Rec_e_vfit",         "ReconstructedParticle::get_e(RecoParticlesPIDAtPVSVV0)")
            .Define("Rec_p_vfit",         "ReconstructedParticle::get_p(RecoParticlesPIDAtPVSVV0)")
            .Define("Rec_pt_vfit",        "ReconstructedParticle::get_pt(RecoParticlesPIDAtPVSVV0)")
            .Define("Rec_px_vfit",        "ReconstructedParticle::get_px(RecoParticlesPIDAtPVSVV0)")
            .Define("Rec_py_vfit",        "ReconstructedParticle::get_py(RecoParticlesPIDAtPVSVV0)")
            .Define("Rec_pz_vfit",        "ReconstructedParticle::get_pz(RecoParticlesPIDAtPVSVV0)")

            .Define("EVT_ThrustInfoNoPointing_vfit",     'Algorithms::minimize_thrust("Minuit2","Migrad")(Rec_px_vfit, Rec_py_vfit, Rec_pz_vfit)') 
            .Define("EVT_ThrustCosThetaNoPointing_vfit", "Algorithms::getAxisCosTheta(EVT_ThrustInfoNoPointing_vfit, Rec_px_vfit, Rec_py_vfit, Rec_pz_vfit)")
            .Define("EVT_ThrustInfo_vfit",               "Algorithms::getThrustPointing(1.)(EVT_ThrustCosThetaNoPointing_vfit, Rec_e_vfit, EVT_ThrustInfoNoPointing_vfit)")
            .Define("Rec_thrustCosTheta_vfit",           "Algorithms::getAxisCosTheta(EVT_ThrustInfo_vfit, Rec_px_vfit, Rec_py_vfit, Rec_pz_vfit)")
            .Define("Rec_in_hemisEmin_vfit",             "myUtils::get_RP_inHemis(1)(Rec_thrustCosTheta_vfit)")

            .Define("EVT_unitThrust_x_vfit",            "myUtils::norm_RVec_x(EVT_ThrustInfo_vfit.at(1),EVT_ThrustInfo_vfit.at(3),EVT_ThrustInfo_vfit.at(5))")
            .Define("EVT_unitThrust_y_vfit",            "myUtils::norm_RVec_x(EVT_ThrustInfo_vfit.at(3),EVT_ThrustInfo_vfit.at(1),EVT_ThrustInfo_vfit.at(5))")
            .Define("EVT_unitThrust_z_vfit",            "myUtils::norm_RVec_x(EVT_ThrustInfo_vfit.at(5),EVT_ThrustInfo_vfit.at(3),EVT_ThrustInfo_vfit.at(1))")
            
            .Define("EVT_sum_Rec_px_vfit",  "ROOT::VecOps::Sum(Rec_px_vfit)")
            .Define("EVT_sum_Rec_py_vfit",  "ROOT::VecOps::Sum(Rec_py_vfit)")
            .Define("EVT_sum_Rec_pz_vfit",  "ROOT::VecOps::Sum(Rec_pz_vfit)")
            
            .Define("EVT_p_vfit",        "sqrt(EVT_sum_Rec_px_vfit*EVT_sum_Rec_px_vfit+EVT_sum_Rec_py_vfit*EVT_sum_Rec_py_vfit+EVT_sum_Rec_pz_vfit*EVT_sum_Rec_pz_vfit)")
            .Define("EVT_e_vfit",        "ROOT::VecOps::Sum(Rec_e_vfit)")


            ########################################################
            # Defining additional variables for case where dont update momentum from vertexing
            ########################################################

            .Define("Rec_e_raw",         "ReconstructedParticle::get_e(RecoParticlesPID)")
            .Define("Rec_p_raw",         "ReconstructedParticle::get_p(RecoParticlesPID)")
            .Define("Rec_pt_raw",        "ReconstructedParticle::get_pt(RecoParticlesPID)")
            .Define("Rec_px_raw",        "ReconstructedParticle::get_px(RecoParticlesPID)")
            .Define("Rec_py_raw",        "ReconstructedParticle::get_py(RecoParticlesPID)")
            .Define("Rec_pz_raw",        "ReconstructedParticle::get_pz(RecoParticlesPID)")

            .Define("EVT_ThrustInfoNoPointing_raw",     'Algorithms::minimize_thrust("Minuit2","Migrad")(Rec_px_raw, Rec_py_raw, Rec_pz_raw)') 
            .Define("EVT_ThrustCosThetaNoPointing_raw", "Algorithms::getAxisCosTheta(EVT_ThrustInfoNoPointing_raw, Rec_px_raw, Rec_py_raw, Rec_pz_raw)")
            .Define("EVT_ThrustInfo_raw",               "Algorithms::getThrustPointing(1.)(EVT_ThrustCosThetaNoPointing_raw, Rec_e_raw, EVT_ThrustInfoNoPointing_raw)")
            .Define("Rec_thrustCosTheta_raw",           "Algorithms::getAxisCosTheta(EVT_ThrustInfo_raw, Rec_px_raw, Rec_py_raw, Rec_pz_raw)")
            .Define("Rec_in_hemisEmin_raw",             "myUtils::get_RP_inHemis(1)(Rec_thrustCosTheta_raw)")

            .Define("EVT_unitThrust_x_raw",            "myUtils::norm_RVec_x(EVT_ThrustInfo_raw.at(1),EVT_ThrustInfo_raw.at(3),EVT_ThrustInfo_raw.at(5))")
            .Define("EVT_unitThrust_y_raw",            "myUtils::norm_RVec_x(EVT_ThrustInfo_raw.at(3),EVT_ThrustInfo_raw.at(1),EVT_ThrustInfo_raw.at(5))")
            .Define("EVT_unitThrust_z_raw",            "myUtils::norm_RVec_x(EVT_ThrustInfo_raw.at(5),EVT_ThrustInfo_raw.at(3),EVT_ThrustInfo_raw.at(1))")
            
            .Define("EVT_sum_Rec_px_raw",  "ROOT::VecOps::Sum(Rec_px_raw)")
            .Define("EVT_sum_Rec_py_raw",  "ROOT::VecOps::Sum(Rec_py_raw)")
            .Define("EVT_sum_Rec_pz_raw",  "ROOT::VecOps::Sum(Rec_pz_raw)")
            
            .Define("EVT_p_raw",        "sqrt(EVT_sum_Rec_px_raw*EVT_sum_Rec_px_raw+EVT_sum_Rec_py_raw*EVT_sum_Rec_py_raw+EVT_sum_Rec_pz_raw*EVT_sum_Rec_pz_raw)")
            .Define("EVT_e_raw",        "ROOT::VecOps::Sum(Rec_e_raw)")





            ######################################################
            # EVT vars that dont care abt vertex fit and momentum update from seeded fit
            ######################################################
            .Define("EVT_ID", "rdfentry_") 
            .Define("EVT_n", "ReconstructedParticles.size()")
            .Define("EVT_nCharged", "Sum(Rec_q != 0)")
            .Define("EVT_nNeutral", "Sum(Rec_q == 0)")



            ###################################################################
            #Reco variables that use RecoParticlesPIDAtVertex but are not affected by momentum update from seeded fit
            ####################################################################
            
            #------------------------------------------------------------------
            ##           IP-like track vars           
            #------------------------------------------------------------------
            .Define("Rec_track_d0",      "ReconstructedParticle2Track::getRP2TRK_D0(RecoParticlesPIDAtVertex, EFlowTrack_1)")
            .Define("Rec_track_d0_sig",  "ReconstructedParticle2Track::getRP2TRK_D0_sig(RecoParticlesPIDAtVertex, EFlowTrack_1)")
            .Define("Rec_track_z0",      "ReconstructedParticle2Track::getRP2TRK_Z0(RecoParticlesPIDAtVertex, EFlowTrack_1)")
            .Define("Rec_track_z0_sig",  "ReconstructedParticle2Track::getRP2TRK_Z0_sig(RecoParticlesPIDAtVertex, EFlowTrack_1)")


            .Define("Rec_track_absd0",      "myUtils::abs_RVec(Rec_track_d0)")
            .Define("Rec_track_absnormd0",  "myUtils::abs_RVec(Rec_track_d0_sig)")
            .Define("Rec_track_absz0",      "myUtils::abs_RVec(Rec_track_z0)")
            .Define("Rec_track_absnormz0",  "myUtils::abs_RVec(Rec_track_z0_sig)")

            # Store total number of tracks
            .Define("Rec_track_n",       "float(ReconstructedParticle2Track::getTK_n(EFlowTrack_1))")

            #Adding IPs corrected so that from PV rather than 000
            .Define("Rec_PV_TLorentz",     "TLorentzVector(Rec_PrimaryVertex.position.x, Rec_PrimaryVertex.position.y, Rec_PrimaryVertex.position.z, 0.)") #time component not used so fill with 0.
            .Define("MC_PV_TLorentz",     "TLorentzVector(MC_PrimaryVertex.X(), MC_PrimaryVertex.Y(), MC_PrimaryVertex.Z(), 0.)")
            .Define("Rec_track_d0_fromRecPV",     "ReconstructedParticle2Track::XPtoPar_dxy(RecoParticlesPIDAtVertex, EFlowTrack_1, Rec_PV_TLorentz, magFieldBz.at(0))")
            .Define("Rec_track_z0_fromRecPV",     "ReconstructedParticle2Track::XPtoPar_dz(RecoParticlesPIDAtVertex, EFlowTrack_1,  Rec_PV_TLorentz, magFieldBz.at(0))")
            .Define("Rec_track_d0_fromMCPV",     "ReconstructedParticle2Track::XPtoPar_dxy(RecoParticlesPIDAtVertex, EFlowTrack_1, MC_PV_TLorentz, magFieldBz.at(0))")
            .Define("Rec_track_z0_fromMCPV",     "ReconstructedParticle2Track::XPtoPar_dz(RecoParticlesPIDAtVertex, EFlowTrack_1,  MC_PV_TLorentz, magFieldBz.at(0))")



            ######################################
            # Association of reco particle to true MC particle
            ######################################

            .Define("MC_fromRP",           "myUtils::get_MCObject_fromRP(MCRecoAssociationsRec, MCRecoAssociationsGen, RecoParticlesPIDAtVertex, Particle)")
            .Define("Rec_true_PDG",        "MCParticle::get_pdg(MC_fromRP)")  # this true ID from MC
            .Define("Rec_true_e",          "MCParticle::get_e(MC_fromRP)")
            .Define("Rec_true_m",          "MCParticle::get_mass(MC_fromRP)")
            .Define("Rec_true_q",          "MCParticle::get_charge(MC_fromRP)")
            .Define("Rec_true_p",          "MCParticle::get_p(MC_fromRP)")
            .Define("Rec_true_pt",         "MCParticle::get_pt(MC_fromRP)")
            .Define("Rec_true_px",         "MCParticle::get_px(MC_fromRP)")
            .Define("Rec_true_py",         "MCParticle::get_py(MC_fromRP)")
            .Define("Rec_true_pz",         "MCParticle::get_pz(MC_fromRP)")
            .Define("Rec_true_eta",        "MCParticle::get_eta(MC_fromRP)")
            .Define("Rec_true_phi",        "MCParticle::get_phi(MC_fromRP)")
            .Define("Rec_true_orivtx_x",   "MCParticle::get_vertex_x(MC_fromRP)")
            .Define("Rec_true_orivtx_y",   "MCParticle::get_vertex_y(MC_fromRP)")
            .Define("Rec_true_orivtx_z",   "MCParticle::get_vertex_z(MC_fromRP)")

            # RecoP true history (mothers and gmothers)
            .Define("True_ParentInfo",     "myUtils::get_MCParentandGParent_fromRP(MCRecoAssociationsRec, MCRecoAssociationsGen, ParticleParents, RecoParticlesPIDAtVertex, Particle)")   # INTERMEDIATE
            .Define("Rec_true_M1",         "True_ParentInfo.at(0)")
            .Define("Rec_true_M2",         "True_ParentInfo.at(1)")
            .Define("Rec_true_M1ofM1",     "True_ParentInfo.at(2)")
            .Define("Rec_true_M2ofM1",     "True_ParentInfo.at(3)")
            .Define("Rec_true_M1ofM2",     "True_ParentInfo.at(4)")
            .Define("Rec_true_M2ofM2",     "True_ParentInfo.at(5)")
            

            ####################################################
            ## Defining additional variables for flavour tagging (agnostic to momentum update)
            ####################################################


            #---------------------------------------------------------------
            # Define qtag (ie. if B0 or B0b and equiv for Bs)
            #----------------------------------------------------------------
            #Thrust axis infor for MC particles
            .Define("MC_thrustCosTheta",           "Algorithms::getAxisCosTheta(EVT_ThrustInfo, MC_px, MC_py, MC_pz)")
            .Define("MC_in_hemisEmin",             "myUtils::get_RP_inHemis(1)(MC_thrustCosTheta)") # just gives 1 or 0 based on whether in min hemisphere or not
            # Get the production flavour of B0/Bs0 mesons in signal hemisphere
            .Define("EVT_sigB_MCprodFlav", "myUtils::get_B_prod_flav_from_nunu(MC_PDG, MC_M1)")
            #use this to define qtag (saved flav as a check)
            .Define("EVT_sigB_MCqTag",  "EVT_sigB_MCprodFlav > 0 ? 1 : (EVT_sigB_MCprodFlav < 0 ? -1 : 0)")

            #also get true momentum
            .Define("EVT_sigB_MCp", "myUtils::get_fsB_MC_var(MC_p,MC_PDG, MC_M1)")
            .Define("EVT_sigB_MCpx", "myUtils::get_fsB_MC_var(MC_px,MC_PDG, MC_M1)")
            .Define("EVT_sigB_MCpy", "myUtils::get_fsB_MC_var(MC_py,MC_PDG, MC_M1)")
            .Define("EVT_sigB_MCpz", "myUtils::get_fsB_MC_var(MC_pz,MC_PDG, MC_M1)")


            #---------------------------------------
            # Number of neutrals
            #---------------------------------------

            #Get the indices of all photons (PDG ID == 22)
            .Define("Rec_photon_indices", "myUtils::sel_PID(22)(RecoParticlesPIDAtVertex)") # intermediate
            # Count photons
            .Define("EVT_nPhotons",       "float(Rec_photon_indices.size())")

            # Get number of KS (various assumption levels)
            #number of MC KS
            .Define("MC_KS", "MCParticle::sel_pdgID(310, false)(Particle)") #m_abs redundant for KS where PID == 310 ##INTERMEDIATE STATE 
            .Define("EVT_nKS_MC", "MC_KS.size()" )

            # Get number of reco KS (still cheating without using combinatorics - ie count if have two pi with same origin vertex and parent as KS). As seen in B2Inv study, some inefficiency to hugh energy KS which decay beyon tracker
            .Define("MC_recParticle_indx", "myUtils::get_RP_idx_from_MC(MCRecoAssociationsRec, MCRecoAssociationsGen, Particle)")
            .Define("Rec_true_KS", "myUtils::get_rec_true_KS(Particle, ParticleChildren , MC_recParticle_indx)")
            .Define("EVT_nKS_recTrue","Rec_true_KS.size()")

            #-----------------------------------------
            # Reco variables for PID
            #-----------------------------------------
            #Define dNdx for PID tool
            # First get the track states  associated with  reconstructed particles
            .Define("Rec_trackStates", "ReconstructedParticle2Track::getRP2TRK(RecoParticlesPIDAtVertex, EFlowTrack_1)")
            # Get the dN/dx for those tracks (maintains the same array indices as d0, z0 etc)
            .Define("Rec_track_dNdx", "ReconstructedTrack::tracks_dNdx(Rec_trackStates, EFlowTrack_1, EFlowTrack, EFlowTrack_2)")
            .Define("Rec_track_dNdx_paddedNeutrals", "ReconstructedParticle2Track::getRP2TRK_dNdX(RecoParticlesPIDAtVertex, Rec_trackStates, EFlowTrack_1, EFlowTrack, EFlowTrack_2)")

            #Also define TOF using TrackerHits: RVec<edm4hep::TrackerHit3DData> object
            .Define("Rec_track_TOF","ReconstructedTrack::tracks_TOF(Rec_trackStates, EFlowTrack_1, EFlowTrack, TrackerHits)")
            .Define("Rec_track_length", "ReconstructedTrack::tracks_length(Rec_trackStates, EFlowTrack_1, EFlowTrack_L)")

            .Define("Rec_track_TOF_paddedNeutrals","ReconstructedParticle2Track::getRP2TRK_TOF(RecoParticlesPIDAtVertex, Rec_trackStates, EFlowTrack_1, EFlowTrack, TrackerHits)")
            .Define("Rec_track_length_paddedNeutrals", "ReconstructedParticle2Track::getRP2TRK_length(RecoParticlesPIDAtVertex, Rec_trackStates, EFlowTrack_1, EFlowTrack_L)")

            #---------------------------------------------
            #Also adding some calorimetry info
            #--------------------------------------------

            #Currently not in use as CalorimeterHits doesn't actually include any energy inormation!!!!!!!
            #.Define("Calo_hits_x", "CaloNtupleizer::getCaloHit_x(CalorimeterHits)")
            #.Define("Calo_hits_y", "CaloNtupleizer::getCaloHit_y(CalorimeterHits)")
            #.Define("Calo_hits_z", "CaloNtupleizer::getCaloHit_z(CalorimeterHits)")
            #.Define("Calo_hits_e", "CaloNtupleizer::getCaloHit_energy(CalorimeterHits)")


        )



        # Add preliminary cut to ensure PV present     
        if cfg.run_mode == 'BflavTag_full':
            df3 = (
                df2 
                #This is clearest way to cut on this as physically impossible to reconstruct a 1 track vertex
                .Filter("Rec_PV_ntracks>1")
                )   
            return df3
 
        
           
    #required output function - actually saves df with selected branches sepcified in yaml
    def output():
         # Get the output branchList from the config YAML file
         with open(cfg.fccana_opts['yamlPath']) as stream:
             yaml = safe_load(stream)
             branchList = yaml[cfg.fccana_opts['outputBranches'][cfg.run_mode]]
             print(f"----> INFO:")
             print(f"            Output branch list used = {cfg.fccana_opts['outputBranches'][cfg.run_mode]}")
 
         return branchList         
        

