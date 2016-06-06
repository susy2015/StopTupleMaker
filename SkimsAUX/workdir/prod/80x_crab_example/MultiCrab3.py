#!/usr/bin/env python
# encoding: utf-8

# File        : MultiCrab3.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2015 Apr 01
#
# Description :


import copy, os, time, sys

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientUtilities import colors
from WMCore.Configuration import saveConfigurationFile
from crab3Config import config as config
from multiprocessing import Process

workArea = 'crabProdv4p0'
outDir =  '/store/group/lpcsusyhad/Spring16_8xX_Jun_2016_Ntp_v1X'
Pubname = 'Spring16_8XX_Dec_2016_Ntp_v1p0'
json_25ns = 'Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt'
# Use the common keyword to select the samples you'd like to submit
# ALL: all of them; NONE: none of them; TEST: test printing out the crab3 config or disable actual submission; STATUS: check job status
# TTJets, WJetsToLNu, ZJetsToNuNu, DYJetsToLL, QCD, TTW, TTZ, ST_tW, SMS, HTMHT, SingleMuon, SingleElectron, DoubleMuon, DoubleEG
# Can be any of the combinations
#selSubmitKey = 'TEST STATUS TTJets' # 'TEST STATUS': no submission of jobs but rather checking crab job status related to the TTJets. If jobs failed, automatically resubmit them.
selSubmitKey = 'TTJets'
doAutoMonitor = False

## Format: keyword : IsData, fulldatasetname, unitperjob
jobslist = {
    # TTbar
#    'TTJets'                                 : [False, '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM', 1],
    'TTJets_SingleLeptFromT'                 : [False, '/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'TTJets_SingleLeptFromT_ext1'            : [False, '/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'TTJets_SingleLeptFromTbar'              : [False, '/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'TTJets_SingleLeptFromTbar_ext1'         : [False, '/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'TTJets_HT-600to800'                     : [False, '/TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'TTJets_HT-800to1200'                    : [False, '/TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'TTJets_HT-1200to2500'                   : [False, '/TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'TTJets_HT-2500toInf'                    : [False, '/TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'TTJets_DiLept'                          : [False, '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'TTJets_DiLept_ext1'                     : [False, '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],

    # WJets,
    'WJetsToLNu'                             : [False, '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM', 1],
    'WJetsToLNu_HT-100To200'                 : [False, '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v2/MINIAODSIM', 1],
    'WJetsToLNu_HT-200To400'                 : [False, '/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM', 1],
    'WJetsToLNu_HT-200To400_ext1'            : [False, '/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'WJetsToLNu_HT-400To600'                 : [False, '/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM', 1],
    'WJetsToLNu_HT-600To800'                 : [False, '/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM', 1],
    'WJetsToLNu_HT-800To1200'                : [False, '/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'WJetsToLNu_HT-800To1200_ext1'           : [False, '/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'WJetsToLNu_HT-1200To2500'               : [False, '/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM', 1],
  'WJetsToLNu_HT-1200To2500_ext1'            : [False, '/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'WJetsToLNu_HT-2500ToInf'                : [False, '/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'WJetsToLNu_HT-600ToInf'                 : [False,'/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    # Zinv,
#    'ZJetsToNuNu_HT-100To200'                : [False, '/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM', 1],
    'ZJetsToNuNu_HT-200To400'                : [False, '/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
#    'ZJetsToNuNu_HT-400To600'                : [False, '/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM', 1],
#    'ZJetsToNuNu_HT-600ToInf'                : [False, '/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/MINIAODSIM', 1],
    'ZJetsToNuNu_HT-600To800'                : [False, '/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'ZJetsToNuNu_HT-800To1200'                : [False, '/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'ZJetsToNuNu_HT-1200To2500'                : [False, '/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'ZJetsToNuNu_HT-2500ToInf'                : [False, '/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    # DYJets,
    'DYJetsToLL-M-50'                        : [False, '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'DYJetsToLL-M-50_extv1'                  : [False, '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'DYJetsToLL_M-50_HT-100to200'            : [False, '/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'DYJetsToLL_M-50_HT-200to400'            : [False, '/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'DYJetsToLL_M-50_HT-400to600'            : [False, '/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
#    'DYJetsToLL_M-50_HT-600toInf'            : [False, '/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM', 1],

    # QCD,
    'QCD_HT100to200'                         : [False, '/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT200to300'                         : [False, '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT200to300_ext1'                    : [False, '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'QCD_HT300to500'                         : [False, '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT300to500_ext1'                         : [False, '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'QCD_HT500to700'                         : [False, '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT500to700_ext1'                         : [False, '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'QCD_HT700to1000'                        : [False, '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT700to1000_ext1'                        : [False, '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'QCD_HT1000to1500'                       : [False, '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT1000to1500_ext1'                       : [False, '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'QCD_HT1500to2000'                       : [False, '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT1500to2000_ext'                       : [False, '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'QCD_HT2000toInf'                        : [False, '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'QCD_HT2000toInf_ext1'                        : [False, '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v2/MINIAODSIM', 1],

    #ttH
    'ttHJetToNonbb'                          :[False, '/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3_ext1-v1/MINIAODSIM', 1],
    'ttHJetTobb'                             :[False, '/ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3_ext3-v1/MINIAODSIM', 1],

    #DiBoson                                 
    'ZZ'                                     :[False, '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM/MINIAODSIM', 1],
    'WZ'                                     :[False, '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'WW'                                     :[False, '/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    #TriBoson
    'WZZ'                                    :[False, '/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'WWZ'                                    :[False, '/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'ZZZ'                                    :[False, '/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    #tt-Gamma
    'TTGJets'                                :[False, '/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    # ttW
    'TTWJetsToQQ'                            : [False, '/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'TTWJetsToLNu'                           : [False, '/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    # ttZ
    'TTZToQQ'                                : [False, '/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'TTZToLLNuNu'                            : [False, '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v3/MINIAODSIM', 1],

    # single top
    'ST_tW_top_5f_inclusiveDecays'           : [False, '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'ST_tW_antitop_5f_inclusiveDecays'       : [False, '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    # Signals
    'SMS-T1tttt_mGluino-1200_mLSP-800'       : [False, '/SMS-T1tttt_mGluino-1200_mLSP-800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'SMS-T1tttt_mGluino-1500_mLSP-100'       : [False, '/SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'SMS-T1bbbb_mGluino-1000_mLSP-900'       : [False, '/SMS-T1bbbb_mGluino-1000_mLSP-900_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'SMS-T1bbbb_mGluino-1500_mLSP-100'       : [False, '/SMS-T1bbbb_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'SMS-T2tt_mStop-500_mLSP-325'            : [False, '/SMS-T2tt_mStop-500_mLSP-325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],
    'SMS-T2tt_mStop-850_mLSP-100'            : [False, '/SMS-T2tt_mStop-850_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM', 1],

    # Data
#    'HTMHT-Run2015B-05Oct2015'               : [True, '/HTMHT/Run2015B-05Oct2015-v1/MINIAOD', 10],
    #'HTMHT-Run2015C-25ns-05Oct2015'          : [True, '/HTMHT/Run2015C_25ns-05Oct2015-v1/MINIAOD', 10],
    #'HTMHT-Run2015D-05Oct2015'               : [True, '/HTMHT/Run2015D-05Oct2015-v1/MINIAOD', 10],
    'HTMHT-Run2016B-PromptReco-v2'           : [True, '/HTMHT/Run2016B-PromptReco-v2/MINIAOD', 1],
    'HTMHT-Run2016B-PromptReco'              : [True, '/HTMHT/Run2016B-PromptReco-v1/MINIAOD', 1],

#    'MET-Run2015B-05Oct2015'                 : [True, '/MET/Run2015B-05Oct2015-v1/MINIAOD', 10],
    #'MET-Run2015C-25ns-05Oct2015'            : [True, '/MET/Run2015C_25ns-05Oct2015-v1/MINIAOD', 10],
    #'MET-Run2015D-05Oct2015'                 : [True, '/MET/Run2015D-05Oct2015-v1/MINIAOD', 10],
    'MET-Run2016B-PromptReco-v2'                : [True, '/MET/Run2016B-PromptReco-v2/MINIAOD', 1],
    'MET-Run2016B-PromptReco-v1'                : [True, '/MET/Run2016B-PromptReco-v1/MINIAOD', 1],

#    'JetHT-Run2015B-05Oct2015'               : [True, '/JetHT/Run2015B-05Oct2015-v1/MINIAOD', 10],
#    'JetHT-Run2015C-25ns-05Oct2015'          : [True, '/JetHT/Run2015C_25ns-05Oct2015-v1/MINIAOD', 10],
#    'JetHT-Run2015D-05Oct2015'               : [True, '/JetHT/Run2015D-05Oct2015-v1/MINIAOD', 10],
#    'JetHT-Run2015D-PromptReco'              : [True, '/JetHT/Run2015D-PromptReco-v4/MINIAOD', 10],

#    'SingleMuon-Run2015B-05Oct2015'          : [True, '/SingleMuon/Run2015B-05Oct2015-v1/MINIAOD', 10],
    #'SingleMuon-Run2015C-25ns-05Oct2015'     : [True, '/SingleMuon/Run2015C_25ns-05Oct2015-v1/MINIAOD', 10],
    #'SingleMuon-Run2015D-05Oct2015'          : [True, '/SingleMuon/Run2015D-05Oct2015-v1/MINIAOD', 10],
    'SingleMuon-Run2016B-PromptReco-v2'         : [True, '/SingleMuon/Run2016B-PromptReco-v2/MINIAOD', 1],
    'SingleMuon-Run2016B-PromptReco-v1'         : [True, '/SingleMuon/Run2016B-PromptReco-v1/MINIAOD', 1],

#    'SingleElectron-Run2015B-05Oct2015'      : [True, '/SingleElectron/Run2015B-05Oct2015-v1/MINIAOD', 10],
#    'SingleElectron-Run2015C-25ns-05Oct2015' : [True, '/SingleElectron/Run2015C_25ns-05Oct2015-v1/MINIAOD', 10],
#    'SingleElectron-Run2015D-05Oct2015'      : [True, '/SingleElectron/Run2015D-05Oct2015-v1/MINIAOD', 10],
#    'SingleElectron-Run2015D-PromptReco'     : [True, '/SingleElectron/Run2015D-PromptReco-v4/MINIAOD', 10],

#    'DoubleMuon-Run2015B-05Oct2015'          : [True, '/DoubleMuon/Run2015B-05Oct2015-v1/MINIAOD', 10],
#    'DoubleMuon-Run2015C-25ns-05Oct2015'     : [True, '/DoubleMuon/Run2015C_25ns-05Oct2015-v1/MINIAOD', 10],
#    'DoubleMuon-Run2015D-05Oct2015'          : [True, '/DoubleMuon/Run2015D-05Oct2015-v1/MINIAOD', 10],
#    'DoubleMuon-Run2015D-PromptReco'         : [True, '/DoubleMuon/Run2015D-PromptReco-v4/MINIAOD', 10],

#    'DoubleEG-Run2015B-05Oct2015'            : [True, '/DoubleEG/Run2015B-05Oct2015-v1/MINIAOD', 10],
#    'DoubleEG-Run2015C-25ns-05Oct2015'       : [True, '/DoubleEG/Run2015C_25ns-05Oct2015-v1/MINIAOD', 10],
#    'DoubleEG-Run2015D-05Oct2015'            : [True, '/DoubleEG/Run2015D-05Oct2015-v1/MINIAOD', 10],
#    'DoubleEG-Run2015D-PromptReco'           : [True, '/DoubleEG/Run2015D-PromptReco-v4/MINIAOD', 10],
}
tasklist = {}

def MonitoringJobs(tasklist):
    while True:
        sumFailed = 0
        sumComp = 0
        for request, name in tasklist.items():
           dirname = './%s/crab_%s' % (workArea, name)
           fulldir = os.path.abspath(dirname)
           try:
              results = crabCommand('status', dir=fulldir)
              if 'FAILED' in results['status']:
                 sumFailed += 1
              if 'COMPLETED' in results['status']:
                 sumComp += 1  
              print "For task", request, "the job states are", results['jobsPerStatus']
              status = results['jobsPerStatus']
              if 'failed' in status:
                 print "failed : ", status['failed']
                 crabCommand('resubmit', dir=fulldir)
           except:
              pass
           time.sleep(2)

        print "\n\n", colors.RED, "sumFailed : ", sumFailed, "  sumComp : ", sumComp, " RE-CHECKING EVERY TASK...\n\n", colors.NORMAL
        if sumFailed == 0 and sumComp == len(tasklist):
           break;

def CreateMonitorList(tasklist):
    monList = open("monList_"+workArea+".txt", 'w')
    for key in sorted(tasklist):
       dirname = './%s/crab_%s' % (workArea, tasklist[key])
       fulldir = os.path.abspath(dirname)
       monList.write("crab status "+fulldir+"\n")
    monList.close()


def SubmitJob(key, value):
    doAll = False 
    doTest = False

    allSelKeys = selSubmitKey.split()

    if selSubmitKey.find('NONE') != -1:
       print "Nothing to be done!"
       sys.exit()
    if selSubmitKey.find('ALL') != -1:
       doAll = True
    if selSubmitKey.find('TEST') != -1:
       doTest = True

    doThis = doAll

    if not doAll:
        for selKey in allSelKeys:
            if key.find(selKey) != -1:
                doThis = True
                break;

    if not doThis:
        return

    tempconfig = copy.deepcopy(config)
    tempconfig.General.requestName = key
    tempconfig.General.workArea = workArea
    tempconfig.Data.outputDatasetTag = Pubname + "_" + key
    tempconfig.Data.outLFNDirBase = outDir

    if len(value) <3:
        print "Not enough argument for %s" % key
        raise  AssertionError()
    if value[0]: # Data
        if key.find('Run2015C') != -1:
            tempconfig.JobType.pyCfgParams = ['mcInfo=0', 'GlobalTag=74X_dataRun2_v4', 'specialFix=JEC', 'jecDBname=Summer15_25nsV6_DATA', 'externalFilterList=csc2015_Dec01.txt.tar.gz,ecalscn1043093_Dec01.txt.tar.gz']
            tempconfig.JobType.inputFiles = [json_25ns, 'Summer15_25nsV6_DATA.db', 'csc2015_Dec01.txt.tar.gz', 'ecalscn1043093_Dec01.txt.tar.gz']
            tempconfig.Data.splitting = 'LumiBased'
            tempconfig.Data.lumiMask = json_25ns
        elif key.find('Run2015D-05Oct2015') != -1:
            tempconfig.JobType.pyCfgParams = ['mcInfo=0', 'GlobalTag=74X_dataRun2_reMiniAOD_v0', 'specialFix=JEC', 'jecDBname=Summer15_25nsV6_DATA', 'externalFilterList=csc2015_Dec01.txt.tar.gz,ecalscn1043093_Dec01.txt.tar.gz']
            tempconfig.JobType.inputFiles = [json_25ns, 'Summer15_25nsV6_DATA.db', 'csc2015_Dec01.txt.tar.gz', 'ecalscn1043093_Dec01.txt.tar.gz']
            tempconfig.Data.splitting = 'LumiBased'
            tempconfig.Data.lumiMask = json_25ns
        elif key.find('Run2015B-PromptReco') != -1:
            tempconfig.JobType.pyCfgParams = ['mcInfo=0', 'GlobalTag=74X_dataRun2_reMiniAOD_v0', 'specialFix=JEC', 'jecDBname=Summer15_25nsV6_DATA', 'externalFilterList=csc2015_Dec01.txt.tar.gz,ecalscn1043093_Dec01.txt.tar.gz']
            tempconfig.JobType.inputFiles = [json_25ns, 'Summer16_25nsV1_DATA.db', 'csc2015_Dec01.txt.tar.gz', 'ecalscn1043093_Dec01.txt.tar.gz']#Not ready in 2016 yes the json_25ns Summer16_25 file
            tempconfig.Data.splitting = 'LumiBased'
            tempconfig.Data.lumiMask = json_25ns
        elif key.find('Run2016B-PromptReco') != -1:
            tempconfig.JobType.pyCfgParams = ['mcInfo=0', 'GlobalTag=80X_dataRun2_Prompt_v8']#, 'specialFix=JEC', 'jecDBname=Summer15_25nsV6_DATA', 'externalFilterList=csc2015_Dec01.txt.tar.gz,ecalscn1043093_Dec01.txt.tar.gz']
           # tempconfig.JobType.inputFiles = [json_25ns, 'Summer16_25nsV1_DATA.db', 'csc2015_Dec01.txt.tar.gz', 'ecalscn1043093_Dec01.txt.tar.gz']# Not ready in 2016 yes
            tempconfig.Data.splitting = 'LumiBased'
            tempconfig.Data.lumiMask = json_25ns
        else:
            pass
    else:
        tempconfig.JobType.pyCfgParams = ['mcInfo=1', 'GlobalTag=80X_mcRun2_asymptotic_2016_v3']#, 'specialFix=JEC', 'jecDBname=Summer15_25nsV6_MC']
        #tempconfig.JobType.inputFiles = ['Spring16_25nsV1_MC.db']#'Summer15_25nsV6_MC.db']
        tempconfig.Data.splitting = 'FileBased'

    tempconfig.Data.inputDataset = value[1].strip()
    tempconfig.Data.unitsPerJob = value[2]

    if value[0] and len(value) > 3:
        tempconfig.Data.lumiMask = value[3]

    # Submitting jobs
    if doTest:
        saveConfigurationFile(tempconfig, workArea+"/test/"+key+"_test_cfg.py")
        tasklist["crab_"+key] = key
    else:
        results = crabCommand('submit', config = tempconfig)
        tasklist[results['uniquerequestname']] = key
    del tempconfig

if __name__ == "__main__":

    if not os.path.exists(workArea):
       os.makedirs(workArea)
       os.makedirs(workArea+"/test")

    allSelKeys = selSubmitKey.split()

    for key, value in jobslist.items():
        p = Process(target=SubmitJob , args=(key, value))
        p.start()
        p.join()

        for selKey in allSelKeys:
           if key.find(selKey) != -1:
              tasklist["crab_"+key] = key
              break;

    CreateMonitorList(tasklist)

    doTest = False
    doCheckStatus = False

    if selSubmitKey.find('TEST') != -1:
       doTest = True
    if selSubmitKey.find('STATUS') != -1:
       doCheckStatus = True

    if doTest:
       print tasklist

    if (not doTest and doAutoMonitor) or doCheckStatus:
       MonitoringJobs(tasklist)
