import argparse
import os
import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd
import sys
from particle import Particle
import copy
from utils.analysis_utils import loadAO2D, applySelections, applySelections2, perform_roofit_fit, get_chi2_significance_sb, set_param

def compute_efficiency(nRec, nGen):
    eff = nRec / nGen
    effErr = eff * np.sqrt((1 - eff) / nRec)
    return eff, effErr

# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('configfile', 
                        metavar='text',
                        help='config file name with inputs and cut configurations')
    args = parser.parse_args()
    with open(args.configfile, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    # Load MC dataframes from files, and apply preliminary selections
    # We can have multiple files corresponding to different MC samples
    mcFileNames = cfg['fileNameMC']

    # initialize output files
    outFolder = cfg['outputFolder']
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    outfile_prefix = cfg['outputFilePrefix']
    # Load and skim MC dataframes
    if not isinstance(mcFileNames, list):
        mcFileNames = [mcFileNames]
    inTreeNameRec = cfg['treeNameRec']
    inTreeNameGen = cfg['treeNameGen']
    weights = cfg['mcWeights']
    dfRecList = [] 
    dfGenList = []
    for fileName in mcFileNames:
        dfRec = loadAO2D(fileName, inTreeNameRec)
        dfRecFiltered = applySelections(dfRec, cfg, isMC=True)
        dfGen = loadAO2D(fileName, inTreeNameGen)
        dfRecList.append(dfRecFiltered)
        dfGenList.append(dfGen)

    # Load data dataframe from file, and apply preliminary selections
    dataFileName = cfg['fileNameData']
    inTreeNameData = cfg['treeNameData']
    dfData = loadAO2D(dataFileName, inTreeNameData)
    dfDataFiltered = applySelections(dfData, cfg, isMC=False)
    
     # apply central cut selections except on BdtScore
    centralCuts = cfg['centralCuts']
    flagsToKeep = cfg['acceptFlags']
    dfDataCentralCut = applySelections2(dfDataFiltered, centralCuts)
    dfDataCentralCut = dfDataCentralCut[['fM','fPt','fPtBach0', 'fMlScoreBkgBach0', 'fMlScoreNonPromptBach0']]

    dfRecList = [applySelections2(dfRec, centralCuts, isMC=True, selectedFlags=flagsToKeep) for irec, dfRec in enumerate(dfRecList)]
   
    # Load configurations to test
    configurations = cfg['nonPromptCutvar']['configurations']
    
    fit_conf = cfg['rawYield']['default_fit_config']
    fit_mass_min = fit_conf['mass_min']
    fit_mass_max = fit_conf['mass_max']
    fit_signal_pdf = fit_conf['signal_pdf']
    fit_bkg_pdf = fit_conf['bkg_pdf']
    fit_params_def = fit_conf['parameters']

    pt_bins_min = cfg['cutVars']['pt']['min']
    pt_bins_max = cfg['cutVars']['pt']['max']
    
    pt_limits = pt_bins_min.copy()
    pt_limits.append(pt_bins_max[-1])
    pt_limits = np.array(pt_limits, np.float64)

    pt_bins_dau_min = cfg['cutVars']['ptBach0']['min']
    pt_bins_dau_max = cfg['cutVars']['ptBach0']['max']

    # Define 2 ROOT.TH2 histograms
    hist_rawyield_cutvar = []
    hist_eff_prompt_cutvar = []
    hist_eff_non_prompt_cutvar = []
    for icut, _ in enumerate(configurations):
        hist_rawyield_cutvar.append(
            ROOT.TH1F("hist_rawyield", ";#it{p}_{T} (GeV/#it{c}); raw yield", len(pt_bins_min), pt_limits))
        hist_eff_prompt_cutvar.append(
            ROOT.TH1F("hist_eff_prompt", ";#it{p}_{T} (GeV/#it{c}); raw yield", len(pt_bins_min), pt_limits))
        hist_eff_non_prompt_cutvar.append(
            ROOT.TH1F("hist_eff_non_prompt", ";#it{p}_{T} (GeV/#it{c}); raw yield", len(pt_bins_min), pt_limits))
        hist_rawyield_cutvar[icut].SetDirectory(0)
        hist_eff_prompt_cutvar[icut].SetDirectory(0)
        hist_eff_non_prompt_cutvar[icut].SetDirectory(0)
    
    centralMeans = [0]*len(pt_bins_min)
    centralSigmas = [0]*len(pt_bins_min)
    for i_conf, conf in enumerate(configurations):
        # 0) retrieve configurations
        bdt_cuts = conf['cuts']
        if len(bdt_cuts) != len(pt_bins_dau_max):
            print(f"ERROR: Number of BDT cuts is not equal to the number of D bacehelor pt bins, {len(bdt_cuts)}, {len(pt_bins_dau_max)}")
            sys.exit(1)
        fit_params_list = [copy.deepcopy(fit_params_def)]*len(pt_bins_min)
        conf_name = conf['name']
        outfile = ROOT.TFile(f"{outFolder}/{outfile_prefix}_{conf_name}.root", "RECREATE")
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_bins_min, pt_bins_max)):
        # split dataframes in pt bins
            dfDataPt = dfDataCentralCut[(dfDataCentralCut['fPt'] >= pt_min) & (dfDataCentralCut['fPt'] < pt_max)]
            mcRecPtList = [dfRec[(dfRec['fPt'] >= pt_min) & (dfRec['fPt'] < pt_max)] for dfRec in dfRecList]
            mcGenPtList = [dfGen[(dfGen['fPt'] >= pt_min) & (dfGen['fPt'] < pt_max)] for dfGen in dfGenList]
            # 1) apply BDT cuts, differential in bachelor pT
            dataSel = pd.DataFrame()
            mcRecSel = [pd.DataFrame(), pd.DataFrame()]
            for i, (pt_min_d, pt_max_d, cut) in enumerate(zip(pt_bins_dau_min, pt_bins_dau_max, bdt_cuts)):
                # apply selections
                dataSel = pd.concat([dataSel, dfDataPt[(dfDataPt['fPtBach0'] >= pt_min_d) & (dfDataPt['fPtBach0'] < pt_max_d) & (dfDataPt['fMlScoreNonPromptBach0'] > cut)]])
                mcRecSel = [pd.concat([mcRecSel[irec], dfRec[(dfRec['fPtBach0'] >= pt_min_d) & (dfRec['fPtBach0'] < pt_max_d) & (dfRec['fMlScoreNonPromptBach0'] > cut)]]) for irec, dfRec in enumerate(mcRecPtList)]
            # 2) compute efficiency
            nRecPrompt, nGenPrompt = 0, 0
            nRecNonPrompt, nGenNonPrompt = 0, 0
            for i, (dfRec, dfGen, w) in enumerate(zip(mcRecSel, mcGenPtList, cfg['mcWeights'])):
                nGenPrompt += w * len(dfGen[dfGen['fOrigin'] == 1])
                nRecPrompt += w * len(dfRec[(dfRec['fOrigin'] == 1)])
                nGenNonPrompt += w * len(dfGen[dfGen['fOrigin'] == 2])
                nRecNonPrompt += w * len(dfRec[(dfRec['fOrigin'] == 2)])
            nRecPrompt = nRecPrompt / np.sum(cfg['mcWeights'])
            nGenPrompt = nGenPrompt / np.sum(cfg['mcWeights'])
            effPrompt, effErrPrompt = compute_efficiency(nRecPrompt, nGenPrompt)
            nRecNonPrompt = nRecNonPrompt / np.sum(cfg['mcWeights'])
            nGenNonPrompt = nGenNonPrompt / np.sum(cfg['mcWeights'])
            effNonPrompt, effErrNonPrompt = compute_efficiency(nRecNonPrompt, nGenNonPrompt)
            hist_eff_prompt_cutvar[i_conf].SetBinContent(ipt + 1, effPrompt)
            hist_eff_prompt_cutvar[i_conf].SetBinError(ipt + 1, effErrPrompt)
            hist_eff_non_prompt_cutvar[i_conf].SetBinContent(ipt + 1, effNonPrompt)
            hist_eff_non_prompt_cutvar[i_conf].SetBinError(ipt + 1, effErrNonPrompt)
            
            # 3) compute raw yield
            if i_conf > 0:
                set_param(fit_params_list[ipt], "sigma", centralSigmas[ipt], centralSigmas[ipt], centralSigmas[ipt])
                set_param(fit_params_list[ipt], "mean", centralMeans[ipt], centralMeans[ipt], centralMeans[ipt])

            ws = perform_roofit_fit(dataSel, fit_signal_pdf, fit_bkg_pdf,  fit_params_list[ipt], fit_mass_min, fit_mass_max)
            chi2, significance, sb = get_chi2_significance_sb(ws, 100, outfile = outfile, name = f"c_{pt_min}_{pt_max}_cut_{conf['name']}")

            S = ws.var("nSig").getVal()
            S_err = ws.var("nSig").getError()
            sigma = ws.var("sigma").getVal()
            sigma_err = ws.var("sigma").getError()
            mean = ws.var("mean").getVal()
            mean_err = ws.var("mean").getError()
            if i_conf == 0:
                centralMeans[ipt] = mean
                centralSigmas[ipt] = sigma
            hist_rawyield_cutvar[i_conf].SetBinContent(ipt + 1, S)
            hist_rawyield_cutvar[i_conf].SetBinError(ipt + 1, S_err)
        hist_rawyield_cutvar[i_conf].Write()
        hist_eff_prompt_cutvar[i_conf].Write()
        hist_eff_non_prompt_cutvar[i_conf].Write()
        outfile.Close()