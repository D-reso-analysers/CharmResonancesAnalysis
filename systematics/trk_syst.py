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
sys.path.insert(0, '..')
from utils.analysis_utils import loadAO2D, applySelections, applySelections2, perform_roofit_fit, get_chi2_significance_sb, set_param

br_dplus = 0.0938
br_dstar = 0.677 * 0.03951
br_k0 = 0.692
br_ds1 = 0.22
br_ds2star = 0.2335

def compute_efficiency(nRec, nGen):
    eff = nRec / nGen
    effErr = eff * np.sqrt((1 - eff) / nRec)
    return eff, effErr

def compute_cross_section (raw_yield, eff, frac, dpt, br, lumi):
    """
    Compute the cross section from the raw yield, efficiency and prompt fraction.

    Args:
        raw_yield (float): Raw yield.
        eff (float): Efficiency.
        frac (float): Fraction.
        dpt (float): Bin width.
        br (float): Branching ratio.
        lumi (float): Luminosity.

    Returns:
        float: Cross section.
    """
    return 1./2 * frac * raw_yield / eff / dpt / br / lumi

def trackingSelections(df, itsNCls0, itsNCls1, tpcRows0, tpcRows1, tpcChi0, tpcChi1):
    dfOut = pd.DataFrame()
    dfOut = df[(df['fItsNClsProngMinBach0'] >= itsNCls0) &
               (df['fItsNClsProngMinBach1'] >= itsNCls1) &
               (df['fTpcNClsCrossedRowsProngMinBach0'] >= tpcRows0) &
               (df['fTpcNClsCrossedRowsProngMinBach1'] >= tpcRows1) &
               (df['fTpcChi2NClProngMaxBach0'] <= tpcChi0) &
               (df['fTpcChi2NClProngMaxBach1'] <= tpcChi1)]
    return dfOut

# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('configfile', 
                        metavar='text',
                        help='config file name with inputs and cut configurations')
    args = parser.parse_args()
    with open(args.configfile, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    # initialize output files
    if cfg['tracking']['optionalRawYieldOutput']:
        rawYieldOutputFile = ROOT.TFile(cfg['tracking']['optionalRawYieldOutput'], "RECREATE")
    else:
        rawYieldOutputFile = None
    outfile = ROOT.TFile(cfg['tracking']['outputFile'], "RECREATE")

     # Load and skim MC dataframes
    mcFileNames = cfg['fileNameMC']
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

    # Load luminosity form input file
    luminosityFile = ROOT.TFile(cfg['luminosityFile'])
    luminosityHist = luminosityFile.Get(cfg['luminosityHist'])
    luminosity = luminosityHist.GetBinContent(1)

    # Load branching ratio
    if cfg['particle'] == "ds1":
        br = br_ds1 * br_dstar * br_k0
    elif cfg['particle'] == "ds2star":
        br = br_ds2star * br_dplus * br_k0

    # Load data dataframe from file, and apply preliminary selections
    dataFileName = cfg['fileNameData']
    inTreeNameData = cfg['treeNameData']
    dfData = loadAO2D(dataFileName, inTreeNameData)
    dfDataFiltered = applySelections(dfData, cfg, isMC=False)

    # apply central cut selections on BdtScore
    centralCuts = cfg['centralCuts']
    flagsToKeep = cfg['acceptFlags']
    configurations = cfg['bdtSelections']['configurations']
    pt_bins_dau_min = cfg['cutVars']['ptBach0']['min']
    pt_bins_dau_max = cfg['cutVars']['ptBach0']['max']
    for conf in configurations:
        if conf['name'] == 'central':
            bdt_cuts = conf['cuts']
            frac_file = conf['frac_file']
            frac_hist = conf['frac_hist']
            break
    dfDataCentralCut = applySelections2(dfDataFiltered, centralCuts)
    dfMcRecCentralCut = [applySelections2(dfRec, centralCuts, isMC=True, selectedFlags=flagsToKeep) for irec, dfRec in enumerate(dfRecList)]
    # Load specific configurations andi initialize output histograms

    # retrieve configurations for yield extraction
    pt_mins = cfg['cutVars']['pt']['min']
    pt_maxs = cfg['cutVars']['pt']['max']
    
    pt_limits = pt_mins.copy()
    pt_limits.append(pt_maxs[-1])
    pt_limits = np.array(pt_limits, np.float64)

    # fit configuration
    fit_params_def = cfg['rawYield']['default_fit_config']['parameters']
    signal_pdf = cfg['rawYield']['default_fit_config']['signal_pdf']
    bkg_pdf = cfg['rawYield']['default_fit_config']['bkg_pdf']
    mass_min = cfg['rawYield']['default_fit_config']['mass_min']
    mass_max = cfg['rawYield']['default_fit_config']['mass_max']

    variations = cfg['tracking']['variations']
    itsNCls0s = variations['itsNCls0']
    itsNCls1s = variations['itsNCls1']
    tpcRows0s = variations['tpcRows0']
    tpcRows1s = variations['tpcRows1']
    tpcChi0s = variations['tpcChi0']
    tpcChi1s = variations['tpcChi1']

    nvars = len(itsNCls0s)*len(itsNCls1s)*len(tpcRows0s)*len(tpcRows1s)*len(tpcChi0s)*len(tpcChi1s)

    h2_raw_yield = ROOT.TH2F("h2_raw_yield", "raw yield; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_chi2 = ROOT.TH2F("h2_chi2", "chi2; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_sigma = ROOT.TH2F("h2_sigma", "sigma; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_signif = ROOT.TH2F("h2_signif", "significance; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_eff = ROOT.TH2F("h2_eff", "efficiency; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_frac = ROOT.TH2F("h2_frac", "prompt fraction; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_sb = ROOT.TH2F("h2_sb", "signal over background; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_xsec = ROOT.TH2F("h2_xsec", "cross section; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    
    central_mean = 0
    central_sigma = 0
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        # split dataframes in pt bins
        dfDataPt = dfDataCentralCut[(dfDataCentralCut['fPt'] >= pt_min) & (dfDataCentralCut['fPt'] < pt_max)]
        mcRecPtList = [dfRec[(dfRec['fPt'] >= pt_min) & (dfRec['fPt'] < pt_max)] for dfRec in dfMcRecCentralCut]
        mcGenPtList = [dfGen[(dfGen['fPt'] >= pt_min) & (dfGen['fPt'] < pt_max)] for dfGen in dfGenList]
        i_conf = 0
        # nested loop on variations
        for itsNCls0 in itsNCls0s:
            for itsNCls1 in itsNCls1s:
                for tpcRows0 in tpcRows0s:
                    for tpcRows1 in tpcRows1s:
                        for tpcChi0 in tpcChi0s:
                            for tpcChi1 in tpcChi1s:
                                # 0) load fraction from file
                                frac_file = ROOT.TFile.Open(conf['frac_file'])
                                frac_hist = frac_file.Get(conf['frac_hist'])
                                frac_hist.SetDirectory(0)
                                frac_file.Close()
                                # 1) specific selections
                                dfDataSelPt = trackingSelections(dfDataPt, itsNCls0, itsNCls1, tpcRows0, tpcRows1, tpcChi0, tpcChi1)
                                mcRecSelPtList = [trackingSelections(dfRec, itsNCls0, itsNCls1, tpcRows0, tpcRows1, tpcChi0, tpcChi1) for dfRec in mcRecPtList]
                                # 2) efficiency
                                nRecPrompt, nGenPrompt = 0, 0
                                for i, (dfRec, dfGen, w) in enumerate(zip(mcRecSelPtList, mcGenPtList, cfg['mcWeights'])):
                                    nGenPrompt += w * len(dfGen[dfGen['fOrigin'] == 1])
                                    nRecPrompt += w * len(dfRec[(dfRec['fOrigin'] == 1)])
                            
                                nRecPrompt = nRecPrompt / np.sum(cfg['mcWeights'])
                                nGenPrompt = nGenPrompt / np.sum(cfg['mcWeights'])
                                eff, effErr = compute_efficiency(nRecPrompt, nGenPrompt)

                                h2_eff.SetBinContent(ipt + 1, i_conf + 1, eff)
                                h2_eff.SetBinError(ipt + 1, i_conf + 1, effErr)

                                # 3)raw yield
                                fit_params = copy.deepcopy(fit_params_def)
                                if i_conf > 0:
                                    set_param(fit_params, "sigma", central_sigma, central_sigma, central_sigma)
                                    set_param(fit_params, "mean", central_mean, central_mean, central_mean)

                                ws = perform_roofit_fit(dfDataSelPt, signal_pdf, bkg_pdf,  fit_params, mass_min, mass_max)
                                chi2, significance, sb = get_chi2_significance_sb(ws, 100, outfile = rawYieldOutputFile, name = f"c_{pt_min}_{pt_max}_{i_conf}")

                                S = ws.var("nSig").getVal()
                                S_err = ws.var("nSig").getError()
                                sigma = ws.var("sigma").getVal()
                                sigma_err = ws.var("sigma").getError()
                                mean = ws.var("mean").getVal()
                                mean_err = ws.var("mean").getError()
                                if i_conf == 0:
                                    central_mean = mean
                                    central_sigma = sigma

                                h2_raw_yield.SetBinContent(ipt + 1, i_conf + 1, S)
                                h2_raw_yield.SetBinError(ipt + 1, i_conf + 1, S_err)
                                h2_chi2.SetBinContent(ipt + 1, i_conf + 1, chi2)
                                h2_chi2.SetBinError(ipt + 1, i_conf + 1, chi2/100)
                                h2_signif.SetBinContent(ipt + 1, i_conf + 1, significance)
                                h2_signif.SetBinError(ipt + 1, i_conf + 1, significance/100)
                                h2_sb.SetBinContent(ipt + 1, i_conf + 1, sb)
                                h2_sb.SetBinError(ipt + 1, i_conf + 1, sb/100)
                                h2_sigma.SetBinContent(ipt + 1, i_conf + 1, sigma)
                                h2_sigma.SetBinError(ipt + 1, i_conf + 1, sigma_err)

                                # 4) cross section
                                dpt = pt_max - pt_min
                                xsec = compute_cross_section(S, eff, frac_hist.GetBinContent(ipt + 1), dpt, br, luminosity)
                                xsec_err = xsec * S_err / S

                                h2_xsec.SetBinContent(ipt + 1, i_conf + 1, xsec)
                                h2_xsec.SetBinError(ipt + 1, i_conf + 1, xsec_err)

                                i_conf+=1




    # Loop over pt bins and create projections
    outfile.cd()
    ROOT.gStyle.SetOptStat(0)
    for ipt in range(1, len(pt_limits)):
        c = ROOT.TCanvas(f"c_{pt_mins[ipt-1]}_{pt_maxs[ipt - 1]}", f"c_{pt_mins[ipt-1]}_{pt_maxs[ipt - 1]}", 900, 600)
        c.Divide(4, 2)
        

        hists = [h2_eff.ProjectionY(f"eff_{ipt}", ipt, ipt),
                 h2_frac.ProjectionY(f"frac_{ipt}", ipt, ipt),
                 h2_raw_yield.ProjectionY(f"raw_yield_{ipt}", ipt, ipt),
                 h2_chi2.ProjectionY(f"chi2_{ipt}", ipt, ipt),
                 h2_signif.ProjectionY(f"signif_{ipt}", ipt, ipt),
                 h2_sb.ProjectionY(f"sb_{ipt}", ipt, ipt),
                 h2_sigma.ProjectionY(f"sigma_{ipt}", ipt, ipt),
                 h2_xsec.ProjectionY(f"xsec_{ipt}", ipt, ipt)]

        for i, hist in enumerate(hists):
            c.cd(i + 1)
            hist.SetMarkerStyle(20)
            hist.Draw("PE")
            if hist.GetEntries() > 0:
                central_value = hist.GetBinContent(1)
                line = ROOT.TLine(hist.GetXaxis().GetXmin(), central_value, hist.GetXaxis().GetXmax(), central_value)
                line.SetLineColor(ROOT.kRed)
                line.Draw("same")
                
                legend = ROOT.TLegend(0.5, 0.7, 0.8, 0.9)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)
                legend.AddEntry("", f" {pt_mins[ipt-1]} < pt < {pt_maxs[ipt - 1]} GeV/c", "")
                legend.AddEntry(line, "Central value", "L")
                if "xsec" in hist.GetName():
                    bin_contents = [hist.GetBinContent(bin) for bin in range(2, hist.GetNbinsX() + 1)]
                    mean = np.mean(bin_contents)
                    rms = np.std(bin_contents)
                    delta = np.abs(central_value - mean)
                    relative_delta = delta / central_value
                    relative_RMS = rms / central_value
                    total_uncertainty = np.sqrt(delta**2 + rms**2)
                    total_relative_uncertainty = total_uncertainty / central_value
                    band_rms = ROOT.TBox(hist.GetXaxis().GetXmin(), central_value - rms, hist.GetXaxis().GetXmax(), central_value + rms)
                    band_rms.SetFillColorAlpha(ROOT.kMagenta, 0.3)
                    band_rms.SetLineColor(ROOT.kMagenta)
                    band_rms.SetLineStyle(2)
                    band_rms.Draw("same")

                    band_total_uncertainty = ROOT.TBox(hist.GetXaxis().GetXmin(), central_value - total_uncertainty, hist.GetXaxis().GetXmax(), central_value + total_uncertainty)
                    band_total_uncertainty.SetFillColorAlpha(ROOT.kRed, 0.3)
                    band_total_uncertainty.SetLineColor(ROOT.kRed)
                    band_total_uncertainty.SetLineStyle(2)
                    band_total_uncertainty.Draw("same")
                    legend.AddEntry(band_rms, f"RMS: {(relative_RMS*100):.1f}%", "LF")
                    legend.AddEntry(band_total_uncertainty, f"RMS + shift: {(total_relative_uncertainty*100):.1f}% ", "LF")
                legend.Draw()
        c.Write()

    h2_eff.Write()
    h2_frac.Write()
    h2_raw_yield.Write()
    h2_chi2.Write()
    h2_signif.Write()
    h2_sb.Write()
    h2_sigma.Write()
    h2_xsec.Write()
    outfile.Close()

    if rawYieldOutputFile:
        rawYieldOutputFile.Close()