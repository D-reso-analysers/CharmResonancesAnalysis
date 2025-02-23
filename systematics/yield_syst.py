import argparse
import os
import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd
import sys
from particle import Particle
sys.path.insert(0, '..')
from utils.analysis_utils import loadAO2D, applySelections, applySelections2, perform_roofit_fit, get_chi2_significance_sb, set_param

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
    if cfg['rawYield']['optionalRawYieldOutput']:
        rawYieldOutputFile = ROOT.TFile(cfg['rawYield']['optionalRawYieldOutput'], "RECREATE")
    else:
        rawYieldOutputFile = None
    outfile = ROOT.TFile(cfg['rawYield']['outputFile'], "RECREATE")

    # Load data dataframe from file, and apply preliminary selections
    dataFileName = cfg['fileNameData']
    inTreeNameData = cfg['treeNameData']
    dfData = loadAO2D(dataFileName, inTreeNameData)
    dfDataFiltered = applySelections(dfData, cfg, isMC=False)

    # apply central cut selections on BdtScore
    centralCuts = cfg['centralCuts']
    dfDataCentralCut = applySelections2(dfDataFiltered, centralCuts)
    dfDataCentralCut = dfDataCentralCut[['fM','fPt','fPtBach0', 'fMlScoreBkgBach0']]

    # retrieve configurations for yield extraction
    pt_mins = cfg['cutVars']['pt']['min']
    pt_maxs = cfg['cutVars']['pt']['max']
    
    pt_limits = pt_mins.copy()
    pt_limits.append(pt_maxs[-1])
    pt_limits = np.array(pt_limits, np.float64)

    fit_params_def = cfg['rawYield']['default_fit_config']['parameters']
    signal_pdf = cfg['rawYield']['default_fit_config']['signal_pdf']
    
    variations = cfg['rawYield']['variations']
    mass_mins = variations['mass_min']
    mass_maxs = variations['mass_max']
    bkg_pdfs = variations["background_pdf"]
    widths = variations['width']

    nvars = len(mass_mins)*len(mass_maxs)*len(bkg_pdfs)*len(widths)

    h2_raw_yield = ROOT.TH2F("h2_raw_yield", "raw yield; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_chi2 = ROOT.TH2F("h2_chi2", "chi2; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_sigma = ROOT.TH2F("h2_sigma", "sigma; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)
    h2_signif = ROOT.TH2F("h2_signif", "significance; pT[GeV/c]; cutset", len(pt_limits) - 1, pt_limits, nvars, -0.5, nvars - 0.5)

    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        # split dataframes in pt bins
        dfDataPt = dfDataCentralCut[(dfDataCentralCut['fPt'] >= pt_min) & (dfDataCentralCut['fPt'] < pt_max)]
        i_conf = 0
        # nested loop on all variations
        for mass_min in mass_mins:
            for mass_max in mass_maxs:
                for bkg_pdf in bkg_pdfs:
                    for width in widths:
                        # set parameters specific for this variation
                        fit_params = fit_params_def
                        set_param(fit_params, "width", width, width, width)
                        conf_name = f"{bkg_pdf}_w{width}_{mass_min}-{mass_max}"
                        ws = perform_roofit_fit(dfDataPt, signal_pdf, bkg_pdf,  fit_params, mass_min, mass_max)
                        chi2, significance, sb = get_chi2_significance_sb(ws, 100, outfile = rawYieldOutputFile, name = f"{pt_min}_{i_conf}_{pt_max}_{conf_name}")

                        s = ws.var("nSig").getVal()
                        s_err = ws.var("nSig").getError()
                        sigma = ws.var("sigma").getVal()
                        sigma_err = ws.var("sigma").getError()

                        h2_raw_yield.SetBinContent(ipt + 1, i_conf + 1, s)
                        h2_raw_yield.SetBinError(ipt + 1, i_conf + 1, s_err)
                        h2_chi2.SetBinContent(ipt + 1, i_conf + 1, chi2)
                        h2_chi2.SetBinError(ipt + 1, i_conf + 1, chi2/100)
                        h2_signif.SetBinContent(ipt + 1, i_conf + 1, significance)
                        h2_signif.SetBinError(ipt + 1, i_conf + 1, significance/100)
                        h2_sigma.SetBinContent(ipt + 1, i_conf + 1, sigma)
                        h2_sigma.SetBinError(ipt + 1, i_conf + 1, sigma_err)

                        i_conf += 1

    # Loop over pt bins and create projections
    outfile.cd()
    ROOT.gStyle.SetOptStat(0)
    for ipt in range(1, len(pt_limits)):
        c = ROOT.TCanvas(f"c_{pt_mins[ipt-1]}_{pt_maxs[ipt - 1]}", f"c_{pt_mins[ipt-1]}_{pt_maxs[ipt - 1]}", 900, 600)
        c.Divide(2, 2)
        
        hists = [h2_chi2.ProjectionY(f"chi2_{ipt}", ipt, ipt),
                 h2_signif.ProjectionY(f"signif_{ipt}", ipt, ipt),
                 h2_sigma.ProjectionY(f"sigma_{ipt}", ipt, ipt),
                 h2_raw_yield.ProjectionY(f"raw_yield_{ipt}", ipt, ipt)]

        for i, hist in enumerate(hists):
            c.cd(i + 1)
            hist.SetMarkerStyle(20)
            hist.Draw("PE")
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
            if "raw_yield" in hist.GetName():
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
    # Save TH2 in output file
    h2_raw_yield.Write()
    h2_chi2.Write()
    h2_signif.Write()
    h2_sigma.Write()
    outfile.Close()

    if rawYieldOutputFile:
        rawYieldOutputFile.Close()
    