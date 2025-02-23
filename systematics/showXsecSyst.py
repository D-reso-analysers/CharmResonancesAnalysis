import argparse
import os
import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('input_file_mb', 
                        metavar='text',
                        help='input_file_with_systematics_output MB')
    parser.add_argument('input_file_trg', 
                        metavar='text',
                        help='input_file_with_systematics_output TRG')
    args = parser.parse_args()

    infile_mb = ROOT.TFile.Open(args.input_file_mb)
    infile_trg = ROOT.TFile.Open(args.input_file_trg)


    h2_xsec_mb = infile_mb.Get("h2_xsec")
    h2_xsec_mb.SetDirectory(0)
    h2_xsec_trg = infile_trg.Get("h2_xsec")
    h2_xsec_trg.SetDirectory(0)
    
    pt_mins = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.] 
    pt_maxs = [3.0, 4.0, 5.0, 6.0, 8.0, 12., 24.]
    
    xsec_mb_var
    xsec_trg_var
    for conf_bin in range(1, h2_raw_yield.GetNbinsY() + 1):
        xsec_conf = []
        xsec_err_conf = []
        pt_conf = []
        for pt_bin in range(1, h2_xsec_mb.GetNbinsX() + 1):
            pt_bin_center = h2_xsec_mb.GetXaxis().GetBinCenter(pt_bin)
            var_ipt_mb = []
            var_ipt_trg = []
        
            xsec_mb = h2_xsec_mb.GetBinContent(pt_bin, conf_bin)
            xsec_mb_err = h2_xsec_mb.GetBinError(pt_bin, conf_bin)
            xsec_trg = h2_xsec_trg.GetBinContent(pt_bin, conf_bin)
            xsec_trg_err = h2_xsec_trg.GetBinError(pt_bin, conf_bin)
            
            if signif < 3 or chi2 > 2 or sigma_err/sigma > 0.5 or abs(sigma - sigma_cent) > sigma:
                continue
            i_eff += 1
            raw_yields. append(raw_yield)
            raw_yields_err. append(raw_yield_err)
            chi2s. append(chi2)
            sigmas. append(sigma)
            sigmas_err. append(sigma_err)
            signifs. append(signif)
        
        h_chi2_skim = ROOT.TH1F("h_chi2", "h_chi2; conf; chi2", i_eff, -0.5, i_eff -0.5)
        h_chi2_skim.SetMarkerStyle(20)
        h_signif_skim = ROOT.TH1F("h_signif", "h_signif; conf; signif", i_eff, -0.5, i_eff -0.5)
        h_signif_skim.SetMarkerStyle(20)
        h_sigma_skim = ROOT.TH1F("h_sigma", "h_sigma; conf; sigma", i_eff, -0.5, i_eff -0.5)
        h_sigma_skim.SetMarkerStyle(20)
        h_raw_yield_skim = ROOT.TH1F("h_raw_yield", "h_raw_yield; conf; raw_yield", i_eff, -0.5, i_eff -0.5)
        h_raw_yield_skim.SetMarkerStyle(20)
        
        for ibin, (raw_yield, raw_yield_err, chi2, sigma, sigma_err, signif) in enumerate(zip (raw_yields, raw_yields_err, chi2s, sigmas, sigmas_err, signifs)):
            h_chi2_skim.SetBinContent(ibin + 1, chi2)
            h_chi2_skim.SetBinError(ibin + 1, chi2/100)
            h_signif_skim.SetBinContent(ibin + 1, signif)
            h_signif_skim.SetBinError(ibin + 1, signif/100)
            h_sigma_skim.SetBinContent(ibin + 1, sigma)
            h_sigma_skim.SetBinError(ibin + 1, sigma_err)
            h_raw_yield_skim.SetBinContent(ibin + 1, raw_yield)
            h_raw_yield_skim.SetBinError(ibin + 1, raw_yield_err)
        
        c2 = ROOT.TCanvas(f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_skim", f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_skim", 900, 600)
        c2.Divide(2, 2)
        c2.cd(1)
        h_chi2_skim.Draw("PE")
        c2.cd(2)
        h_signif_skim.Draw("PE")
        c2.cd(3)
        h_sigma_skim.Draw("PE")
        c2.cd(4)
        h_raw_yield_skim.Draw("PE")
        line2 = ROOT.TLine(h_raw_yield_skim.GetXaxis().GetXmin(), raw_yield_cent, h_raw_yield_skim.GetXaxis().GetXmax(), raw_yield_cent)
        line2.SetLineColor(ROOT.kRed)
        line2.Draw("same")
        
        legend2 = ROOT.TLegend(0.5, 0.7, 0.8, 0.9)
        legend2.SetBorderSize(0)
        legend2.SetFillStyle(0)
        legend2.SetTextSize(0.03)
        legend2.AddEntry("", f" {pt_mins[pt_bin-1]} < pt < {pt_maxs[pt_bin - 1]} GeV/c", "")
        legend2.AddEntry(line2, "Central value", "L")
        mean = np.mean(raw_yields)
        rms = np.std(raw_yields)
        delta = np.abs(raw_yield_cent - mean)
        relative_delta = delta / raw_yield_cent
        relative_RMS = rms / raw_yield_cent
        total_uncertainty = np.sqrt(delta**2 + rms**2)
        total_relative_uncertainty = total_uncertainty / raw_yield_cent
        band_rms2 = ROOT.TBox(hist.GetXaxis().GetXmin(), raw_yield_cent - rms, hist.GetXaxis().GetXmax(), raw_yield_cent + rms)
        band_rms2.SetFillColorAlpha(ROOT.kMagenta, 0.3)
        band_rms2.SetLineColor(ROOT.kMagenta)
        band_rms2.SetLineStyle(2)
        band_rms2.Draw("same")

        band_total_uncertainty2 = ROOT.TBox(hist.GetXaxis().GetXmin(), raw_yield_cent - total_uncertainty, hist.GetXaxis().GetXmax(), raw_yield_cent + total_uncertainty)
        band_total_uncertainty2.SetFillColorAlpha(ROOT.kRed, 0.3)
        band_total_uncertainty2.SetLineColor(ROOT.kRed)
        band_total_uncertainty2.SetLineStyle(2)
        band_total_uncertainty2.Draw("same")
        legend2.AddEntry(band_rms2, f"RMS: {(relative_RMS*100):.1f}%", "LF")
        legend2.AddEntry(band_total_uncertainty2, f"RMS + shift: {(total_relative_uncertainty*100):.1f}% ", "LF")
        legend2.Draw()
        c2.Write()

    outfile.Close()
    