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
    parser.add_argument('input_file', 
                        metavar='text',
                        help='input_file_with_systematics_output')
    parser.add_argument('output_file', 
                        metavar='text',
                        default='output/sys_clened.root',
                        help='name of output file, default: output/sys_clened.root')
    parser.add_argument('syst_type', 
                        metavar='text',
                        default='yield',
                        help='type of systematic to clean (yield/mass/track)')
    args = parser.parse_args()

    outfile = ROOT.TFile(args.output_file, "RECREATE")
    infile = ROOT.TFile.Open(args.input_file)

    if args.syst_type == 'yield':
        h2_raw_yield = infile.Get("h2_raw_yield")
        h2_raw_yield.SetDirectory(0)
        h2_chi2 = infile.Get("h2_chi2")
        h2_chi2.SetDirectory(0)
        h2_sigma = infile.Get("h2_sigma")
        h2_sigma.SetDirectory(0)
        h2_signif = infile.Get("h2_signif")
        h2_signif.SetDirectory(0)


        pt_mins = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.] 
        pt_maxs = [3.0, 4.0, 5.0, 6.0, 8.0, 12., 24.]
        outfile.cd()
        ROOT.gStyle.SetOptStat(0)
        for pt_bin in range(1, h2_raw_yield.GetNbinsX() + 1):
            # default plot before skimming
            c = ROOT.TCanvas(f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_all", f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_all", 900, 600)
            c.Divide(2, 2)
            
            hists = [h2_chi2.ProjectionY(f"chi2_{pt_bin}", pt_bin, pt_bin),
                        h2_signif.ProjectionY(f"signif_{pt_bin}", pt_bin, pt_bin),
                        h2_sigma.ProjectionY(f"sigma_{pt_bin}", pt_bin, pt_bin),
                        h2_raw_yield.ProjectionY(f"raw_yield_{pt_bin}", pt_bin, pt_bin)]

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
                legend.AddEntry("", f" {pt_mins[pt_bin-1]} < pt < {pt_maxs[pt_bin - 1]} GeV/c", "")
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

            # Skimming values with large Chi2, or absurd sigmas
            i_eff = 0
            
            raw_yields = []
            raw_yields_err = [] 
            chi2s = []
            sigmas = []
            sigmas_err = []
            signifs = []

            sigma_cent = h2_sigma.GetBinContent(pt_bin, 1)
            raw_yield_cent = h2_raw_yield.GetBinContent(pt_bin, 1)

            for conf_bin in range(1, h2_raw_yield.GetNbinsY() + 1):
                raw_yield = h2_raw_yield.GetBinContent(pt_bin, conf_bin)
                raw_yield_err = h2_raw_yield.GetBinError(pt_bin, conf_bin)
                chi2 = h2_chi2.GetBinContent(pt_bin, conf_bin)
                sigma = h2_sigma.GetBinContent(pt_bin, conf_bin)
                sigma_err = h2_sigma.GetBinError(pt_bin, conf_bin)
                signif = h2_signif.GetBinContent(pt_bin, conf_bin)

                if signif < 3 or chi2 > 2 or sigma_err/sigma > 0.5 or abs(sigma - sigma_cent) > 0.2*sigma:
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
    else:
        h2_raw_yield = infile.Get("h2_raw_yield")
        h2_raw_yield.SetDirectory(0)
        h2_chi2 = infile.Get("h2_chi2")
        h2_chi2.SetDirectory(0)
        h2_sigma = infile.Get("h2_sigma")
        h2_sigma.SetDirectory(0)
        h2_signif = infile.Get("h2_signif")
        h2_signif.SetDirectory(0)
        h2_eff = infile.Get("h2_eff")
        h2_eff.SetDirectory(0)
        h2_frac = infile.Get("h2_frac")
        h2_frac.SetDirectory(0)
        h2_raw_yield = infile.Get("h2_raw_yield")
        h2_raw_yield.SetDirectory(0)
        h2_xsec = infile.Get("h2_xsec")
        h2_xsec.SetDirectory(0)


        pt_mins = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.] 
        pt_maxs = [3.0, 4.0, 5.0, 6.0, 8.0, 12., 24.]
        outfile.cd()
        ROOT.gStyle.SetOptStat(0)
        for pt_bin in range(1, h2_raw_yield.GetNbinsX() + 1):
            # default plot before skimming
            c = ROOT.TCanvas(f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_all", f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_all", 900, 600)
            c.Divide(4, 2)
            
            hists = [h2_eff.ProjectionY(f"eff_{pt_bin}", pt_bin, pt_bin),
                    h2_frac.ProjectionY(f"frac_{pt_bin}", pt_bin, pt_bin),
                    h2_raw_yield.ProjectionY(f"raw_yield_{pt_bin}", pt_bin, pt_bin),
                    h2_chi2.ProjectionY(f"chi2_{pt_bin}", pt_bin, pt_bin),
                    h2_signif.ProjectionY(f"signif_{pt_bin}", pt_bin, pt_bin),
                    h2_sigma.ProjectionY(f"sigma_{pt_bin}", pt_bin, pt_bin),
                    h2_xsec.ProjectionY(f"xsec_{pt_bin}", pt_bin, pt_bin)]
                 
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
                legend.AddEntry("", f" {pt_mins[pt_bin-1]} < pt < {pt_maxs[pt_bin - 1]} GeV/c", "")
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

            # Skimming values with large Chi2, or absurd sigmas
            i_eff = 0
            
            raw_yields = []
            raw_yields_err = [] 
            chi2s = []
            sigmas = []
            sigmas_err = []
            signifs = []
            effs = []
            effs_err = []
            fracs = []
            fracs_err = []
            xsecs = []
            xsecs_err = []

            sigma_cent = h2_sigma.GetBinContent(pt_bin, 1)
            raw_yield_cent = h2_raw_yield.GetBinContent(pt_bin, 1)
            xsec_cent = h2_xsec.GetBinContent(pt_bin, 1)

            for conf_bin in range(1, h2_raw_yield.GetNbinsY() + 1):
                raw_yield = h2_raw_yield.GetBinContent(pt_bin, conf_bin)
                raw_yield_err = h2_raw_yield.GetBinError(pt_bin, conf_bin)
                chi2 = h2_chi2.GetBinContent(pt_bin, conf_bin)
                sigma = h2_sigma.GetBinContent(pt_bin, conf_bin)
                sigma_err = h2_sigma.GetBinError(pt_bin, conf_bin)
                signif = h2_signif.GetBinContent(pt_bin, conf_bin)
                eff = h2_eff.GetBinContent(pt_bin, conf_bin)
                eff_err = h2_eff.GetBinError(pt_bin, conf_bin)
                frac = h2_frac.GetBinContent(pt_bin, conf_bin)
                frac_err = h2_frac.GetBinError(pt_bin, conf_bin)
                xsec = h2_xsec.GetBinContent(pt_bin, conf_bin)
                xsec_err = h2_xsec.GetBinError(pt_bin, conf_bin)

           

                if signif < 3 or chi2 > 2 or sigma_err/sigma > 0.5 or abs(sigma - sigma_cent) > 0.5*sigma:
                    continue
                i_eff += 1
                raw_yields. append(raw_yield)
                raw_yields_err. append(raw_yield_err)
                chi2s. append(chi2)
                sigmas. append(sigma)
                sigmas_err. append(sigma_err)
                signifs. append(signif)
                effs. append(eff)
                effs_err. append(eff_err)
                fracs. append(frac)
                fracs_err. append(frac_err)
                xsecs. append(xsec)
                xsecs_err. append(xsec_err)
            
            h_chi2_skim = ROOT.TH1F("h_chi2", "h_chi2; conf; chi2", i_eff, -0.5, i_eff -0.5)
            h_chi2_skim.SetMarkerStyle(20)
            h_signif_skim = ROOT.TH1F("h_signif", "h_signif; conf; signif", i_eff, -0.5, i_eff -0.5)
            h_signif_skim.SetMarkerStyle(20)
            h_sigma_skim = ROOT.TH1F("h_sigma", "h_sigma; conf; sigma", i_eff, -0.5, i_eff -0.5)
            h_sigma_skim.SetMarkerStyle(20)
            h_raw_yield_skim = ROOT.TH1F("h_raw_yield", "h_raw_yield; conf; raw_yield", i_eff, -0.5, i_eff -0.5)
            h_raw_yield_skim.SetMarkerStyle(20)
            h_eff_skim = ROOT.TH1F("h_eff", "h_eff; conf; eff", i_eff, -0.5, i_eff -0.5)
            h_eff_skim.SetMarkerStyle(20)
            h_frac_skim = ROOT.TH1F("h_frac", "h_frac; conf; frac", i_eff, -0.5, i_eff -0.5)
            h_frac_skim.SetMarkerStyle(20)
            h_xsec_skim = ROOT.TH1F("h_xsec", "h_xsec; conf; xsec", i_eff, -0.5, i_eff -0.5)
            h_xsec_skim.SetMarkerStyle(20)
            
            
            for ibin, (raw_yield, raw_yield_err, chi2, sigma, sigma_err, signif, eff, eff_err, frac, frac_err, xsec, xsec_err) in enumerate(zip (raw_yields, raw_yields_err, chi2s, sigmas, sigmas_err, signifs, effs, effs_err, fracs, fracs_err, xsecs, xsecs_err)):
                h_chi2_skim.SetBinContent(ibin + 1, chi2)
                h_chi2_skim.SetBinError(ibin + 1, chi2/100)
                h_signif_skim.SetBinContent(ibin + 1, signif)
                h_signif_skim.SetBinError(ibin + 1, signif/100)
                h_sigma_skim.SetBinContent(ibin + 1, sigma)
                h_sigma_skim.SetBinError(ibin + 1, sigma_err)
                h_raw_yield_skim.SetBinContent(ibin + 1, raw_yield)
                h_raw_yield_skim.SetBinError(ibin + 1, raw_yield_err)
                h_frac_skim.SetBinContent(ibin + 1, frac)
                h_frac_skim.SetBinError(ibin + 1, frac_err)
                h_eff_skim.SetBinContent(ibin + 1, eff)
                h_eff_skim.SetBinError(ibin + 1, eff_err)
                h_xsec_skim.SetBinContent(ibin + 1, xsec)
                h_xsec_skim.SetBinError(ibin + 1, xsec_err)
            
            c2 = ROOT.TCanvas(f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_skim", f"c_{pt_mins[pt_bin-1]}_{pt_maxs[pt_bin - 1]}_skim", 900, 600)
            c2.Divide(2, 2)
            c2.cd(3)
            h_signif_skim.Draw("PE")
            c2.cd(2)
            h_raw_yield_skim.Draw("PE")
            c2.cd(1)
            h_eff_skim.Draw("PE")
            c2.cd(4)
            h_xsec_skim.Draw("PE")
            
            line2 = ROOT.TLine(h_xsec_skim.GetXaxis().GetXmin(), xsec_cent, h_xsec_skim.GetXaxis().GetXmax(), xsec_cent)
            line2.SetLineColor(ROOT.kRed)
            line2.Draw("same")
            
            legend2 = ROOT.TLegend(0.5, 0.7, 0.8, 0.9)
            legend2.SetBorderSize(0)
            legend2.SetFillStyle(0)
            legend2.SetTextSize(0.03)
            legend2.AddEntry("", f" {pt_mins[pt_bin-1]} < pt < {pt_maxs[pt_bin - 1]} GeV/c", "")
            legend2.AddEntry(line2, "Central value", "L")
            mean = np.mean(xsecs)
            rms = np.std(xsecs)
            delta = np.abs(xsec_cent - mean)
            relative_delta = delta / xsec_cent
            relative_RMS = rms / xsec_cent
            total_uncertainty = np.sqrt(delta**2 + rms**2)
            total_relative_uncertainty = total_uncertainty / xsec_cent
            band_rms2 = ROOT.TBox(hist.GetXaxis().GetXmin(), xsec_cent - rms, hist.GetXaxis().GetXmax(), xsec_cent + rms)
            band_rms2.SetFillColorAlpha(ROOT.kMagenta, 0.3)
            band_rms2.SetLineColor(ROOT.kMagenta)
            band_rms2.SetLineStyle(2)
            band_rms2.Draw("same")

            band_total_uncertainty2 = ROOT.TBox(hist.GetXaxis().GetXmin(), xsec_cent - total_uncertainty, hist.GetXaxis().GetXmax(), xsec_cent + total_uncertainty)
            band_total_uncertainty2.SetFillColorAlpha(ROOT.kRed, 0.3)
            band_total_uncertainty2.SetLineColor(ROOT.kRed)
            band_total_uncertainty2.SetLineStyle(2)
            band_total_uncertainty2.Draw("same")
            legend2.AddEntry(band_rms2, f"RMS: {(relative_RMS*100):.1f}%", "LF")
            legend2.AddEntry(band_total_uncertainty2, f"RMS + shift: {(total_relative_uncertainty*100):.1f}% ", "LF")
            legend2.Draw()
            c2.Write()

    outfile.Close()
    