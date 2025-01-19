import ROOT
import argparse
import numpy as np
import array
from utils.style_formatter import set_object_style
import sys

def rebin_to_common(hist, common_bins, sum_errors=False, sum_in_quadrature = False):
    '''
    Function to rebin a histogram to a given set of bin edges
    '''
    bin_edges = []
    for i in range(1, hist.GetNbinsX() + 2):  # +2 to include the upper edge of the last bin
        bin_edges.append(hist.GetBinLowEdge(i))
    if bin_edges == common_bins:
        print("The histograms already have the same binning")
        return hist
    ibin_new = 0
    new_bin_lims = []
    nbins = hist.GetNbinsX()
    rebinned_hist = ROOT.TH1F(
        f"{hist.GetName()}_rebinned", f"{hist.GetTitle()};{hist.GetXaxis().GetTitle()};{hist.GetYaxis().GetTitle()}", 
        len(common_bins) - 1, np.array(common_bins, np.float64))
    for ibin in range(1, nbins+2):
        low_edge = hist.GetBinLowEdge(ibin)
        if low_edge not in common_bins or low_edge in new_bin_lims:
            continue
        print(f"found low_edge: {low_edge}")
        new_bin_lims.append(low_edge)
        ibin_new += 1
        for ibin_high in range(ibin + 1, nbins + 2):
            high_edge = hist.GetBinLowEdge(ibin_high)
            if high_edge in common_bins:
                print(f"found high_edge: {hist.GetBinLowEdge(ibin_high)}")
                print(f'ibin: {ibin}, ibin_high: {ibin_high}, ibin_new: {ibin_new}')
                y, sy = 0, 0
                for ibins_merged in range (ibin, ibin_high):
                    dx = hist.GetBinLowEdge(ibins_merged + 1) - hist.GetBinLowEdge(ibins_merged)
                    if sum_in_quadrature:
                        y += (dx * hist.GetBinContent(ibins_merged))**2
                    else:
                        y += dx * hist.GetBinContent(ibins_merged)
                    print(f'content of bin {ibins_merged}: {hist.GetBinContent(ibins_merged)}')
                    if sum_errors: 
                        sy += (dx * hist.GetBinError(ibins_merged))**2
                new_width = high_edge - low_edge
                if sum_in_quadrature:
                    y = np.sqrt(y)
                rebinned_hist.SetBinContent(ibin_new, y/new_width)
                print(f'Filling bin {ibin_new}, from {low_edge} to {high_edge}')
                print(f'New content {y/new_width}, new error {sy}')
                if sum_errors:
                    rebinned_hist.SetBinError(ibin_new, np.sqrt(sy)/new_width)
                break
    return rebinned_hist

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_reso", "-ire", metavar="text",
                        default="xsec_Ds1_trg_syst.root",
                        help="input file with resonance cross section")
    parser.add_argument("--input_ds", "-ids", metavar="text",
                        default="HEPData-ins2697877-v1-Figure_4_Ds_pp_13_TeV.root",
                        help="Hepdata file with Ds cross section")
    parser.add_argument("--input_run2", "-ir2", metavar="text",
                        default=None, 
                        help="input file with run2 ratios")
    parser.add_argument("--particle", "-p", metavar="text",
                        default="ds1",
                        help="particle to be analyzed")
    args = parser.parse_args()

    if args.particle == "ds1":
        br = 0.22
        flag = "D_{s1}^{+}"
    elif args.particle == "ds2star":
        br = 0.2335
        flag = "D_{s2}^{*+}"
    else:
        print(f"ERROR: particle {args.particle} not supported, exit")
        sys.exit()
    data_reso = ROOT.TFile.Open(args.input_reso)
    data_ds = ROOT.TFile.Open(args.input_ds)

    h_reso_stat = data_reso.Get('hist_xsec_stat')
    h_reso_syst = data_reso.Get('hist_xsec_syst')

    h_ds_y = data_ds.Get('Figure 4 Ds pp 13 TeV/Hist1D_y1')
    h_ds_stat = data_ds.Get('Figure 4 Ds pp 13 TeV/Hist1D_y1_e1')
    fd_ds_syst = data_ds.Get('Figure 4 Ds pp 13 TeV/Hist1D_y1_e2plus')

    bin_edges_reso, bin_edges_ds = [], []
    n_bins_reso = h_reso_stat.GetNbinsX()
    for i in range(1, n_bins_reso + 2):  # +2 to include the upper edge of the last bin
        bin_edges_reso.append(h_reso_stat.GetBinLowEdge(i))
    n_bins_ds = h_ds_y.GetNbinsX()
    for i in range(1, n_bins_ds + 2):  # +2 to include the upper edge of the last bin
        bin_edges_ds.append(h_ds_y.GetBinLowEdge(i))
    
    # Finiding common bin edges
    common_bins = list(set(bin_edges_reso).intersection(set(bin_edges_ds)))
    
    h_reso_stat_rebinned = rebin_to_common(h_reso_stat, common_bins, sum_errors=True)
    h_reso_syst_rebinned = rebin_to_common(h_reso_syst, common_bins, sum_errors=True)
    h_ds_y_rebinned = rebin_to_common(h_ds_y, common_bins)
    h_ds_stat_rebinned = rebin_to_common(h_ds_stat, common_bins, sum_in_quadrature=True)
    h_ds_syst_rebinned = rebin_to_common(fd_ds_syst, common_bins, sum_in_quadrature=True)

    # Make true ratio:
    h_ratio_stat = ROOT.TH1F(
    f"h_ratio_stat", f"{flag}/Ds ratio;{h_reso_stat_rebinned.GetXaxis().GetTitle()};{flag}/Ds", 
    len(common_bins) - 1, np.array(common_bins, np.float64))
    h_ratio_syst = ROOT.TH1F(
    f"h_ratio_syst", f"{flag}/Ds ratio;{h_reso_stat_rebinned.GetXaxis().GetTitle()};{flag}/Ds", 
    len(common_bins) - 1, np.array(common_bins, np.float64))

    for ibin in range(1, h_reso_stat_rebinned.GetNbinsX() + 1):
        r = h_reso_stat_rebinned.GetBinContent(ibin) / h_ds_y_rebinned.GetBinContent(ibin)
        r_stat = r * ((h_reso_stat_rebinned.GetBinError(ibin) / h_reso_stat_rebinned.GetBinContent(ibin))**2 + (h_ds_stat_rebinned.GetBinContent(ibin) / h_ds_y_rebinned.GetBinContent(ibin))**2)**0.5
        r_syst = r * ((h_reso_syst_rebinned.GetBinError(ibin) / h_reso_syst_rebinned.GetBinContent(ibin))**2 + (h_ds_syst_rebinned.GetBinContent(ibin) / h_ds_y_rebinned.GetBinContent(ibin))**2)**0.5
        # r_stat = r * ((h_reso_stat_rebinned.GetBinError(ibin) / h_reso_stat_rebinned.GetBinContent(ibin))**2 )**0.5
        # r_syst = r * ((h_reso_syst_rebinned.GetBinError(ibin) / h_reso_syst_rebinned.GetBinContent(ibin))**2 )**0.5


        h_ratio_stat.SetBinContent(ibin, r)
        h_ratio_stat.SetBinError(ibin, r_stat)
        h_ratio_syst.SetBinContent(ibin, r)
        h_ratio_syst.SetBinError(ibin, r_syst)
    
    outfile = ROOT.TFile.Open(f"ratio_{args.particle}.root", "RECREATE")
    h_ratio_stat.Write()
    h_ratio_syst.Write()

    # make plot 
    ROOT.gStyle.SetOptStat(0)
    set_object_style(h_ratio_stat, color=ROOT.kBlack)
    set_object_style(h_ratio_syst, color=ROOT.kBlack)

    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadTopMargin(0.065)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.TGaxis.SetMaxDigits(3)

    h_ratio_stat.SetMarkerStyle(ROOT.kFullCircle)
    h_ratio_stat.SetMarkerColor(ROOT.kBlack)
    h_ratio_stat.SetMarkerSize(1.5)
    h_ratio_stat.SetLineWidth(2)
    h_ratio_stat.SetLineColor(ROOT.kBlack)

    g_ratio_syst = ROOT.TGraphErrors(h_ratio_syst)
    g_ratio_syst.SetMarkerStyle(ROOT.kFullCircle)
    g_ratio_syst.SetMarkerColor(ROOT.kBlack)
    g_ratio_syst.SetMarkerSize(1.5)
    g_ratio_syst.SetLineWidth(2)
    g_ratio_syst.SetLineColor(ROOT.kBlack)
    g_ratio_syst.SetFillStyle(0)

    for i in range(0, g_ratio_syst.GetN()):
        g_ratio_syst.SetPointError(i, g_ratio_syst.GetErrorX(i)/2, g_ratio_syst.GetErrorY(i))

    c = ROOT.TCanvas('c', 'c', 700, 600)
    h_frame = c.DrawFrame(1, 5.e4, 23.5, 1.e8, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (pb #kern[-0.5]{#it{c}} / GeV)')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)
    h_ratio_stat.GetYaxis().SetRangeUser( 0.0, 0.25)
    h_ratio_stat.Draw('pE')
    g_ratio_syst.Draw('same 5')

    # optionally add run2 ratios
    if args.input_run2:
        data_run2 = ROOT.TFile.Open(args.input_run2)
        g_run2_stat = data_run2.Get(f'gratio')
        g_run2_syst = data_run2.Get(f'gratio_sys')
        
        r_MB = g_run2_stat.GetY()[0]/br
        r_HM = g_run2_stat.GetY()[1]/br

        r_stat_MB = g_run2_stat.GetErrorYhigh(0)/br
        r_stat_HM = g_run2_stat.GetErrorYhigh(1)/br

        r_syst_MB = (g_run2_syst.GetErrorYhigh(0)+0.09*g_run2_stat.GetY()[0])/br
        r_syst_HM = (g_run2_syst.GetErrorYhigh(1)+0.09*g_run2_stat.GetY()[1])/br

        h_run2_stat = ROOT.TH1F("h_run2_stat", "", 1, 1, 23.5)
        h_run2_stat.SetBinContent(1, r_MB)
        h_run2_stat.SetBinError(1, r_stat_MB)
        h_run2_stat.SetLineColor(ROOT.kRed)
        h_run2_stat.SetLineWidth(0)
        h_run2_stat.SetMarkerStyle(0)
        h_run2_stat.SetMarkerColor(ROOT.kRed)
        h_run2_stat.SetMarkerSize(1.5)
        h_run2_stat.SetFillStyle(3006)
        h_run2_stat.SetFillColorAlpha(ROOT.kRed, 0.6)

        h_run2_stat_HM = ROOT.TH1F("h_run2_stat_HM", "", 1, 1, 23.5)
        h_run2_stat_HM.SetBinContent(1, r_HM)
        h_run2_stat_HM.SetBinError(1, r_stat_HM)
        h_run2_stat_HM.SetLineColor(ROOT.kBlue)
        h_run2_stat_HM.SetLineWidth(0)
        h_run2_stat_HM.SetMarkerStyle(0)
        h_run2_stat_HM.SetMarkerColor(ROOT.kBlue)
        h_run2_stat_HM.SetMarkerSize(1.5)
        h_run2_stat_HM.SetFillStyle(23)
        h_run2_stat_HM.SetFillStyle(3004)
        h_run2_stat_HM.SetFillColorAlpha(ROOT.kBlue, 0.6)

        h_run2_syst = ROOT.TH1F("h_run2_syst", "", 1, 1, 23.5)
        h_run2_syst.SetBinContent(1, r_MB)
        h_run2_syst.SetBinError(1, r_syst_MB)
        h_run2_syst.SetLineColor(ROOT.kRed)
        h_run2_syst.SetLineWidth(0)
        h_run2_syst.SetLineStyle(0)
        h_run2_syst.SetFillStyle(3007)
        h_run2_syst.SetFillColorAlpha(ROOT.kRed, 0.6)

        h_run2_syst_HM = ROOT.TH1F("h_run2_syst_HM", "", 1, 1, 23.5)
        h_run2_syst_HM.SetBinContent(1, r_HM)
        h_run2_syst_HM.SetBinError(1, r_syst_HM)
        h_run2_syst_HM.SetLineColor(ROOT.kBlue)
        h_run2_syst_HM.SetFillStyle(3005)
        h_run2_syst_HM.SetFillColorAlpha(ROOT.kBlue, 0.6)
        h_run2_syst_HM.SetLineWidth(0)
        h_run2_syst_HM.SetLineStyle(0)

        h_run2_stat.Draw('same E2')
        h_run2_stat_HM.Draw('same E2')
        h_run2_syst.Draw('same E2')
        h_run2_syst_HM.Draw('same E2')

        legend = ROOT.TLegend(0.15, 0.7, 0.5, 0.85)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(h_ratio_stat, "Run 3 Stat", "lep")
        legend.AddEntry(g_ratio_syst, "Run 3 Syst", "f")
        legend.AddEntry(h_run2_stat, "Run 2 MB Stat", "f")
        legend.AddEntry(h_run2_syst, "Run 2 MB Syst", "f")
        legend.AddEntry(h_run2_stat_HM, "Run 2 HM Stat", "f")
        legend.AddEntry(h_run2_syst_HM, "Run 2 HM Syst", "f")
        legend.Draw()
        

    c.SaveAs(f"ratio_{args.particle}.png")