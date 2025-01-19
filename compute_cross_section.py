import argparse
import ROOT
import sys
import yaml
from utils.style_formatter import set_object_style

br_dplus = 0.0938
br_dstar = 0.677 * 0.03951
br_k0 = 0.692
br_ds1 = 0.22
br_ds2star = 0.2335

def have_same_binning(hist1, hist2):
    """
    Check if two ROOT histograms have the same binning edges.

    Args:
        hist1 (ROOT.TH1): First histogram.
        hist2 (ROOT.TH1): Second histogram.

    Returns:
        bool: True if the histograms have the same binning edges, False otherwise.
    """
    if hist1.GetNbinsX() != hist2.GetNbinsX():
        return False

    for i in range(hist1.GetNbinsX() + 1):  # Include underflow/overflow edges
        if hist1.GetBinLowEdge(i) != hist2.GetBinLowEdge(i):
            return False
    
    return True

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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('configfile', 
                        metavar='text',
                        help='config file name with inputs and cut configurations')
    args = parser.parse_args()
    with open(args.configfile, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    if cfg['particle'] == "ds1":
        br = br_ds1 * br_dstar * br_k0
    elif cfg['particle'] == "ds2star":
        br = br_ds2star * br_dplus * br_k0
    else:
        print(f"ERROR: particle {cfg['particle']} not supported, exit")
        sys.exit()

    infile_yield = ROOT.TFile.Open(cfg['rawYieldFile'])
    infile_eff = ROOT.TFile.Open(cfg['effFile'])
    infile_frac = ROOT.TFile.Open(cfg['fracFile'])
    infile_norm = ROOT.TFile.Open(cfg['luminosityFile'])

    hist_yield = infile_yield.Get(cfg["rawYieldHist"])
    hist_yield.SetDirectory(0)
    hist_eff = infile_eff.Get(cfg["effHist"])
    hist_eff.SetDirectory(0)
    hist_frac = infile_frac.Get(cfg["fracHist"])
    hist_frac.SetDirectory(0)
    hist_norm = infile_norm.Get(cfg["luminosityHist"])
    hist_norm.SetDirectory(0)
    lumi = hist_norm.GetBinContent(1)

    syst = cfg['systematics']['relSyst']

    if not have_same_binning(hist_yield, hist_eff):
        print(f"ERROR: efficiency and raw yield histograms have different binnings!")
        sys.exit()
    if not have_same_binning(hist_yield, hist_frac):
        print(f"ERROR: fraction and raw yield histograms have different binnings!")
        sys.exit()

    hist_xsec_stat = hist_yield.Clone()
    hist_xsec_stat.SetDirectory(0)
    hist_xsec_stat.Reset()
    hist_xsec_stat.SetName("hist_xsec_stat")
    hist_xsec_stat.GetYaxis().SetTitle("d#sigma/d#it{p}_{T} (#mub GeV^{-1} #it{c})")
    hist_xsec_stat.Sumw2(0)
    set_object_style(hist_xsec_stat, color=ROOT.kBlack)

    hist_xsec_syst = hist_xsec_stat.Clone()
    hist_xsec_syst.SetName("hist_xsec_syst")

    for ipt in range(hist_yield.GetNbinsX()):
        dpt = hist_yield.GetBinLowEdge(ipt + 2) - hist_yield.GetBinLowEdge(ipt + 1)
        raw_yield = hist_yield.GetBinContent(ipt + 1)
        eff = hist_eff.GetBinContent(ipt + 1)
        prompt_frac = hist_frac.GetBinContent(ipt + 1)
        xsec = compute_cross_section(raw_yield, eff, prompt_frac, dpt, br, lumi)
        xsec_err = xsec * hist_yield.GetBinError(ipt + 1)/hist_yield.GetBinContent(ipt + 1)
        hist_xsec_stat.SetBinContent(ipt + 1, xsec)
        hist_xsec_stat.SetBinError(ipt + 1, xsec_err)
        hist_xsec_syst.SetBinContent(ipt + 1, xsec)
        syst_err = xsec*syst[ipt]
        hist_xsec_syst.SetBinError(ipt + 1, syst_err)

    outfile = ROOT.TFile(cfg['outfile'], 'RECREATE')
    hist_xsec_stat.Write()
    hist_xsec_syst.Write()
    outfile.Close()

    # -------- Make pretty plot-------------
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadTopMargin(0.065)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.TGaxis.SetMaxDigits(3)
    
    hist_xsec_stat.SetMarkerStyle(ROOT.kFullCircle)
    hist_xsec_stat.SetMarkerColor(ROOT.kBlack)
    # hist_xsec_stat.SetMarkerSize(1.5)
    hist_xsec_stat.SetLineWidth(2)
    hist_xsec_stat.SetLineColor(ROOT.kBlack)

    g_xsec_syst = ROOT.TGraphErrors(hist_xsec_syst)

    g_xsec_syst.SetMarkerStyle(ROOT.kFullCircle)
    g_xsec_syst.SetMarkerColor(ROOT.kBlack)
    # g_xsec_syst.SetMarkerSize(1.5)
    g_xsec_syst.SetLineWidth(2)
    g_xsec_syst.SetLineColor(ROOT.kBlack)
    g_xsec_syst.SetFillStyle(0)
    for i in range(0, g_xsec_syst.GetN()):
        g_xsec_syst.SetPointError(i, g_xsec_syst.GetErrorX(i)/2, g_xsec_syst.GetErrorY(i))

    c = ROOT.TCanvas('c', 'c', 700, 600)
    c.SetLogy()
    h_frame = c.DrawFrame(2, 0.06, 24.0, 0.3, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (pb #kern[-0.5]{#it{c}} / GeV)')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)
    hist_xsec_stat.Draw('same pE')
    g_xsec_syst.Draw('same 5')

    # Add the text
    if cfg['particle'] == 'ds1':
        text_decay = ROOT.TLatex(0.55, 0.85, 'D_{s1}^{+} #rightarrow D*^{+}K_{S}^{0}')
    text_decay.SetNDC()
    text_decay.SetTextSize(0.04)
    text_decay.SetTextFont(42)
    text_decay.Draw()

    text_conj = ROOT.TLatex(0.55, 0.8, 'and charge conjugate')
    text_conj.SetNDC()
    text_conj.SetTextSize(0.04)
    text_conj.SetTextFont(42)
    text_conj.Draw()

    text_ALICE = ROOT.TLatex(0.55, 0.7, 'Work in Progress')
    text_ALICE.SetNDC()
    text_ALICE.SetTextSize(0.06)
    text_ALICE.SetTextFont(42)
    text_ALICE.Draw()

    text_pp = ROOT.TLatex(0.55, 0.65, 'pp collisions, #sqrt{#it{s}} = 13.6 TeV')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.04)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    text_rapidity = ROOT.TLatex(0.55, 0.6, '|y| < 0.5')
    text_rapidity.SetNDC()
    text_rapidity.SetTextSize(0.04)
    text_rapidity.SetTextFont(42)
    text_rapidity.Draw()
    
    ROOT.gPad.RedrawAxis()

    c.SaveAs(cfg['plotFile'])