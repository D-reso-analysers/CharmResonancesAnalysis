"""
Script for the (non-)prompt fraction calculation of Ds resonances from the one of the D daughter obtained with the cut-variation method
"""

import argparse
import sys
import yaml
import numpy as np
import uproot
from alive_progress import alive_bar

import ROOT  # pylint: disable=import-error
sys.path.insert(0, '..')
from cut_variation import CutVarMinimiser
from utils.style_formatter import set_object_style, set_global_style

# pylint: disable=no-member,too-many-locals,too-many-statements,too-many-branches

def compute_fraction(config):
    """
    Main function
    """

    ROOT.gROOT.SetBatch(True)
    set_global_style(palette=ROOT.kRainBow, padbottommargin=0.14, padleftmargin=0.16)

    with open(config, "r") as f:
        cfg = yaml.safe_load(f)

    pdg_reso = 0
    if cfg["hadron"] == "ds1":
        pdg_reso = 10433
    elif cfg["hadron"] == "ds2star":
        pdg_reso = 435
    else:
        print(f"ERROR: hadron {cfg['hadron']} not supported! Choose ds1 or ds2star")
        sys.exit()

    # compute raw fraction for D daughter
    infile_cutvar = ROOT.TFile.Open(cfg["inputs"]["cutvariation"])
    hist_corry_prompt = infile_cutvar.Get("hCorrYieldsPrompt")
    hist_corry_nonprompt = infile_cutvar.Get("hCorrYieldsNonPrompt")
    hist_covariance = infile_cutvar.Get("hCovPromptNonPrompt")
    hist_corry_prompt.SetDirectory(0)
    hist_corry_nonprompt.SetDirectory(0)
    hist_covariance.SetDirectory(0)
    infile_cutvar.Close()

    infile_eff = ROOT.TFile.Open(cfg["inputs"]["efficiencies"]["file"])
    hist_eff_prompt = infile_eff.Get(cfg["inputs"]["efficiencies"]["histonames"]["prompt"])
    hist_eff_nonprompt = infile_eff.Get(cfg["inputs"]["efficiencies"]["histonames"]["nonprompt"])
    hist_eff_prompt.SetDirectory(0)
    hist_eff_nonprompt.SetDirectory(0)
    infile_eff.Close()

    rel_sys_unc = cfg["rel_sys_unc"]
    graph_frac_prompt_dau_sys = ROOT.TGraphAsymmErrors(1)
    graph_frac_prompt_dau_sys.SetName("graph_frac_prompt_dau_sys")

    hist_frac_prompt_dau = hist_eff_prompt.Clone("hist_frac_prompt_Ddau")
    for ipt in range(1, hist_frac_prompt_dau.GetNbinsX() + 1):
        corry_p = hist_corry_prompt.GetBinContent(ipt)
        corry_np = hist_corry_nonprompt.GetBinContent(ipt)
        unc_corry_p = hist_corry_prompt.GetBinError(ipt)
        unc_corry_np = hist_corry_nonprompt.GetBinError(ipt)
        cov_p_np = hist_covariance.GetBinContent(ipt)
        eff_p = hist_eff_prompt.GetBinContent(ipt)
        eff_np = hist_eff_nonprompt.GetBinContent(ipt)
        minimiser = CutVarMinimiser()
        frac_p, unc_frac_p = minimiser.get_raw_prompt_fraction_ext(corry_p, corry_np, unc_corry_p,
                                                                   unc_corry_np, cov_p_np, eff_p, eff_np)
        hist_frac_prompt_dau.SetBinContent(ipt, frac_p)
        hist_frac_prompt_dau.SetBinError(ipt, unc_frac_p)
        ptcent = hist_frac_prompt_dau.GetBinCenter(ipt)
        ptwidth = hist_frac_prompt_dau.GetBinWidth(ipt)
        graph_frac_prompt_dau_sys.SetPoint(ipt - 1, ptcent, frac_p)
        graph_frac_prompt_dau_sys.SetPointError(ipt - 1, ptwidth / 4, ptwidth / 4,
                                                frac_p * rel_sys_unc[ipt-1], frac_p * rel_sys_unc[ipt-1])

    # convert from pT(D) -> pT(Dreso) using the decay kinematics
    df_kine = uproot.open(cfg["inputs"]["kinematics"])["treeD"].arrays(library="pd")
    df_kine.query(f"pdgDreso == {pdg_reso}", inplace=True)

    pt_maxs = cfg["ptbins"]["maxs"]
    pt_mins = cfg["ptbins"]["mins"]
    pt_max = pt_maxs[-1]
    pt_lims = pt_mins.copy()
    pt_lims.append(pt_max)

    hist_kinematics = ROOT.TH2F("hist_kinematics",
                                ";#it{p}_{T}(D dau) (GeV/#it{c});#it{p}_{T}(D reso) (GeV/#it{c})",
                                int(pt_max * 10), 0., pt_max, int(pt_max * 10), 0., pt_max)
    hist_frac_prompt_reso_vs_pt = {}
    for hist in ["cent", "statmin", "statmax", "sysmin", "sysmax", "scalemin", "scalemax"]:
        hist_frac_prompt_reso_vs_pt[hist] = ROOT.TH2F(f"hist_frac_prompt_reso_vs_pt_{hist}",
                                                      ";#it{p}_{T}(D reso) (GeV/#it{c});#it{f}_{prompt}",
                                                      len(pt_mins), np.array(pt_lims, np.float64),
                                                      1000, 0., 1.)

    scale_factor = cfg["scale_factor"]["value"]
    scale_factor_unc = cfg["scale_factor"]["uncertainty"]

    with alive_bar(len(df_kine), title='Processing kinematics') as bar:
        for pt_dau, pt_reso in zip(df_kine["ptDdau"].to_numpy(), df_kine["ptDreso"].to_numpy()):
            hist_kinematics.Fill(pt_dau, pt_reso)
            pt_dau_bin = hist_frac_prompt_dau.GetXaxis().FindBin(pt_dau)
            if pt_dau_bin < 1:
                pt_dau_bin = 1
            elif pt_dau_bin > hist_frac_prompt_dau.GetNbinsX():
                pt_dau_bin = hist_frac_prompt_dau.GetNbinsX()
            frac_dau = hist_frac_prompt_dau.GetBinContent(pt_dau_bin)
            frac = 1 - ((1 - frac_dau) * scale_factor)
            frac_scalefactmax = 1 - ((1 - frac_dau) * (scale_factor - scale_factor_unc))
            frac_scalefactmin = 1 - ((1 - frac_dau) * (scale_factor + scale_factor_unc))
            statunc_frac_dau = hist_frac_prompt_dau.GetBinError(pt_dau_bin)
            statunc_frac = statunc_frac_dau * scale_factor
            frac_statmin = frac - statunc_frac
            frac_statmax = frac + statunc_frac
            frac_sysmin = frac * (1 - rel_sys_unc[pt_dau_bin - 1])
            frac_sysmax = frac * (1 + rel_sys_unc[pt_dau_bin - 1])

            hist_frac_prompt_reso_vs_pt["cent"].Fill(pt_reso, frac)
            hist_frac_prompt_reso_vs_pt["statmin"].Fill(pt_reso, frac_statmin)
            hist_frac_prompt_reso_vs_pt["statmax"].Fill(pt_reso, frac_statmax)
            hist_frac_prompt_reso_vs_pt["scalemin"].Fill(pt_reso, frac_scalefactmin)
            hist_frac_prompt_reso_vs_pt["scalemax"].Fill(pt_reso, frac_scalefactmax)
            hist_frac_prompt_reso_vs_pt["sysmin"].Fill(pt_reso, frac_sysmin)
            hist_frac_prompt_reso_vs_pt["sysmax"].Fill(pt_reso, frac_sysmax)
            bar()

    hist_frac_prompt_reso = ROOT.TH1F("hist_frac_prompt_reso_vs_pt",
                                      ";#it{p}_{T}(D reso) (GeV/#it{c});#it{f}_{prompt}",
                                      len(pt_mins), np.array(pt_lims, np.float64))
    graph_frac_prompt_reso_sys = ROOT.TGraphAsymmErrors(1)
    graph_frac_prompt_reso_sysscale = ROOT.TGraphAsymmErrors(1)
    graph_frac_prompt_reso_systot = ROOT.TGraphAsymmErrors(1)
    graph_frac_prompt_reso_sys.SetName("graph_frac_prompt_reso_sys")
    graph_frac_prompt_reso_sysscale.SetName("graph_frac_prompt_reso_sysscale")
    graph_frac_prompt_reso_systot.SetName("graph_frac_prompt_reso_systot")
    for ipt in range(1, hist_frac_prompt_reso.GetNbinsX() + 1):
        htmp = hist_frac_prompt_reso_vs_pt["cent"].ProjectionY(f"htmp_{ipt}", ipt, ipt)
        htmp_statmin = hist_frac_prompt_reso_vs_pt["statmin"].ProjectionY(f"htmp_statmin_{ipt}", ipt, ipt)
        htmp_statmax = hist_frac_prompt_reso_vs_pt["statmax"].ProjectionY(f"htmp_statmax_{ipt}", ipt, ipt)
        htmp_sysmin = hist_frac_prompt_reso_vs_pt["sysmin"].ProjectionY(f"htmp_sysmin_{ipt}", ipt, ipt)
        htmp_sysmax = hist_frac_prompt_reso_vs_pt["sysmax"].ProjectionY(f"htmp_sysmax_{ipt}", ipt, ipt)
        htmp_scalemin = hist_frac_prompt_reso_vs_pt["scalemin"].ProjectionY(f"htmp_scalemin_{ipt}", ipt, ipt)
        htmp_scalemax = hist_frac_prompt_reso_vs_pt["scalemax"].ProjectionY(f"htmp_scalemax_{ipt}", ipt, ipt)
        average_frac = htmp.GetMean()
        # we symmetrise all uncertainties for simplicity
        deltas_stat = [average_frac - htmp_statmin.GetMean(), htmp_statmax.GetMean() - average_frac]
        deltas_sys = [average_frac - htmp_sysmin.GetMean(), htmp_sysmax.GetMean() - average_frac]
        deltas_scale = [average_frac - htmp_scalemin.GetMean(), htmp_scalemax.GetMean() - average_frac]
        average_frac_statunc = max(deltas_stat)
        average_frac_sysunc = max(deltas_sys)
        average_frac_scaleunc = max(deltas_scale)
        average_frac_systot = np.sqrt(average_frac_sysunc**2 + average_frac_scaleunc**2)
        hist_frac_prompt_reso.SetBinContent(ipt, average_frac)
        hist_frac_prompt_reso.SetBinError(ipt, average_frac_statunc)
        ptcent = hist_frac_prompt_reso.GetBinCenter(ipt)
        ptwidth = hist_frac_prompt_reso.GetBinWidth(ipt)
        graph_frac_prompt_reso_sys.SetPoint(ipt - 1, ptcent, average_frac)
        graph_frac_prompt_reso_sysscale.SetPoint(ipt - 1, ptcent, average_frac)
        graph_frac_prompt_reso_systot.SetPoint(ipt - 1, ptcent, average_frac)
        graph_frac_prompt_reso_sys.SetPointError(ipt - 1, ptwidth / 4, ptwidth / 4,
                                                 average_frac_sysunc, average_frac_sysunc)
        graph_frac_prompt_reso_sysscale.SetPointError(ipt - 1, ptwidth / 4, ptwidth / 4,
                                                      average_frac_scaleunc, average_frac_scaleunc)
        graph_frac_prompt_reso_systot.SetPointError(ipt - 1, ptwidth / 4, ptwidth / 4,
                                                    average_frac_systot, average_frac_systot)

    set_object_style(hist_frac_prompt_dau, color=ROOT.kAzure+4)
    set_object_style(graph_frac_prompt_dau_sys, color=ROOT.kAzure+4, fillstyle=0)
    set_object_style(hist_frac_prompt_reso, color=ROOT.kRed+1)
    set_object_style(graph_frac_prompt_reso_sys, color=ROOT.kRed+1, fillstyle=0)
    set_object_style(graph_frac_prompt_reso_sysscale, color=ROOT.kRed+1, fillstyle=0)
    set_object_style(graph_frac_prompt_reso_systot, color=ROOT.kRed+1, fillstyle=0)

    canv = ROOT.TCanvas("canv", "", 1500, 500)
    canv.Divide(3, 1)
    canv.cd(1)
    hist_kinematics.DrawCopy("colz")
    pt_min_dau = hist_frac_prompt_dau.GetBinLowEdge(1)
    pt_max_dau = hist_frac_prompt_dau.GetXaxis().GetBinUpEdge(hist_frac_prompt_dau.GetNbinsX())
    canv.cd(2).DrawFrame(pt_min_dau, 0., pt_max_dau, 1., ";#it{p}_{T}(D dau) (GeV/#it{c});#it{f}_{prompt}")
    graph_frac_prompt_dau_sys.Draw("2")
    hist_frac_prompt_dau.DrawCopy("same")
    canv.cd(3).DrawFrame(pt_mins[0], 0., pt_max, 1., ";#it{p}_{T}(D reso) (GeV/#it{c});#it{f}_{prompt}")
    graph_frac_prompt_reso_systot.Draw("2")
    hist_frac_prompt_reso.DrawCopy("same")

    output = ROOT.TFile(cfg["output"], "recreate")
    hist_kinematics.Write()
    for hist in hist_frac_prompt_reso_vs_pt.values():
        hist.Write()
    hist_frac_prompt_dau.Write()
    graph_frac_prompt_dau_sys.Write()
    hist_frac_prompt_reso.Write()
    graph_frac_prompt_reso_sys.Write()
    graph_frac_prompt_reso_sysscale.Write()
    graph_frac_prompt_reso_systot.Write()
    canv.Write()
    output.Close()

    canv.SaveAs(cfg["output"].replace(".root", ".pdf"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument(
        "--config", "-c",
        metavar="text",
        default="config_cutvar_example.yml",
        help="yaml config file",
    )
    args = parser.parse_args()

    compute_fraction(args.config)
