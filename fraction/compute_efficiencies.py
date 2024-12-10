"""
Script for the computation of the efficiencies for several BDT output scores for cut-variation method
"""

import os
import sys
import argparse
import numpy as np
import yaml
import ROOT



def calculate_efficiency(num, den):
    """

    """

    eff = num / den
    unc = np.sqrt(eff * (1 - eff) / den)

    return eff, unc

def compute_efficiencies(input_config):
    """

    """

    with open(input_config, "r") as f:
        cfg = yaml.safe_load(f)

    if cfg["hadron"] not in ["dstar", "dplus"]:
        print(f"ERROR: {cfg['hadron']} not supported, exit")
        sys.exit()

    outputdir = cfg["output"]["efficiencies"]["directory"]
    suffix = cfg["output"]["efficiencies"]["suffix"]

    pt_mins = cfg["pt_mins"]
    pt_maxs = cfg["pt_maxs"]
    pt_limits = pt_mins.copy()
    pt_limits.append(pt_maxs[-1])
    pt_limits = np.array(pt_limits, np.float64)

    bdt_np_mins = cfg["bdt_cuts"]["nonprompt"]

    infile_name = os.path.join(outputdir, f"hist_mcpt{suffix}.root")

    hist_eff_p_nocut = ROOT.TH1F(
        "hist_eff_prompt", ";#it{p}_{T} (GeV/#it{c}); #font[152]{e} #times Acc",
        len(pt_mins), pt_limits)
    hist_eff_np_nocut = ROOT.TH1F(
        "hist_eff_nonprompt", ";#it{p}_{T} (GeV/#it{c}); #font[152]{e} #times Acc",
        len(pt_mins), pt_limits)
    hist_eff_p_nocut.SetDirectory(0)
    hist_eff_np_nocut.SetDirectory(0)

    hist_eff_p, hist_eff_np = [], []
    for icut, _ in enumerate(bdt_np_mins):
        hist_eff_p.append(
            ROOT.TH1F("hist_eff_prompt", ";#it{p}_{T} (GeV/#it{c}); #font[152]{e} #times Acc",
                      len(pt_mins), pt_limits))
        hist_eff_np.append(
            ROOT.TH1F("hist_eff_nonprompt", ";#it{p}_{T} (GeV/#it{c}); #font[152]{e} #times Acc",
                      len(pt_mins), pt_limits))
        hist_eff_p[icut].SetDirectory(0)
        hist_eff_np[icut].SetDirectory(0)

    infile = ROOT.TFile.Open(infile_name)
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        hist_gen_p = infile.Get(f"hist_genpt_p_pt{pt_min:.1f}_{pt_max:.1f}")
        hist_gen_np = infile.Get(f"hist_genpt_np_pt{pt_min:.1f}_{pt_max:.1f}")
        n_gen_p = hist_gen_p.Integral()
        n_gen_np = hist_gen_np.Integral()
        hist_rec_p = infile.Get(
            f"hist_recopt_p_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp")
        hist_rec_np = infile.Get(
            f"hist_recopt_np_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp")
        n_rec_p = hist_rec_p.Integral()
        n_rec_np = hist_rec_np.Integral()
        eff_p, unc_p = calculate_efficiency(n_rec_p, n_gen_p)
        eff_np, unc_np = calculate_efficiency(n_rec_np, n_gen_np)
        hist_eff_p_nocut.SetBinContent(ipt+1, eff_p)
        hist_eff_p_nocut.SetBinError(ipt+1, unc_p)
        hist_eff_np_nocut.SetBinContent(ipt+1, eff_np)
        hist_eff_np_nocut.SetBinError(ipt+1, unc_np)
        for icut, bdt_np_min in enumerate(bdt_np_mins):
            hist_rec_p = infile.Get(
                f"hist_recopt_p_pt{pt_min:.1f}_{pt_max:.1f}_bdtnp{bdt_np_min:.2f}")
            hist_rec_np = infile.Get(
                f"hist_recopt_np_pt{pt_min:.1f}_{pt_max:.1f}_bdtnp{bdt_np_min:.2f}")
            n_rec_p = hist_rec_p.Integral()
            n_rec_np = hist_rec_np.Integral()
            eff_p, unc_p = calculate_efficiency(n_rec_p, n_gen_p)
            eff_np, unc_np = calculate_efficiency(n_rec_np, n_gen_np)
            hist_eff_p[icut].SetBinContent(ipt+1, eff_p)
            hist_eff_p[icut].SetBinError(ipt+1, unc_p)
            hist_eff_np[icut].SetBinContent(ipt+1, eff_np)
            hist_eff_np[icut].SetBinError(ipt+1, unc_np)

    outfile_name_cutvar = os.path.join(
        outputdir, f"efficiencies_nocutnp{suffix}.root")
    outfile = ROOT.TFile(outfile_name_cutvar, "recreate")
    hist_eff_p_nocut.Write()
    hist_eff_np_nocut.Write()
    outfile.Close()

    for icut, bdt_np_min in enumerate(bdt_np_mins):
        outfile_name_cutvar = os.path.join(
            outputdir, f"efficiencies_bdtnp{bdt_np_min:0.2f}{suffix}.root")
        outfile = ROOT.TFile(outfile_name_cutvar, "recreate")
        hist_eff_p[icut].Write()
        hist_eff_np[icut].Write()
        outfile.Close()


# function to project the sparse
def project(input_config):
    """

    """

    with open(input_config, "r") as f:
        cfg = yaml.safe_load(f)

    if cfg["hadron"] not in ["dstar", "dplus"]:
        print(f"ERROR: {cfg['hadron']} not supported, exit")
        sys.exit()

    infile_names = cfg["input"]["mc"]
    mc_weights = cfg["input"]["mc_weights"]
    outputdir = cfg["output"]["efficiencies"]["directory"]
    suffix = cfg["output"]["efficiencies"]["suffix"]
    include_pitomu = cfg["efficiency"]["include_pitomu"]

    sparse_recop, sparse_reconp, sparse_genp, sparse_gennp = ([] for _ in range(4))
    for infile_name in infile_names:
        infile = ROOT.TFile.Open(infile_name)
        sparse_recop.append(infile.Get("hRecoPrompt"))
        sparse_reconp.append(infile.Get("hRecoNonPrompt"))
        sparse_genp.append(infile.Get("hGenPrompt"))
        sparse_gennp.append(infile.Get("hGenNonPrompt"))
    infile.Close()

    pt_mins = cfg["pt_mins"]
    pt_maxs = cfg["pt_maxs"]
    bdt_bkg_cuts = cfg["bdt_cuts"]["bkg"]
    if not isinstance(bdt_bkg_cuts, list):
        bdt_bkg_cuts = [bdt_bkg_cuts]*len(pt_mins)

    outfile = ROOT.TFile(os.path.join(outputdir, f"hist_mcpt{suffix}.root"), "recreate")

    histos_recpt_p, histos_recpt_np = [], []
    histos_recpt_p_nocut, histos_recpt_np_nocut = [], []
    histos_genpt_p, histos_genpt_np = [], []
    axis_y_gen = 2 if cfg["hadron"] == "dstar" else 1
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        # generated
        pt_bin_min = sparse_genp[0].GetAxis(0).FindBin(pt_min*1.001)
        pt_bin_max = sparse_genp[0].GetAxis(0).FindBin(pt_max*0.999)
        # |y|<0.5 by default
        y_bin_min = sparse_genp[0].GetAxis(axis_y_gen).FindBin(-0.4999)
        y_bin_max = sparse_genp[0].GetAxis(axis_y_gen).FindBin(0.4999)
        for ifile, _ in enumerate(infile_names):
            sparse_genp[ifile].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
            sparse_gennp[ifile].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
            sparse_genp[ifile].GetAxis(axis_y_gen).SetRange(y_bin_min, y_bin_max)
            sparse_gennp[ifile].GetAxis(axis_y_gen).SetRange(y_bin_min, y_bin_max)
            if not include_pitomu: # exlcude pi->mu decays from signal
                sparse_recop[ifile].GetAxis(4).SetRange(1, 1)
                sparse_reconp[ifile].GetAxis(4).SetRange(1, 1)
            if ifile == 0:
                histos_genpt_p.append(sparse_genp[ifile].Projection(0))
                histos_genpt_np.append(sparse_gennp[ifile].Projection(0))
                histos_genpt_p[ipt].Sumw2()
                histos_genpt_np[ipt].Sumw2()
                if len(infile_names) > 1:
                    histos_genpt_p[ipt].Scale(mc_weights[ifile])
                    histos_genpt_np[ipt].Scale(mc_weights[ifile])
            else:
                histos_genpt_p[ipt].Add(sparse_genp[ifile].Projection(0), mc_weights[ifile])
                histos_genpt_np[ipt].Add(sparse_gennp[ifile].Projection(0), mc_weights[ifile])
        histos_genpt_p[ipt].SetName(f"hist_genpt_p_pt{pt_min:.1f}_{pt_max:.1f}")
        histos_genpt_np[ipt].SetName(f"hist_genpt_np_pt{pt_min:.1f}_{pt_max:.1f}")
        outfile.cd()
        histos_genpt_p[ipt].Write()
        histos_genpt_np[ipt].Write()

        # reco (no np cut)
        bdt_bkg_bin_max = sparse_recop[0].GetAxis(2).FindBin(bdt_bkg_cuts[ipt]*0.999)
        pt_bin_min = sparse_recop[0].GetAxis(1).FindBin(pt_min*1.001)
        pt_bin_max = sparse_recop[0].GetAxis(1).FindBin(pt_max*0.999)
        for ifile, _ in enumerate(infile_names):
            sparse_recop[ifile].GetAxis(2).SetRange(1, bdt_bkg_bin_max)
            sparse_reconp[ifile].GetAxis(2).SetRange(1, bdt_bkg_bin_max)
            sparse_recop[ifile].GetAxis(1).SetRange(pt_bin_min, pt_bin_max)
            sparse_reconp[ifile].GetAxis(1).SetRange(pt_bin_min, pt_bin_max)
            if ifile == 0:
                histos_recpt_p_nocut.append(sparse_recop[ifile].Projection(1))
                histos_recpt_np_nocut.append(sparse_reconp[ifile].Projection(1))
                histos_recpt_p_nocut[ipt].Sumw2()
                histos_recpt_np_nocut[ipt].Sumw2()
                if len(infile_names) > 1:
                    histos_recpt_p_nocut[ipt].Scale(mc_weights[ifile])
                    histos_recpt_np_nocut[ipt].Scale(mc_weights[ifile])
            else:
                histos_recpt_p_nocut[ipt].Add(sparse_recop[ifile].Projection(1), mc_weights[ifile])
                histos_recpt_np_nocut[ipt].Add(sparse_reconp[ifile].Projection(1), mc_weights[ifile])
        histos_recpt_p_nocut[ipt].SetName(
            f"hist_recopt_p_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp")
        histos_recpt_np_nocut[ipt].SetName(
            f"hist_recopt_np_pt{pt_min:.1f}_{pt_max:.1f}_nocutnp")
        outfile.cd()
        histos_recpt_p_nocut[ipt].Write()
        histos_recpt_np_nocut[ipt].Write()

        # reco (np cuts)
        histos_recpt_p.append([])
        histos_recpt_np.append([])
        for icut, bdt_np_min in enumerate(cfg["bdt_cuts"]["nonprompt"]):
            bdt_np_bin_min = sparse_recop[0].GetAxis(3).FindBin(bdt_np_min*1.001)
            for ifile, _ in enumerate(infile_names):
                sparse_recop[ifile].GetAxis(3).SetRange(bdt_np_bin_min, sparse_recop[0].GetAxis(3).GetNbins()+1)
                sparse_reconp[ifile].GetAxis(3).SetRange(bdt_np_bin_min, sparse_reconp[0].GetAxis(3).GetNbins()+1)
                if ifile == 0:
                    histos_recpt_p[ipt].append(sparse_recop[ifile].Projection(1))
                    histos_recpt_np[ipt].append(sparse_reconp[ifile].Projection(1))
                    histos_recpt_p[ipt][icut].Sumw2()
                    histos_recpt_np[ipt][icut].Sumw2()
                    if len(infile_names) > 1:
                        histos_recpt_p[ipt][icut].Scale(mc_weights[ifile])
                        histos_recpt_np[ipt][icut].Scale(mc_weights[ifile])
                else:
                    histos_recpt_p[ipt][icut].Add(sparse_recop[ifile].Projection(1), mc_weights[ifile])
                    histos_recpt_np[ipt][icut].Add(sparse_reconp[ifile].Projection(1), mc_weights[ifile])
            histos_recpt_p[ipt][icut].SetName(
                f"hist_recopt_p_pt{pt_min:.1f}_{pt_max:.1f}_bdtnp{bdt_np_min:.2f}")
            histos_recpt_np[ipt][icut].SetName(
                f"hist_recopt_np_pt{pt_min:.1f}_{pt_max:.1f}_bdtnp{bdt_np_min:.2f}")
            outfile.cd()
            histos_recpt_p[ipt][icut].Write()
            histos_recpt_np[ipt][icut].Write()
            for ifile, _ in enumerate(infile_names):
                sparse_recop[ifile].GetAxis(3).SetRange(-1, -1)
                sparse_reconp[ifile].GetAxis(3).SetRange(-1, -1)
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--cfg_file", "-c", metavar="text",
                        default="config.yml", help="config file")
    parser.add_argument("--project", "-p", action="store_true",
                        default=False, help="enable projection")
    parser.add_argument("--doeff", "-e", action="store_true",
                        default=False, help="enable efficiency computation")
    args = parser.parse_args()

    if args.project:
        project(args.cfg_file)

    if args.doeff:
        compute_efficiencies(args.cfg_file)
