import argparse
import os
import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

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


    h2_eff_mb = infile_mb.Get("h2_eff")
    h2_eff_mb.SetDirectory(0)
    h2_eff_trg = infile_trg.Get("h2_eff")
    h2_eff_trg.SetDirectory(0)
    h2_frac_mb = infile_mb.Get("h2_frac")
    h2_frac_mb.SetDirectory(0)
    h2_frac_trg = infile_trg.Get("h2_frac")
    h2_frac_trg.SetDirectory(0)
    h2_raw_yield_mb = infile_mb.Get("h2_raw_yield")
    h2_raw_yield_mb.SetDirectory(0)
    h2_raw_yield_trg = infile_trg.Get("h2_raw_yield")
    h2_raw_yield_trg.SetDirectory(0)
    h2_xsec_mb = infile_mb.Get("h2_xsec")
    h2_xsec_mb.SetDirectory(0)
    h2_xsec_trg = infile_trg.Get("h2_xsec")
    h2_xsec_trg.SetDirectory(0)
    
    print(f"MB binning X: {h2_raw_yield_mb.GetNbinsX()} Y:{h2_raw_yield_mb.GetNbinsY()}")
    print(f"Trigger binning X: {h2_raw_yield_trg.GetNbinsX()} Y:{h2_raw_yield_trg.GetNbinsY()}")

    pt_centers = []
    for pt_bin in range(1, h2_xsec_mb.GetNbinsX() + 1):
        pt_bin_center = h2_xsec_mb.GetXaxis().GetBinCenter(pt_bin)
        pt_bin_low_edge = h2_xsec_mb.GetXaxis().GetBinLowEdge(pt_bin)
        pt_bin_up_edge = h2_xsec_mb.GetXaxis().GetBinUpEdge(pt_bin)
        print(f"Pt bin {pt_bin}: center={pt_bin_center}, low edge={pt_bin_low_edge}, up edge={pt_bin_up_edge}")
        pt_centers.append(pt_bin_center)
        eff_mb_ipt, eff_trg_ipt, eff_err_mb_ipt, eff_err_trg_ipt = [], [], [], []
        frac_mb_ipt, frac_trg_ipt, frac_err_mb_ipt, frac_err_trg_ipt = [], [], [], []
        raw_yield_mb_ipt, raw_yield_trg_ipt, raw_yield_err_mb_ipt, raw_yield_err_trg_ipt = [], [], [], []
        xsec_mb_ipt, xsec_trg_ipt, xsec_err_mb_ipt, xsec_err_trg_ipt = [], [], [], []
        for conf_bin in range(1, h2_xsec_mb.GetNbinsY() + 1):
            eff_mb_ipt.append( h2_eff_mb.GetBinContent(pt_bin, conf_bin))
            eff_err_mb_ipt.append( h2_eff_mb.GetBinError(pt_bin, conf_bin))
            eff_trg_ipt.append( h2_eff_trg.GetBinContent(pt_bin, conf_bin))
            eff_err_trg_ipt.append( h2_eff_trg.GetBinError(pt_bin, conf_bin))

            frac_mb_ipt.append( h2_frac_mb.GetBinContent(pt_bin, conf_bin))
            frac_err_mb_ipt.append( h2_frac_mb.GetBinError(pt_bin, conf_bin))
            frac_trg_ipt.append( h2_frac_trg.GetBinContent(pt_bin, conf_bin))
            frac_err_trg_ipt.append( h2_frac_trg.GetBinError(pt_bin, conf_bin))

            raw_yield_mb_ipt.append( h2_raw_yield_mb.GetBinContent(pt_bin, conf_bin))
            raw_yield_err_mb_ipt.append( h2_raw_yield_mb.GetBinError(pt_bin, conf_bin))
            raw_yield_trg_ipt.append( h2_raw_yield_trg.GetBinContent(pt_bin, conf_bin))
            raw_yield_err_trg_ipt.append( h2_raw_yield_trg.GetBinError(pt_bin, conf_bin))

            xsec_mb_ipt.append( h2_xsec_mb.GetBinContent(pt_bin, conf_bin))
            xsec_err_mb_ipt.append( h2_xsec_mb.GetBinError(pt_bin, conf_bin))
            xsec_trg_ipt.append( h2_xsec_trg.GetBinContent(pt_bin, conf_bin))
            xsec_err_trg_ipt.append( h2_xsec_trg.GetBinError(pt_bin, conf_bin))
        # print(eff_mb_ipt, frac_mb_ipt, raw_yield_mb_ipt, xsec_mb_ipt)
        fig, axs = plt.subplots(2, 2, figsize=(15, 10))
        pt_low_edge = h2_xsec_mb.GetXaxis().GetBinLowEdge(pt_bin)
        pt_high_edge = h2_xsec_mb.GetXaxis().GetBinUpEdge(pt_bin)
        fig.suptitle(f'{pt_low_edge} < pT < {pt_high_edge}')
        # Efficiency plot
        axs[0, 0].errorbar(range(len(eff_trg_ipt)), eff_trg_ipt, yerr=eff_err_trg_ipt, fmt='o', label=f'TRG (RMS: {np.std(eff_trg_ipt):.4f}, Shift: {eff_trg_ipt[0] - np.mean(eff_trg_ipt):.4f})', color='blue')
        axs[0, 0].errorbar(range(len(eff_trg_ipt), len(eff_trg_ipt) + len(eff_mb_ipt)), eff_mb_ipt, yerr=eff_err_mb_ipt, fmt='o', label=f'MB (RMS: {np.std(eff_mb_ipt):.4f}, Shift: {eff_mb_ipt[0] - np.mean(eff_mb_ipt):.4f})', color='red')
        axs[0, 0].axhline(y=eff_trg_ipt[0], color='blue', linestyle='--')
        axs[0, 0].axhline(y=eff_mb_ipt[0], color='red', linestyle='--')
        axs[0, 0].set_title('Efficiency')
        # axs[0, 0].legend()
        # combined_eff = eff_trg_ipt + eff_mb_ipt
        # axs[0, 0].legend(title=f'Combined RMS: {np.std(combined_eff):.4f}')

        # Fraction plot
        axs[0, 1].errorbar(range(len(frac_trg_ipt)), frac_trg_ipt, yerr=frac_err_trg_ipt, fmt='o', label=f'TRG (RMS: {np.std(frac_trg_ipt):.4f}, Shift: {frac_trg_ipt[0] - np.mean(frac_trg_ipt):.4f})', color='blue')
        axs[0, 1].errorbar(range(len(frac_trg_ipt), len(frac_trg_ipt) + len(frac_mb_ipt)), frac_mb_ipt, yerr=frac_err_mb_ipt, fmt='o', label=f'MB (RMS: {np.std(frac_mb_ipt):.4f}, Shift: {frac_mb_ipt[0] - np.mean(frac_mb_ipt):.4f})', color='red')
        axs[0, 1].axhline(y=frac_trg_ipt[0], color='blue', linestyle='--')
        axs[0, 1].axhline(y=frac_mb_ipt[0], color='red', linestyle='--')
        axs[0, 1].set_title('Fraction')
        # axs[0, 1].legend()
        # combined_frac = frac_trg_ipt + frac_mb_ipt
        # axs[0, 1].legend(title=f'Combined RMS: {np.std(combined_frac):.4f}')

        # Raw Yield plot
        axs[1, 0].errorbar(range(len(raw_yield_trg_ipt)), raw_yield_trg_ipt, yerr=raw_yield_err_trg_ipt, fmt='o', label=f'TRG (RMS: {np.std(raw_yield_trg_ipt):.4f}, Shift: {raw_yield_trg_ipt[0] - np.mean(raw_yield_trg_ipt):.4f})', color='blue')
        axs[1, 0].errorbar(range(len(raw_yield_trg_ipt), len(raw_yield_trg_ipt) + len(raw_yield_mb_ipt)), raw_yield_mb_ipt, yerr=raw_yield_err_mb_ipt, fmt='o', label=f'MB (RMS: {np.std(raw_yield_mb_ipt):.4f}, Shift: {raw_yield_mb_ipt[0] - np.mean(raw_yield_mb_ipt):.4f})', color='red')
        axs[1, 0].axhline(y=raw_yield_trg_ipt[0], color='blue', linestyle='--')
        axs[1, 0].axhline(y=raw_yield_mb_ipt[0], color='red', linestyle='--')
        axs[1, 0].set_title('Raw Yield')
        # axs[1, 0].legend()
        # combined_raw_yield = raw_yield_trg_ipt + raw_yield_mb_ipt
        # axs[1, 0].legend(title=f'Combined RMS: {np.std(combined_raw_yield):.4f}')

        # Cross Section plot
        rms_trg = (np.std(xsec_trg_ipt)/xsec_trg_ipt[0]*100)
        shift_trg = ((xsec_trg_ipt[0] - np.mean(xsec_trg_ipt))/xsec_trg_ipt[0]*100)
        rmsshift_trg = np.sqrt(np.std(xsec_trg_ipt)**2 + (xsec_trg_ipt[0] - np.mean(xsec_trg_ipt))**2)/ xsec_trg_ipt[0]*100
        rms_mb = (np.std(xsec_mb_ipt)/xsec_mb_ipt[0]*100)
        shift_mb = ((xsec_mb_ipt[0] - np.mean(xsec_mb_ipt))/xsec_mb_ipt[0]*100)
        rmsshift_mb = np.sqrt(np.std(xsec_mb_ipt)**2 + (xsec_mb_ipt[0] - np.mean(xsec_mb_ipt))**2)/ xsec_mb_ipt[0]*100
        axs[1, 1].errorbar(range(len(xsec_trg_ipt)), xsec_trg_ipt, yerr=xsec_err_trg_ipt, fmt='o', label=f'TRG (RMS: {rms_trg:.1f}, Shift: {shift_trg:.1f}, RMS + SHIFT: {rmsshift_trg:.1f})', color='blue')
        axs[1, 1].errorbar(range(len(xsec_trg_ipt), len(xsec_trg_ipt) + len(xsec_mb_ipt)), xsec_mb_ipt, yerr=xsec_err_mb_ipt, fmt='o', label=f'MB (RMS: {rms_mb:.1f}, Shift: {shift_mb:.1f}, RMS + SHIFT: {rmsshift_mb:.1f}', color='red')
        axs[1, 1].axhline(y=xsec_trg_ipt[0], color='blue', linestyle='--')
        axs[1, 1].axhline(y=xsec_mb_ipt[0], color='red', linestyle='--')
        axs[1, 1].set_title('Cross Section')
        axs[1, 1].legend()
        xsec_comb_ipt = xsec_trg_ipt + xsec_mb_ipt
        rms_comb = (np.std(xsec_comb_ipt)/xsec_comb_ipt[0]*100)
        shift_comb = ((xsec_comb_ipt[0] - np.mean(xsec_comb_ipt))/xsec_comb_ipt[0]*100)
        rmsshift_comb = np.sqrt(np.std(xsec_comb_ipt)**2 + (xsec_comb_ipt[0] - np.mean(xsec_comb_ipt))**2)/ xsec_comb_ipt[0]*100
        axs[1, 1].legend(title=f'Combined (RMS: {rms_comb:.1f}, Shift: {shift_comb:.1f}, RMS + SHIFT: {rmsshift_comb:.1f})')
        plt.tight_layout()
        fig.savefig(f"output/BDTsys{pt_bin}.png")