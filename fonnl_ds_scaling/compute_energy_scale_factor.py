"""
Script to get energy scale factor for Ds (or in general, a FONLL ratio)
"""

import argparse
import pandas as pd
import ROOT


# pylint: disable=too-many-arguments,too-many-locals
def get_ratio_fonll(df_num, df_den, graph_name, graph_color):
    """
    Function to compute FONLL ratio
    """

    graph_ratio = ROOT.TGraphAsymmErrors(len(df_num))
    graph_ratio.SetNameTitle(graph_name,
                             ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} ratio")
    graph_ratio.SetLineColorAlpha(graph_color, 0.5)
    graph_ratio.SetFillColorAlpha(graph_color, 0.5)
    graph_ratio.SetMarkerColor(graph_color)
    graph_ratio.SetLineWidth(2)
    graph_ratio.SetMarkerStyle(ROOT.kFullCircle)
    # set_object_style(graph_ratio, color=graph_color, alpha=0.5)
    for ipt, (ptmin, ptmax, cent_num, min_sc_num, max_sc_num,
              min_mass_num, max_mass_num, min_pdf_num, max_pdf_num,
              fr_dot5_dot5_num, fr_2_2_num, fr_2_1_num, fr_1_2_num,
              fr_1_dot5_num, fr_dot5_1_num, cent_den, min_sc_den, max_sc_den,
              min_mass_den, max_mass_den, min_pdf_den, max_pdf_den,
              fr_dot5_dot5_den, fr_2_2_den, fr_2_1_den, fr_1_2_den,
              fr_1_dot5_den, fr_dot5_1_den) in enumerate(zip(df_num["ptmin"].to_numpy(),
                                                             df_num["ptmax"].to_numpy(),
                                                             df_num["central"].to_numpy(),
                                                             df_num["min_sc"].to_numpy(),
                                                             df_num["max_sc"].to_numpy(),
                                                             df_num["min_mass"].to_numpy(),
                                                             df_num["max_mass"].to_numpy(),
                                                             df_num["min_pdf"].to_numpy(),
                                                             df_num["max_pdf"].to_numpy(),
                                                             df_num["fr_dot5_dot5"].to_numpy(),
                                                             df_num["fr_2_2"].to_numpy(),
                                                             df_num["fr_2_1"].to_numpy(),
                                                             df_num["fr_1_2"].to_numpy(),
                                                             df_num["fr_1_dot5"].to_numpy(),
                                                             df_num["fr_dot5_1"].to_numpy(),
                                                             df_den["central"].to_numpy(),
                                                             df_den["min_sc"].to_numpy(),
                                                             df_den["max_sc"].to_numpy(),
                                                             df_den["min_mass"].to_numpy(),
                                                             df_den["max_mass"].to_numpy(),
                                                             df_den["min_pdf"].to_numpy(),
                                                             df_den["max_pdf"].to_numpy(),
                                                             df_den["fr_dot5_dot5"].to_numpy(),
                                                             df_den["fr_2_2"].to_numpy(),
                                                             df_den["fr_2_1"].to_numpy(),
                                                             df_den["fr_1_2"].to_numpy(),
                                                             df_den["fr_1_dot5"].to_numpy(),
                                                             df_den["fr_dot5_1"].to_numpy())):
        ptcent = (ptmax + ptmin) / 2
        ptunc = (ptmax - ptmin) / 2
        ratio = cent_num / cent_den
        ratios_vars = [
            min_sc_num / min_sc_den,
            max_sc_num / max_sc_den,
            min_mass_num / min_mass_den,
            max_mass_num / max_mass_den,
            min_pdf_num / min_pdf_den,
            max_pdf_num / max_pdf_den,
            fr_dot5_dot5_num / fr_dot5_dot5_den,
            fr_2_2_num / fr_2_2_den,
            fr_2_1_num / fr_2_1_den,
            fr_1_2_num / fr_1_2_den,
            fr_1_dot5_num / fr_1_dot5_den,
            fr_dot5_1_num / fr_dot5_1_den,
        ]

        graph_ratio.SetPoint(ipt, ptcent, ratio)
        graph_ratio.SetPointError(ipt, ptunc, ptunc,
                                  ratio - min(ratios_vars), max(ratios_vars) - ratio)

    return graph_ratio

def get_corr_factor(infile_fonll_num, infile_fonll_den, outfile_name):
    """
    """

    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetTitleOffset(1.5, "y")
    ROOT.gStyle.SetTitleSize(0.045, "xy")
    ROOT.gStyle.SetLabelSize(0.045, "xy")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gROOT.SetBatch(True)

    col_names = ["ptmin", "ptmax", "central", "min", "max", "min_sc", "max_sc",
                 "min_mass", "max_mass", "min_pdf", "max_pdf", "fr_dot5_dot5", "fr_2_2",
                 "fr_2_1", "fr_1_2", "fr_1_dot5", "fr_dot5_1"]
    df_num = pd.read_csv(infile_fonll_num, names=col_names, comment="#", sep=" ")
    df_den = pd.read_csv(infile_fonll_den, names=col_names, comment="#", sep=" ")

    graph_ratio = get_ratio_fonll(df_num, df_den, "graph_fonll_dmix_13dot6tev_13tev", ROOT.kRed+1)

    canv = ROOT.TCanvas("canv", "", 500, 500)
    frame = canv.DrawFrame(
        min(df_num["ptmin"])-1, 0.9, max(df_num["ptmax"])+1, 1.1,
        ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} FONLL 13.6 TeV / 13 TeV"
    )
    frame.GetYaxis().SetDecimals()
    graph_ratio.Draw("p2")
    canv.SaveAs(outfile_name.replace(".root", ".pdf"))

    outfile = ROOT.TFile(outfile_name, "recreate")
    graph_ratio.Write()
    canv.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infile_num", "-n", metavar="text",
                        default="fonll_dmix_y05_13dot6tev_dsbins.txt",
                        help="fonll file for numerator", required=False)
    parser.add_argument("--infile_den", "-d", metavar="text",
                        default="fonll_dmix_y05_13tev_dsbins.txt",
                        help="fonll file for denominator", required=False)
    parser.add_argument("--outfile", "-o", metavar="text",
                        default="fonll_ratio_dmix_13dot6tev_13tev_dsbins.root",
                        help="output file name", required=False)
    args = parser.parse_args()

    get_corr_factor(args.infile_num, args.infile_den, args.outfile)
