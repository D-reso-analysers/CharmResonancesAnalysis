"""
Compute cross section from output of cut-variation method
"""

import sys
import argparse
import ROOT
sys.path.insert(0, '..')
from utils.style_formatter import set_object_style


def compute_xsec(input_cutvar, input_norm, had, suffix):
    """

    """

    if had == "dplus":
        br = 0.0938
    elif had == "dstar":
        br = 0.677 * 0.03951
    else:
        print(f"ERROR: particle {had} not supported, exit")
        sys.exit()

    infile_cutvar = ROOT.TFile.Open(input_cutvar)
    hist_xsec_p = infile_cutvar.Get("hCorrYieldsPrompt")
    hist_xsec_np = infile_cutvar.Get("hCorrYieldsNonPrompt")
    hist_xsec_p.SetDirectory(0)
    hist_xsec_np.SetDirectory(0)
    hist_xsec_p.SetName("hist_xsec_p")
    hist_xsec_np.SetName("hist_xsec_np")
    hist_xsec_p.Sumw2(0)
    hist_xsec_np.Sumw2(0)
    set_object_style(hist_xsec_p, color=ROOT.kRed+1)
    set_object_style(hist_xsec_np, color=ROOT.kAzure+4)
    infile_cutvar.Close()

    infile_norm = ROOT.TFile.Open(input_norm)
    hist_norm = infile_norm.Get("hist_lumi")
    lumi = hist_norm.GetBinContent(1)
    infile_cutvar.Close()

    hist_xsec_p.Scale(1./2/lumi/br, "width")
    hist_xsec_np.Scale(1./2/lumi/br, "width")

    outfile = ROOT.TFile(f"{had}_xsec_pp13dot6TeV{suffix}.root", "recreate")
    hist_xsec_p.Write()
    hist_xsec_np.Write()
    hist_norm.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_cutvar", "-ic", metavar="text",
                        default="cutvariation/cutvar_dstar_pp13dot6tev_MB_triggerCuts.root",
                        help="input file with cut-variation output")
    parser.add_argument("--input_normalisation", "-in", metavar="text",
                        default="luminosity_LHC22_LHC23.root",
                        help="input file with normalisation")
    parser.add_argument("--particle", "-p", metavar="text",
                        default="dstar",
                        help="Particle species, options [dplus, dstar]")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output file")
    args = parser.parse_args()

    compute_xsec(args.input_cutvar, args.input_normalisation, args.particle, args.suffix)
