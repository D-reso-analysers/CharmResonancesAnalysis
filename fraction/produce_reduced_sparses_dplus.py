"""
Script to produce simplified THnSparses from taskCharmPolarisation.cxx for D*
"""

import argparse
import numpy as np
import ROOT

def th2_to_sparse(th2, histname):
    """
    Convert a ROOT TH2 histogram into a THnSparse.
    
    Parameters:
        th2 (ROOT.TH2): Input 2D histogram.
    
    Returns:
        ROOT.THnSparseF: A sparse N-dimensional histogram equivalent to the input TH2.
    """
    # Extract binning information from the TH2
    nbins_x = th2.GetNbinsX()
    nbins_y = th2.GetNbinsY()
    
    # Define bin edges for each axis
    xbins = [th2.GetXaxis().GetBinLowEdge(i) for i in range(1, nbins_x + 2)]
    ybins = [th2.GetYaxis().GetBinLowEdge(i) for i in range(1, nbins_y + 2)]
    
    # Create THnSparse axes arrays
    axes_array = [
        (len(xbins) - 1, np.array(xbins, dtype='d')),
        (len(ybins) - 1, np.array(ybins, dtype='d'))
    ]
    
    # Create the THnSparse
    thnsparse = ROOT.THnSparseF(histname, histname, 2, np.array([nbins_x, nbins_y], dtype="i"), np.array([th2.GetXaxis().GetBinLowEdge(nbins_x + 1), th2.GetXaxis().GetBinLowEdge(nbins_y + 1)], dtype="d"))
    thnsparse.SetBinEdges(0, np.array(xbins, dtype='d'))
    thnsparse.SetBinEdges(1, np.array(ybins, dtype='d'))

    # Fill the THnSparse with data from the TH2
    for xbin in range(1, nbins_x + 1):
        for ybin in range(1, nbins_y + 1):
            content = th2.GetBinContent(xbin, ybin)
            if content != 0:
                xcenter = th2.GetXaxis().GetBinCenter(xbin)
                ycenter = th2.GetYaxis().GetBinCenter(ybin)
                coords = np.array([xcenter, ycenter], dtype=np.double)
                thnsparse.Fill(coords, content)

    return thnsparse

def reduce(infile_name, outfile_name, is_mc):
    """

    """

    axes_to_keep = np.array([0, 1, 2, 4], dtype=np.int32)

    infile = ROOT.TFile.Open(infile_name)
    outfile = ROOT.TFile(outfile_name, "recreate")

    if not is_mc:
        sparse = infile.Get(f"hf-task-dplus/hSparseMass")
        outfile.cd()
        sparse_reduced = sparse.Projection(len(axes_to_keep), axes_to_keep)
        sparse_reduced.SetName("hData")
        sparse_reduced.Write()
    else:
        sparse_recop = infile.Get(f"hf-task-dplus/hSparseMassPrompt")
        sparse_reconp = infile.Get(f"hf-task-dplus/hSparseMassFD")
        outfile.cd()
        sparse_recop_reduced = sparse_recop.Projection(len(axes_to_keep), axes_to_keep)
        sparse_reconp_reduced = sparse_reconp.Projection(len(axes_to_keep), axes_to_keep)
        sparse_recop_reduced.SetName("hRecoPrompt")
        sparse_reconp_reduced.SetName("hRecoNonPrompt")
        th2_genp = infile.Get(f"hf-task-dplus/hPtVsYGenPrompt")
        sparse_genp = th2_to_sparse(th2_genp, "hGenPrompt" )
        th2_gennp = infile.Get(f"hf-task-dplus/hPtVsYGenPrompt")
        sparse_gennp = th2_to_sparse(th2_gennp, "hGenNonPrompt" )
        sparse_recop_reduced.Write()
        sparse_reconp_reduced.Write()
        sparse_genp.Write()
        sparse_gennp.Write()

    outfile.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text",
                        default="AnalysisResults_LHC22o_pass7.root",
                        help="input root file", required=False)
    parser.add_argument("--output_file", "-o", metavar="text",
                        default="sparse_dstar_LHC22o_pass7.root",
                        help="output file", required=False)
    parser.add_argument("--is_mc", "-mc", action="store_true",
                        default=False, help="flag for MC")
    args = parser.parse_args()
    reduce(args.input_file, args.output_file, args.is_mc)
