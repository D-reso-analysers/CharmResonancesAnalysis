"""
Script to produce simplified THnSparses from taskCharmPolarisation.cxx for D*
"""

import argparse
import numpy as np
import ROOT

def reduce(infile_name, outfile_name, is_mc):
    """

    """

    axes_to_keep = np.array([0, 1, 6, 7], dtype=np.int32)

    infile = ROOT.TFile.Open(infile_name)
    outfile = ROOT.TFile(outfile_name, "recreate")
    found_sparse = False

    for axis in ["Helicity", "Production", "Random", "Beam"]:
        if found_sparse: # one is enough
            break

        if not is_mc:
            sparse = infile.Get(f"task-polarisation-charm-hadrons/h{axis}")
            if not isinstance(sparse, ROOT.THnSparse):
                continue

            found_sparse = True
            y_bin_min = sparse.GetAxis(3).FindBin(-0.7999)
            y_bin_max = sparse.GetAxis(3).FindBin(0.7999)
            sparse.GetAxis(3).SetRange(y_bin_min, y_bin_max)

            outfile.cd()
            sparse_reduced = sparse.Projection(len(axes_to_keep), axes_to_keep)
            sparse_reduced.SetName("hData")
            sparse_reduced.Write()
        else:
            sparse_recop = infile.Get("task-polarisation-charm-hadrons/hRecoPrompt")
            if not isinstance(sparse_recop, ROOT.THnSparse):
                continue

            found_sparse = True
            sparse_reconp = infile.Get("task-polarisation-charm-hadrons/hRecoNonPrompt")
            y_bin_min = sparse_recop.GetAxis(3).FindBin(-0.7999)
            y_bin_max = sparse_recop.GetAxis(3).FindBin(0.7999)
            sparse_recop.GetAxis(3).SetRange(y_bin_min, y_bin_max)
            sparse_reconp.GetAxis(3).SetRange(y_bin_min, y_bin_max)

            sparse_genp = infile.Get("task-polarisation-charm-hadrons/hGenPrompt")
            sparse_gennp = infile.Get("task-polarisation-charm-hadrons/hGenNonPrompt")

            outfile.cd()
            sparse_recop_reduced = sparse_recop.Projection(len(axes_to_keep), axes_to_keep)
            sparse_recop_reduced.SetName("hRecoPrompt")
            sparse_recop_reduced.Write()
            sparse_reconp_reduced = sparse_reconp.Projection(len(axes_to_keep), axes_to_keep)
            sparse_reconp_reduced.SetName("hRecoNonPrompt")
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
