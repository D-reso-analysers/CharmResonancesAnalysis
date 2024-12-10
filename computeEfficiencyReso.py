import argparse
import os
import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd
from utils.analysis_utils import loadAO2D, applySelections

def computeEfficiency(dfRec, dfGen, cfg):
    def eff(nRec, nGen):
        return nRec / nGen
    def effErr(nRec, nGen):
        return nRec/nGen * np.sqrt((1-nRec/nGen)/nRec)
    cutVars = cfg['cutVars']
    ptLabel = cutVars['pt']['varname']
    originLabel = 'fOrigin'
    effPrompt , effErrPrompt, effNonPrompt, effErrNonPrompt = [], [], [], []
    for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['pt']['min'], cutVars['pt']['max'])):
        dfRecCut = dfRec[(dfRec[f"{ptLabel}"] >= ptMin) & (dfRec[f"{ptLabel}"] < ptMax)]
        dfGenCut = dfGen[(dfGen[f"{ptLabel}"] >= ptMin) & (dfGen[f"{ptLabel}"] < ptMax)]
        nRecPrompt = len(dfRecCut[dfRecCut[f"{originLabel}"] == 1])
        nRecNonPrompt = len(dfRecCut[dfRecCut[f"{originLabel}"] == 2])
        nGenPrompt = len(dfGenCut[dfGenCut[f"{originLabel}"] == 1])
        nGenNonPrompt = len(dfGenCut[dfGenCut[f"{originLabel}"] == 2])
        effPrompt.append(eff(nRecPrompt, nGenPrompt))
        effErrPrompt.append(effErr(nRecPrompt, nGenPrompt))
        effNonPrompt.append(eff(nRecNonPrompt, nGenNonPrompt))
        effErrNonPrompt.append(effErr(nRecNonPrompt, nGenNonPrompt))
    return effPrompt , effErrPrompt, effNonPrompt, effErrNonPrompt

def fillHisto(eff, effUnc, cfg):
    lower_edges = cfg['cutVars']['pt']['min']
    upper_edges = cfg['cutVars']['pt']['max']
    bin_edges = np.array(lower_edges + [upper_edges[-1]], dtype=np.float64)
    hist = ROOT.TH1F("hEff", "Efficiency histogram", len(bin_edges) - 1, bin_edges)
    hist.GetXaxis().SetTitle("pT [GeV/c]")
    hist.GetYaxis().SetTitle("acc x efficiency")

    for bin in range(1, hist.GetNbinsX() + 1):
        iBin = bin -1
        hist.SetBinContent(bin, eff[iBin])
        hist.SetBinError(bin, effUnc[iBin])
    return hist
    
# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('configfile', 
                        metavar='text',
                        help='config file name with inputs and cut configurations')
    args = parser.parse_args()
    with open(args.configfile, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    inFileNames = cfg['fileNameMC']
    inTreeNameRec = cfg['treeNameRec']
    inTreeNameGen = cfg['treeNameGen']
    # Load trees and apply cuts on Reconstructed candidates
    dfRec = loadAO2D(inFileNames, inTreeNameRec)
    dfRecFiltered = applySelections(dfRec, cfg, isMC=True)
    dfGen = loadAO2D (inFileNames, inTreeNameGen)
    # Compute efficiencies
    effPrompt , effErrPrompt, effNonPrompt, effErrNonPrompt = computeEfficiency(dfRecFiltered, dfGen, cfg)
    # Make histograms 
    hEffPrompt = fillHisto(effPrompt, effErrPrompt, cfg)
    hEffNonPrompt = fillHisto(effNonPrompt, effErrNonPrompt, cfg)
    hEffPrompt.SetName('hEffPrompt')
    hEffPrompt.SetTitle('hEffPrompt')
    hEffNonPrompt.SetName('hEffNonPrompt')
    hEffNonPrompt.SetTitle('hEffNonPrompt')
    # Drawing options and canvas
    hEffPrompt.SetMarkerColor(2)
    hEffPrompt.SetLineColor(2)
    hEffPrompt.SetMarkerStyle(20)
    hEffNonPrompt.SetMarkerStyle(20)
    hEffNonPrompt.SetMarkerColor(4)
    hEffNonPrompt.SetLineColor(4)
    legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)
    legend.SetHeader("hEff", "C")  # Optional header for the legend
    legend.AddEntry(hEffPrompt, "prompt", "pl") 
    legend.AddEntry(hEffNonPrompt, f"non-prompt", "pl") 
    c = ROOT.TCanvas('cEff', 'cEff', 900, 400)
    c.SetLogy()
    c.SetGridy()
    # c.SetTitle("Ds2Star Efficiency vs pT")
    hEffPrompt.Draw()
    hEffNonPrompt.Draw('same')
    legend.Draw('same')
    # Save histograms and canvas in a rootFile
    outFileName = cfg['outputEff']
    outdir = os.path.dirname(outFileName)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = ROOT.TFile(outFileName, 'recreate')
    hEffPrompt.Write()
    hEffNonPrompt.Write()
    c.Write()