import argparse
import os
import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd
# from flarefly.data_handler import DataHandler
# from flarefly.fitter import F2MassFitter
from particle import Particle
from utils.analysis_utils import loadAO2D, applySelections, perform_roofit_fit

def skim_tree(cfg):
    inFileNames = cfg['fileNameData']
    inTreeName = cfg['treeNameData']
    outdir = os.path.dirname(cfg['outputRawYield'])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    dfData = loadAO2D(inFileNames, inTreeName)
    dfFiltered = applySelections(dfData, cfg)
    
    dfFiltered.to_parquet(os.path.join(outdir, cfg['filteredDf']))
    print(dfFiltered)


def compute_significance(workspace):
    # significance computation
    data = workspace.data("data")

    meanR = workspace.var("meanR")
    sigmaR = workspace.var("sigmaR")
    widthR = workspace.var("widthR")
    meanR.setConstant(1)
    sigmaR.setConstant(1)
    widthR.setConstant(1)
    nSigR = workspace.var("nSigR")

    totPdf = workspace.pdf("totPdf")

    model = ROOT.RooStats.ModelConfig()
    model.SetWorkspace(workspace)
    model.SetPdf("totPdf")
    poi = ROOT.RooArgSet(nSigR)
    nullParams = poi.snapshot()
    nullParams.setRealValue("nSigR",0.)

    plc = ROOT.RooStats.ProfileLikelihoodCalculator()

    plc.SetData(data)
    plc.SetModel(model)
    plc.SetParameters(poi)
    plc.SetNullParameters(nullParams)

    #We get a HypoTestResult out of the calculator, and we can query it.
    hypo_test_result = plc.GetHypoTest()
    significance = hypo_test_result.Significance()
    print ("-------------------------------------------------")
    print ("The p-value for the null is ", hypo_test_result.NullPValue())
    print ("Corresponding to a signifcance of ", significance)
    print ("-------------------------------------------------")
    #input()
    del plc
    return significance


def fit_with_roofit(cfg):
    pt_mins = cfg["cutVars"]["pt"]["min"]
    pt_maxs = cfg["cutVars"]["pt"]["max"]
    outdir = os.path.dirname(cfg['outputRawYield'])
    outFileName = cfg['outputRawYield']
    outfile = ROOT.TFile(outFileName, "recreate")
    df = pd.read_parquet(os.path.join(outdir, cfg['filteredDf']))
    pdgId = cfg['pdgId']
    
    # Define output histograms
    pt_limits = pt_mins.copy()
    pt_limits.append(pt_maxs[-1])
    pt_limits = np.array(pt_limits, np.float64)

    hist_rawyield = ROOT.TH1F(
        "hist_rawyield", ";#it{p}_{T} (GeV/#it{c}); raw yield", len(pt_mins), pt_limits)
    hist_sigma = ROOT.TH1F(
        "hist_sigma", ";#it{p}_{T} (GeV/#it{c}); #sigma (GeV/#it{c}^{2})", len(pt_mins), pt_limits)
    hist_mean = ROOT.TH1F(
        "hist_mean", ";#it{p}_{T} (GeV/#it{c}); #mu (GeV/#it{c}^{2})", len(pt_mins), pt_limits)
    
    hist_rawyield.SetDirectory(0)
    hist_sigma.SetDirectory(0)
    hist_mean.SetDirectory(0)

    if pdgId == 10433:
            resoName = "Ds1"
    if pdgId == 435:
        resoName = "Ds2Star"
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        df_pt = df[(df['fPt'] >= pt_min) & (df['fPt'] < pt_max)]
        df_pt = df_pt[(df_pt['fM'] >= cfg["fit"]["mass_mins"][ipt]) & (df_pt['fM'] < cfg["fit"]["mass_maxs"][ipt])]
        # Perform fit
        fit_params = cfg["fit"]["parameters"]
        workspace = perform_roofit_fit(df_pt, "voigtian", "thr3",  fit_params, 0.5, 0.6)
        # Retrieve fit results
        mass = workspace.var("mass")
        data = workspace.data("data")
        meanR = workspace.var("mean")
        sigmaR = workspace.var("sigma")
        nSigR = workspace.var("nSig")
        nBkg = workspace.var("nBkg")
        totPdf = workspace.pdf("totPdf")

        hist_rawyield.SetBinContent(ipt+1, nSigR.getVal())
        hist_rawyield.SetBinError(ipt+1, nSigR.getError())
        hist_mean.SetBinContent(ipt+1, meanR.getVal())
        hist_mean.SetBinError(ipt+1, meanR.getError())
        hist_sigma.SetBinContent(ipt+1, sigmaR.getVal())
        hist_sigma.SetBinError(ipt+1, sigmaR.getError())

        # Create a frame to draw the fit result and data
        mass_frame = mass.frame( cfg["fit"]["mass_mins"][ipt], cfg["fit"]["mass_maxs"][ipt], cfg["plot"]["nbins"][ipt])
        mass_frame.SetTitle(f"{resoName} mass pt {pt_min}-{pt_max}")
        data.plotOn(mass_frame)
        totPdf.plotOn(mass_frame, ROOT.RooFit.Components("sigRPdf"), ROOT.RooFit.LineColor(ROOT.kAzure + 2), ROOT.RooFit.FillColor(ROOT.kAzure-2))
        totPdf.plotOn(mass_frame, ROOT.RooFit.Components("bkgPdf"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
        totPdf.plotOn(mass_frame, ROOT.RooFit.LineColor(ROOT.kBlue))
    
        hresid = mass_frame.residHist()

        chi2 = mass_frame.chiSquare()  # Compute chi-squared/ndf
        legend = ROOT.TLegend(0.3, 0.15, 0.8, 0.45)  # Create a legend at a specified position
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)  # Transparent background
        legend.SetTextSize(0.035)
        canvas_yield = ROOT.TCanvas("canvas_yield", "canvas_yield", 900, 900)
        legend.AddEntry("", f"#chi^{{2}}/ndf = {chi2:.2f}", "")  # Add chi-squared
        legend.AddEntry("", f"Mean = {(meanR.getVal()*1000):.2f} #pm {(meanR.getError()*1000):.2f} MeV/c^{{2}}", "")  # Add mean
        legend.AddEntry("", f"Sigma = {(sigmaR.getVal()*1000):.2f} #pm {(sigmaR.getError()*1000):.2f} MeV/c^{{2}}", "")  # Add sigma
        legend.AddEntry("", f"Nsig = {nSigR.getVal():.0f} #pm {nSigR.getError():.0f}", "")
        # signifcance = compute_significance(workspace)

        # legend.AddEntry("", f"Significance = {signifcance:.1f}", "")
        print(f"#chi^{{2}}/ndf = {chi2:.2f}")
        print(f"Mean = {meanR.getVal():.4f} #pm {meanR.getError():.4f} GeV/c^2")  # Add mean
        print(f"Sigma = {sigmaR.getVal():.4f} #pm {sigmaR.getError():.4f} GeV/c^2")  # Add sigma
        print(f"Nsig = {nSigR.getVal():.0f} #pm {nSigR.getError():.0f}")
        # input()
        # Draw the frame
        
        mass_frame_res = mass.frame( cfg["fit"]["mass_mins"][ipt], cfg["fit"]["mass_maxs"][ipt], cfg["plot"]["nbins"][ipt])
        data.plotOn(mass_frame_res)
        hist = mass_frame.getHist("data")  # 0 if it's the first histogram plotted

        # Check if the histogram exists and get the maximum value
        if hist:
            max_value = hist.GetMaximum()
            y_max = 1.5 * max_value  # Scale the maximum value
            mass_frame.SetMaximum(y_max)  # Set the y-axis maximum
        else:
            print("No histogram found in the RooPlot. Ensure data is plotted first.")
        totPdf.plotOn(mass_frame_res, ROOT.RooFit.Components("sigRPdf"), ROOT.RooFit.LineColor(ROOT.kGreen))
        totPdf.plotOn(mass_frame_res, ROOT.RooFit.LineColor(ROOT.kBlue))
        totPdf.plotOn(mass_frame_res, ROOT.RooFit.Components("bkgPdf"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
    
        hresid = mass_frame_res.residHist()


        # mf2 =  mass.frame( cfg["fit"]["mass_mins"][ipt], cfg["fit"]["mass_maxs"][ipt], cfg["plot"]["nbins"][ipt])
        # mass_frame_res.SetTitle(f"Residual {resoName} mass pt {pt_min}-{pt_max}")
        # data.plotOn(mf2, ROOT.RooFit.ResLevel(RooAbsData.Residual))
        # totPdf.plotOn(mass_frame_res, ROOT.RooFit.Components("sigRPdf"), ROOT.RooFit.LineColor(ROOT.kGreen))

        xaxis = mass_frame.GetXaxis()
        x_min = xaxis.GetXmin()
        x_max = xaxis.GetXmax()
        n_bins = xaxis.GetNbins()

        # Calculate the bin width
        bin_width = (x_max - x_min) / n_bins

        # Format the y-axis title with bin width
        y_title = f"Counts per {(bin_width*1000):.0f} MeV/c^{{2}}"
        mass_frame.GetYaxis().SetTitle(y_title)
        mass_frame.GetXaxis().SetTitle(f"M(K#pi#piK^{{0}}_{{S}})-M(K#pi#pi) (GeV/c^{{2}})")

        # Optionally adjust text size and positioning
        mass_frame.GetYaxis().SetTitleSize(0.035)
        mass_frame.GetYaxis().SetTitleOffset(1.5)

        mass_frame.Draw()
        legend.Draw("same")
        # canvas_res = ROOT.TCanvas("canvas_res", "canvas_res", 800, 600)
        # mf2.Draw()
        
       
        canvas_yield.SaveAs(os.path.join(outdir, "plots", f"{cfg['plotLabel']}_pt{pt_min}-{pt_max}.pdf"))
        # canvas_res.SaveAs(os.path.join(outdir, f"{cfg['plotLabel']}_res_pt{pt_min}-{pt_max}.pdf"))
    outfile.cd()
    hist_rawyield.Write()
    hist_sigma.Write()
    hist_mean.Write()
    outfile.Close()
    
# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('configfile', 
                        metavar='text',
                        help='config file name with inputs and cut configurations')
    parser.add_argument("--skim", "-s", action="store_true",
                        default=False, help="apply cuts and save df to parquet")
    parser.add_argument("--fit", "-f", action="store_true",
                        default=False, help="enable fit w/o cut from parquet")
    args = parser.parse_args()
    with open(args.configfile, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    if args.skim:
        skim_tree(cfg)
    if args.fit:
        fit_with_roofit(cfg)
    
