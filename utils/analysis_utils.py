import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd
import sys

def loadAO2D(inFileNames, inTreeName):
    '''
    utility to load tree from AO2D and return its content in a dataframe
    '''
    dataframe = pd.DataFrame()
    if not isinstance(inFileNames, list):
            inFileNames = [inFileNames]
    for file in inFileNames:
        print(f'Loading trees form {file}')
        f = uproot.open(file)
        dirnames = f.keys(recursive = False)
        print(f"keys found in file: {dirnames}")
        for dirname in dirnames:
            if 'DF' in dirname:
                path = f'{file}:{dirname}/{inTreeName}'
                dataframe = pd.concat([dataframe, uproot.open(path).arrays(library='pd')], ignore_index=True)
                print(f"processed key: {dirname}")
            else:
                continue
    return dataframe

def applyParametrizedMassCut(df, cfg):
    '''
    utility to apply parametrized invariant mass cut (for D+ bachelor)
    '''
    cutVars = cfg['cutVars']
    ptLabel = cutVars['ptBach0']['varname']
    params = cutVars['invMassProng0']['params']
    massLabel = cutVars['invMassProng0']['varname']

    pdgMass = 1.8697
    # Compute peakMean
    peakMean = np.where(
        df[ptLabel] < params['ptMaxDelta'],
        (pdgMass + params['deltaMassConst'] +
        params['deltaMassLinear'] * df[ptLabel]),
        pdgMass
    )

    # Compute peakWidth
    peakWidth = (params['sigmaConst'] +
                params['sigmaLinear'] * df[ptLabel])

    # Final condition
    df['pass_cut'] = ~(np.abs(df[massLabel] - peakMean) >
                        params['nSigma'] * peakWidth)

    # Filter rows that pass the cut
    dfOut = df[df['pass_cut']]

    return dfOut

def applySelections(dfIn, cfg, isMC = False):
    '''
    utility to apply pt differential cuts on dataframes
    '''
    # TO DO: add BDT selections as a function of pt of D-meson not resonance
    cutVars = cfg['cutVars']
    ptLabel = cutVars['pt']['varname']
    selectedFlags = cfg['acceptFlags']
    dfOut = pd.DataFrame()
    dfFiltered = pd.DataFrame()
    # Selection on candidate type:
    if isMC:
        for flag in selectedFlags:
            dfPart = dfIn[np.abs(dfIn["fFlagMcMatchRec"]) == flag]
            dfFiltered = pd.concat([dfFiltered, dfPart])
    else:
        dfFiltered = dfIn
    # cuts differential in bachelor pT
    print(f"Starting Number of candidates: {len(dfFiltered)}")
    dfBachCut = pd.DataFrame()
    alreadyCut = []
    if cutVars['invMassProng0']['parametrized']:
        print("Applying parametrized mass cut for D+")
        dfProv = pd.DataFrame()
        dfProv = applyParametrizedMassCut(dfFiltered, cfg)
        alreadyCut.append('invMassProng0')
        dfFiltered = dfProv
    for iPtBach0, (ptMinBach0, ptMaxBach0) in enumerate(zip(cutVars['ptBach0']['min'], cutVars['ptBach0']['max'])):
        dfCut = dfFiltered[(dfFiltered[f"{cutVars['ptBach0']['varname']}"] >= ptMinBach0) & (dfFiltered[f"{cutVars['ptBach0']['varname']}"] < ptMaxBach0)]
        for var in cutVars:
            if not cutVars[var]['useBachPt'] or var == 'ptBach0':
                continue
            if var == 'invMassProng0':
                if cutVars[var]['parametrized']:
                    continue
            varMin = cutVars[var]['min'][iPtBach0]
            varMax = cutVars[var]['max'][iPtBach0]
            varLabel = cutVars[var]['varname']
            if varMin is not None:
                # print(f'cutting on {var} > {varMin} for pT Bachelor {ptMinBach0}-{ptMaxBach0}')
                dfCut = dfCut[dfCut[f"{varLabel}"] >= varMin]
                alreadyCut.append(var)
            if varMax is not None:
                # print(f'cutting on {var} < {varMax} for pT Bachelor {ptMinBach0}-{ptMaxBach0}')
                dfCut = dfCut[dfCut[f"{varLabel}"] < varMax]
                alreadyCut.append(var)
        dfBachCut = pd.concat([dfBachCut, dfCut])
    # cuts differential in candidate pT
    for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['pt']['min'], cutVars['pt']['max'])):
        dfCut = dfBachCut[(dfBachCut[f"{ptLabel}"] >= ptMin) & (dfBachCut[f"{ptLabel}"] < ptMax)]
        for var in cutVars:
            if var == 'pt' or var in alreadyCut or var == 'ptBach0':
                continue
            varMin = cutVars[var]['min'][iPt]
            varMax = cutVars[var]['max'][iPt]
            varLabel = cutVars[var]['varname']
            if varMin is not None:
                # print(f'cutting on {var} > {varMin} for pT Candidate {ptMin}-{ptMax}')
                dfCut = dfCut[dfCut[f"{varLabel}"] >= varMin]
            if varMax is not None:
                # print(f'cutting on {var} < {varMax} for pT Candidate {ptMin}-{ptMax}')
                dfCut = dfCut[dfCut[f"{varLabel}"] < varMax]
        dfOut = pd.concat([dfOut, dfCut])
    print(f"Kept {len(dfOut)} candidates")
    return dfOut

def applySelections2(dfIn, cutVars, varsToSkip = [], isMC = False, selectedFlags = None):
    '''
    utility to apply pt differential cuts on dataframes
    '''
    # TO DO: add BDT selections as a function of pt of D-meson not resonance
    ptLabel = cutVars['pt']['varname']
    dfOut = pd.DataFrame()
    dfFiltered = pd.DataFrame()
    if isMC and not selectedFlags:
        print("for MC a list of MC flags to keep should be provided")
        sys.exit()
    # Selection on candidate type:
    if isMC:
        for flag in selectedFlags:
            dfPart = dfIn[np.abs(dfIn["fFlagMcMatchRec"]) == flag]
            dfFiltered = pd.concat([dfFiltered, dfPart])
    else:
        dfFiltered = dfIn
    # cuts differential in bachelor pT
    print(f"Starting Number of candidates: {len(dfFiltered)}")
    dfBachCut = pd.DataFrame()
    alreadyCut = []
    if cutVars['invMassProng0']['parametrized']:
        print("Applying parametrized mass cut for D+")
        dfProv = pd.DataFrame()
        dfProv = applyParametrizedMassCut(dfFiltered, cfg)
        alreadyCut.append('invMassProng0')
        dfFiltered = dfProv
    for iPtBach0, (ptMinBach0, ptMaxBach0) in enumerate(zip(cutVars['ptBach0']['min'], cutVars['ptBach0']['max'])):
        dfCut = dfFiltered[(dfFiltered[f"{cutVars['ptBach0']['varname']}"] >= ptMinBach0) & (dfFiltered[f"{cutVars['ptBach0']['varname']}"] < ptMaxBach0)]
        for var in cutVars:
            if var in varsToSkip:
                print(f"skipping {var}")
                continue
            if not cutVars[var]['useBachPt'] or var == 'ptBach0':
                continue
            if var == 'invMassProng0':
                if cutVars[var]['parametrized']:
                    continue
            varMin = cutVars[var]['min'][iPtBach0]
            varMax = cutVars[var]['max'][iPtBach0]
            varLabel = cutVars[var]['varname']
            if varMin is not None:
                # print(f'cutting on {var} > {varMin} for pT Bachelor {ptMinBach0}-{ptMaxBach0}')
                dfCut = dfCut[dfCut[f"{varLabel}"] >= varMin]
                alreadyCut.append(var)
            if varMax is not None:
                # print(f'cutting on {var} < {varMax} for pT Bachelor {ptMinBach0}-{ptMaxBach0}')
                dfCut = dfCut[dfCut[f"{varLabel}"] < varMax]
                alreadyCut.append(var)
        dfBachCut = pd.concat([dfBachCut, dfCut])
    # cuts differential in candidate pT
    for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['pt']['min'], cutVars['pt']['max'])):
        dfCut = dfBachCut[(dfBachCut[f"{ptLabel}"] >= ptMin) & (dfBachCut[f"{ptLabel}"] < ptMax)]
        for var in cutVars:
            if var in varsToSkip:
                print(f"skipping {var}")
                continue
            if var == 'pt' or var in alreadyCut or var == 'ptBach0':
                continue
            varMin = cutVars[var]['min'][iPt]
            varMax = cutVars[var]['max'][iPt]
            varLabel = cutVars[var]['varname']
            if varMin is not None:
                # print(f'cutting on {var} > {varMin} for pT Candidate {ptMin}-{ptMax}')
                dfCut = dfCut[dfCut[f"{varLabel}"] >= varMin]
            if varMax is not None:
                # print(f'cutting on {var} < {varMax} for pT Candidate {ptMin}-{ptMax}')
                dfCut = dfCut[dfCut[f"{varLabel}"] < varMax]
        dfOut = pd.concat([dfOut, dfCut])
    print(f"Kept {len(dfOut)} candidates")
    return dfOut

def set_param(param_list, name, value, min_value, max_value):
    """
    Replaces the value, min, and max of a parameter in the paramList.
    If a parameter with the given name does not exist is added to the list instead

    :param param_list: List of parameter dictionaries.
    :param name: The name of the parameter to update.
    :param value: The new value for the parameter.
    :param min_value: The new minimum value for the parameter.
    :param max_value: The new maximum value for the parameter.
    :return: None. The param_list is updated in place.
    """
    for param in param_list:
        if param["name"] == name:
            param["value"] = value
            param["min"] = min_value
            param["max"] = max_value
            break
    else:
        par_dict = {
            "name": name,
            "value": value,
            "min": min_value,
            "max": max_value 
        }
        param_list.append(par_dict)

def initPars(parNameList, parListInput):
    '''
    Function to initialize the parameters for the fit

    Args:
    parNameList (list): List of parameter names to be initialized.
    parListInput (list): List of dictionaries with fit parameters to initialize of the form:
        {   
            "name": "mean",
            "value": "1.87",
            "min": "1.86",
            "max": "1.88" 
        }

    Returns:
    parDict: dictionary of RooRealVars with the parameters initialized
    '''
    parDict = {}
    alreadySet = False
    for pName in parNameList:
        for p in parListInput:
            if p['name'] == pName:
                print(f"Setting parameter {p['name']} with value {p['value']}")
                if p['min'] == p['max']:
                    parDict[p['name']] = ROOT.RooRealVar(p['name'], p['name'], float(p['value']))
                    parDict[p['name']].setConstant(ROOT.kTRUE)
                else:
                    parDict[p['name']] = ROOT.RooRealVar(p['name'], p['name'], float(p['value']), float(p['min']), float(p['max']))
                alreadySet = True
        if not alreadySet:
            print(f"Setting parameter {pName} with default values")
            parDict[pName] = ROOT.RooRealVar(pName, pName, 0, -1000, 1000)
        alreadySet = False
    return parDict

def perform_roofit_fit(df, signalPdfName, bkgPdfName, paramList, massMin, massMax, massName = 'fM'):
    '''
    Function to perform a RooFit fit on a dataset

    Args:
    df: pandas DataFrame with the data to fit
    signalPdfName: string with the name of the signal PDF
    bkgPdfName: string with the name of the background PDF
    paramList: list of dictionaries with fit parmaeters of the form:
        {   "name": "mean",
            "value": "1.87",
            "min": "1.86",
            "max": "1.88" 
        }
    massMin: float with the minimum mass value
    massMax: float with the maximum mass value
    massName: string with the name of the mass variable in the DataFrame

    Returns:
    workspace: ROOT.RooWorkspace with the fit results and data
    '''
    # Create the mass variable and dataset
    mass = ROOT.RooRealVar("mass", "Invariant Mass [GeV/c^2]", massMin, massMax)
    data = ROOT.RooDataSet("data", "Dataset from Pandas", ROOT.RooArgSet(mass))
    for m in df[massName]:
        if m < massMin or m > massMax:
            continue
        mass.setVal(m)
        data.add(ROOT.RooArgSet(mass))

    # signal PDFs
    if signalPdfName not in ["voigtian", "gaussian"]:
        print(f"Signal PDF {signalPdfName} not implemented")
        return
    if signalPdfName == "voigtian":
        parNamesSignal = ["mean", "sigma", "width"]
        parsDictSignal = initPars(parNamesSignal, paramList)
        sigPdf = ROOT.RooVoigtian("sigPdf", "Voigtian Pdf", mass, parsDictSignal['mean'], parsDictSignal['width'], parsDictSignal['sigma'])
    elif signalPdfName == "gaussian":
        parNamesSignal = ["mean", "sigma"]
        parsDictSignal = initPars(parNamesSignal, paramList)
        sigPdf = ROOT.RooGaussian("sigPdf", "Gaussian Pdf", mass, parsDictSignal['mean'], parsDictSignal['sigma'])

    # background PDFs
    if bkgPdfName not in ["exp", "thrn", "thr1", "thr2", "thr3", "cheb2", "cheb3", "pol1"]:
        print(f"Background PDF {bkgPdfName} not implemented")
        return
    if bkgPdfName == "thrn":
        parNamesBkg = ["mTh", "l", "alpha", "n"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh)^n)"
        bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, parsDictBkg['mTh'], parsDictBkg['alpha'], parsDictBkg['l'], parsDictBkg['n']))
    elif bkgPdfName == "thr1":
        parNamesBkg = ["mTh", "l", "alpha"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh))"
        bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, parsDictBkg['mTh'], parsDictBkg['alpha'], parsDictBkg['l']))
    elif bkgPdfName == "thr2":
        parNamesBkg = ["mTh", "l", "alpha", "beta"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh) + beta * (mass - mTh) * (mass - mTh))"
        bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, parsDictBkg['mTh'], parsDictBkg['alpha'], parsDictBkg['beta'], parsDictBkg['l']))
    elif bkgPdfName == "thr3":
        parNamesBkg = ["mTh", "l", "alpha", "beta", "gamma"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh) + beta * (mass - mTh) * (mass - mTh)+ gamma * (mass - mTh)* (mass - mTh)* (mass - mTh))"
        bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, parsDictBkg['mTh'], parsDictBkg['alpha'], parsDictBkg['beta'], parsDictBkg['gamma'], parsDictBkg['l']))
    elif bkgPdfName == "cheb2":
        parNamesBkg = ["a0", "a1", "a2"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        bkgPdf = ROOT.RooChebychev("bkgPdf", "Chebychev Pdf", mass, ROOT.RooArgList(parsDictBkg['a0'], parsDictBkg['a1'], parsDictBkg['a2']))
    elif bkgPdfName == "cheb3":   
        parNamesBkg = ["a0", "a1", "a2", "a3"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        bkgPdf = ROOT.RooChebychev("bkgPdf", "Chebychev Pdf", mass, ROOT.RooArgList(parsDictBkg['a0'], parsDictBkg['a1'], parsDictBkg['a2'], parsDictBkg['a3']))
    elif bkgPdfName == "pol1":
        parNamesBkg = ["a0", "a1"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        bkgPdf = ROOT.RooPolynomial("bkgPdf", "Polynomial Pdf", mass, ROOT.RooArgList(parsDictBkg['a0'], parsDictBkg['a1']))
    elif bkgPdfName == "exp":
        parNamesBkg = ["a0"]
        parsDictBkg = initPars(parNamesBkg, paramList)
        bkgPdf = ROOT.RooExponential("bkgPdf", "Exponential Pdf", mass, parsDictBkg['a0'])

    nSig = ROOT.RooRealVar("nSig","nSig",0.1*len(df), 0 ,len(df))
    nBkg = ROOT.RooRealVar("nBkg","nBkg",0.5*len(df), 0 ,len(df))
    
    totPdf = ROOT.RooAddPdf("totPdf", "Total PDF with Nsig as parameter", ROOT.RooArgList(sigPdf, bkgPdf), ROOT.RooArgList(nSig, nBkg))

    fit_result = totPdf.fitTo(data, ROOT.RooFit.Save())
    
    workspace = ROOT.RooWorkspace(f"workspace", f"workspace")
    getattr(workspace, 'import')(totPdf)
    getattr(workspace, 'import')(data)
    getattr(workspace, 'import')(fit_result, "fitResults")

    return workspace

def get_chi2_significance_sb(ws, nbins, outfile = None, name = None):
    '''
    Function to get the chi2 and significance of the fit

    Parameters:
    ws: RooWorkspace
        Workspace containing the fit results
    nbins: int
        Number of bins in the histogram
    outfile: ROOT.TFile
        Output file to save the canvas
    name: str

    Returns:
    chi2: float
        Chi2 of the fit
    significance: float
        Significance of the fit
    '''
    mass = ws.var("mass")
    data = ws.data("data")
    model = ws.pdf("totPdf")
    bkgPdf = ws.pdf("bkgPdf")

    frame = mass.frame(ROOT.RooFit.Bins(nbins))
    data.plotOn(frame)
    model.plotOn(frame, ROOT.RooFit.Components("bkgPdf"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
    model.plotOn(frame)

    chi2 = frame.chiSquare()

    s = ws.var("nSig").getVal()
    sigma = ws.var("sigma").getVal()
    mean = ws.var("mean").getVal()
    b = ws.var("nBkg").getVal()
    
    mass.setRange("signal_range", mean - 3 * sigma, mean + 3 * sigma)
    integral_bkg = bkgPdf.createIntegral(ROOT.RooArgSet(mass), ROOT.RooFit.NormSet(ROOT.RooArgSet(mass)), ROOT.RooFit.Range("signal_range"))
    integral_value = integral_bkg.getVal()

    significance = s / np.sqrt(s + b * integral_value)
    sb = s / (b * integral_value)
    print(f" Signal: {s}, Background: {b*integral_value}, Chi2: {chi2}, Significance: {significance}")

    if outfile and name:
        c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 800, 800)
        
        legend = ROOT.TLegend(0.4, 0.55, 0.9, 0.85)  # Create a legend at a specified position
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)  # Transparent background
        legend.SetTextSize(0.035)
        legend.AddEntry("", f"#chi^{{2}}/ndf = {chi2:.2f}", "")  # Add chi-squared
        legend.AddEntry("", f"Mean = {(mean*1000):.2f} #pm {(ws.var("mean").getError()*1000):.2f} MeV/c^{{2}}", "")  # Add mean
        legend.AddEntry("", f"#sigma = {(sigma*1000):.2f} #pm {(ws.var("sigma").getError()*1000):.2f} MeV/c^{{2}}", "")  # Add sigma
        legend.AddEntry("", f"S = {s:.0f} #pm {ws.var("nSig").getError():.0f}", "")
        legend.AddEntry("", f"B(3 #sigma) = {(b * integral_value):.0f} #pm {(ws.var("nBkg").getError() * np.sqrt(integral_value)):.0f}", "")
        legend.AddEntry("", f"Significance = {significance:.1f}", "")  # Add significance
        legend.AddEntry("", f"S/B = {s/(b * integral_value):.3f}", "")  # Add s/b

        frame.SetTitle(f"{name}")
        frame.Draw()
        legend.Draw()
        outfile.cd()
        c.Write()

    return chi2, significance, sb