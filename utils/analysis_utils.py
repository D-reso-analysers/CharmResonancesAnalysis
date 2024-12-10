import yaml
import ROOT
import uproot
import numpy as np
import pandas as pd

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
                print(f'cutting on {var} > {varMin} for pT Bachelor {ptMinBach0}-{ptMaxBach0}')
                dfCut = dfCut[dfCut[f"{varLabel}"] >= varMin]
                alreadyCut.append(var)
            if varMax is not None:
                print(f'cutting on {var} < {varMax} for pT Bachelor {ptMinBach0}-{ptMaxBach0}')
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
                print(f'cutting on {var} > {varMin} for pT Candidate {ptMin}-{ptMax}')
                dfCut = dfCut[dfCut[f"{varLabel}"] >= varMin]
            if varMax is not None:
                print(f'cutting on {var} < {varMax} for pT Candidate {ptMin}-{ptMax}')
                dfCut = dfCut[dfCut[f"{varLabel}"] < varMax]
        dfOut = pd.concat([dfOut, dfCut])
    print(f"Kept {len(dfOut)} candidates")
    return dfOut