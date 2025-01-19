import ROOT
import sys
import argparse

def load_histogram(input_file, histogram_name):
    # Open the ROOT file
    file = ROOT.TFile.Open(input_file)
    if not file or file.IsZombie():
        print("Error: Cannot open file", input_file)
        return

    # Load the histogram (assuming the histogram name is "hist")
    hist = file.Get(histogram_name)
    if not hist:
        print("Error: Cannot find histogram 'hist' in file", input_file)
        return

   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load a histogram from a ROOT file")
    parser.add_argument("input_file", help="Path to the input ROOT file")
    args = parser.parse_args()

    file = ROOT.TFile.Open(args.input_file)
    gFracSyst = file.Get("graph_frac_prompt_reso_systot")
    x_values = []
    y_values = []
    x_errors_low = []
    x_errors_high = []
    y_errors_low = []
    y_errors_high = []
    rel_y_errors_high = []

    for i in range(gFracSyst.GetN()):
        x = gFracSyst.GetPointX(i)
        y = gFracSyst.GetPointY(i)
        x_values.append(float(x))
        y_values.append(float(y))
        x_errors_low.append(float(gFracSyst.GetErrorXlow(i)))
        x_errors_high.append(float(gFracSyst.GetErrorXhigh(i)))
        y_errors_low.append(float(gFracSyst.GetErrorYlow(i)))
        y_errors_high.append(float(gFracSyst.GetErrorYhigh(i)))
        rel_y_errors_high.append(float(gFracSyst.GetErrorYhigh(i))/float(y)*100)

    print("x values:", x_values)
    print("y values:", y_values)
    print("x errors low:", x_errors_low)
    print("x errors high:", x_errors_high)
    print("y errors low:", y_errors_low)
    print("y errors high:", y_errors_high)

    print("relative systematic uncertainty:", rel_y_errors_high)
