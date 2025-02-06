import ROOT
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_frac.py file1.root file2.root")
        sys.exit(1)

    file1 = ROOT.TFile.Open(sys.argv[1])
    file2 = ROOT.TFile.Open(sys.argv[2])

    hist1 = file1.Get("hist_frac_prompt_reso_vs_pt")
    hist2 = file2.Get("hRawFracPrompt")
    hist3 = file1. Get("hist_frac_prompt_Ddau")
    for i in range(1, hist1.GetNbinsX() + 1):
        content = hist1.GetBinContent(i)
        hist1.SetBinError(i, 0.05 * content)
    c = ROOT.TCanvas("c", "compare_frac", 800, 600)
    hist1.GetYaxis().SetRangeUser(0.7, 1)
    hist1.GetXaxis().SetRangeUser(5, 24)
    hist2.GetXaxis().SetRangeUser(5, 24)

    hist1.SetLineColor(ROOT.kRed)
    hist2.SetLineColor(ROOT.kBlue)
    hist3.SetLineColor(ROOT.kBlack)

    hist1.SetMarkerColor(ROOT.kRed)
    hist2.SetMarkerColor(ROOT.kBlue)
    hist3.SetMarkerColor(ROOT.kBlack)

    legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
    legend.AddEntry(hist2, "direct cut variation", "lpe")
    legend.AddEntry(hist1, "fraction propagation", "lpe")
    legend.AddEntry(hist3, "D* fraction", "lpe")


    hist1.Draw('PE')
    hist2.Draw("PE SAME")
    hist3.Draw("PE SAME")
    legend.Draw()


    c.SaveAs("compare_frac.pdf")

if __name__ == "__main__":
    main()
    