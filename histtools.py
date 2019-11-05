import ROOT

MYCOLZ = [ROOT.kBlack,
       ROOT.TColor.GetColor("#00B5E2"), # light blue
       ROOT.TColor.GetColor("#F68D2E"), # light orange
       ROOT.TColor.GetColor("#78BE21"), # light green
       ROOT.TColor.GetColor("#AF272F"), # light red
       ROOT.TColor.GetColor("#EAAA00"), # dark gold
       ROOT.TColor.GetColor("#004C97"), # blue
       ROOT.TColor.GetColor("#4C8C2C"), # dark green
]
MYCOLZ=[ROOT.kBlack,ROOT.kBlue,ROOT.kRed,ROOT.kViolet]

def DrawHistsEmSim(em_hists, sim_hists, plotname,
                   ytitle="object multiplicity",order=["Tracks","EM calo","Calo","Muons"]):
    #configure the canvas, pad
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("canv","",750,750)
    pad = ROOT.TPad("pad", "pad", .005, .01, .995, .995)
    pad.Draw()
    pad.cd()

    ROOT.gPad.SetLeftMargin(0.17)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetTopMargin(0.05)

    #draw this guy first
    b = em_hists[0].Clone("tmp")
    b.Reset()
    b.Draw()
    mymax = max(
        max([h.GetMaximum() for h in em_hists]),
        max([h.GetMaximum() for h in sim_hists])
    )
    b.SetMaximum(1.4*mymax)
    b.Draw()
    
    b.GetXaxis().SetTitleOffset(1.2)
    b.GetYaxis().SetTitleOffset(1.5)
    b.GetYaxis().SetLabelSize(0.04)
    b.GetYaxis().SetTitleSize(0.04)
    b.GetYaxis().SetTitle(ytitle)
    b.Draw("axis")

    #col=[ROOT.kBlack,ROOT.kBlue,ROOT.kRed,ROOT.kViolet]
    for i,h in enumerate(em_hists):
        h.SetLineWidth(2)
        h.SetLineColor(MYCOLZ[i])
        h.Draw("hist same")
    for i,h in enumerate(sim_hists):
        h.SetMarkerStyle(20)
        h.SetMarkerColor(MYCOLZ[i])
        h.Draw("p0 same")
    
    ll = ROOT.TLatex()
    ll.SetNDC()
    ll.SetTextFont(62)
    ll.DrawLatex(.22,.89,"#scale[0.9]{CMS}")
    ll.SetTextFont(42)
    ll.DrawLatex(.32,.89,"#scale[0.9]{Preliminary}")
    #ll.DrawLatex(.30,.89,"#scale[0.9]{Internal}")

    ll.SetTextFont(42)
    ll.DrawLatex(.22,.84,"#scale[0.7]{t#bar{t} events, #LT#mu#GT=200}")

    #create and fill legend
    leg = ROOT.TLegend(.22,.73, 0.87,.88)
    leg.SetTextFont(42)
    leg.SetHeader("")
    leg.SetNColumns(3)
    #leg.AddEntry(sim_hists[0], "Simulation","pel")
    leg.AddEntry(sim_hists[0], "Firmware","p")
    leg.AddEntry(em_hists[0],  order[0],"l")
    leg.AddEntry(em_hists[1],  order[1],"l")
    #leg.AddEntry(em_hists[0],  "Emulation","l")
    leg.AddEntry(em_hists[0],  "Emulation","f")
    leg.AddEntry(em_hists[2],  order[2],"l")
    leg.AddEntry(em_hists[3],  order[3],"l")

    #draw the legend
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.Draw()
    
    c.SaveAs(plotname+".pdf")

def DrawHistsLinkClk(hists, plotname,
                   ytitle="# objects",order=["Tracks","EM calo","Calo","Muons"]):
    #configure the canvas, pad
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("canv","",750,350)
    pad = ROOT.TPad("pad", "pad", .005, .01, .995, .995)
    pad.Draw()
    pad.cd()

    ROOT.gPad.SetLeftMargin(0.07)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetTopMargin(0.05)

    #draw this guy first
    b = hists[0].Clone("tmp")
    b.Reset()
    b.Draw()
    mymax = max([h.GetMaximum() for h in hists])
    b.SetMaximum(1.4*mymax)
    b.Draw()
    
    b.GetXaxis().SetTitleOffset(1.2)
    b.GetYaxis().SetTitleOffset(0.7)
    b.GetYaxis().SetLabelSize(0.04)
    b.GetYaxis().SetTitleSize(0.04)
    b.GetYaxis().SetTitle(ytitle)
    b.Draw("axis")

    #col=[ROOT.kBlack,ROOT.kBlue,ROOT.kRed,ROOT.kViolet]


    for i,h in enumerate(hists):
        h.SetLineWidth(2)
        h.SetLineColor(MYCOLZ[i])
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.6)
        h.SetMarkerColor(MYCOLZ[i])
        h.Draw("hist same")
        h.Draw("p same")
    b.Draw("axis same")
    
    ll = ROOT.TLatex()
    ll.SetNDC()
    ll.SetTextFont(72)
    ll.DrawLatex(.12,.89,"#scale[0.9]{CMS}")
    ll.SetTextFont(42)
    ll.DrawLatex(.165,.89,"#scale[0.9]{Internal}")

    #create and fill legend
    leg = ROOT.TLegend(.12,.81, 0.52,.96)
    leg.SetTextFont(42)
    leg.SetHeader("")
    leg.SetNColumns(4)
    leg.AddEntry(hists[0],  order[0],"l")
    leg.AddEntry(hists[1],  order[1],"l")
    leg.AddEntry(hists[2],  order[2],"l")
    leg.AddEntry(hists[3],  order[3],"l")

    #draw the legend
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.Draw()
    
    c.SaveAs(plotname+".pdf")



def DrawHistsLinkClk2D(hists, plotname,
                   ytitle="# objects"):
    #configure the canvas, pad
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("canv","",750,750)
    pad = ROOT.TPad("pad", "pad", .005, .01, .995, .995)
    pad.Draw()
    pad.cd()

    ROOT.gPad.SetLeftMargin(0.17)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetTopMargin(0.05)

    #draw this guy first
    b = hists[0].Clone("tmp")
    b.Reset()
    for h in hists: b.Add(h)
    #    b.Divide(b)
    b.Draw()
    #b.Draw("colz same")
    
    b.GetXaxis().SetTitleOffset(1.2)
    b.GetYaxis().SetTitleOffset(1.5)
    b.GetYaxis().SetLabelSize(0.04)
    b.GetYaxis().SetTitleSize(0.04)
#    b.GetYaxis().SetTitle(ytitle)
    b.Draw()

    # col=[ROOT.kBlack,ROOT.kBlue,ROOT.kRed,ROOT.kViolet]
    # for i,h in enumerate(hists):
    #     h.SetLineWidth(2)
    #     h.SetLineColor(col[i])
    #     h.SetMarkerStyle(20)
    #     h.SetMarkerColor(col[i])
    #     h.Draw("hist same")
    #     h.Draw("p same")
    
    # ll = ROOT.TLatex()
    # ll.SetNDC()
    # ll.SetTextFont(72)
    # ll.DrawLatex(.20,.89,"#scale[0.9]{CMS}")
    # ll.SetTextFont(42)
    # ll.DrawLatex(.30,.89,"#scale[0.9]{Internal}")

    # #create and fill legend
    # leg = ROOT.TLegend(.20,.76, 0.85,.91)
    # leg.SetTextFont(42)
    # leg.SetHeader("")
    # leg.SetNColumns(4)
    # leg.AddEntry(hists[0],  order[0],"l")
    # leg.AddEntry(hists[1],  order[1],"l")
    # leg.AddEntry(hists[2],  order[2],"l")
    # leg.AddEntry(hists[3],  order[3],"l")

    # #draw the legend
    # leg.SetFillStyle(0)
    # leg.SetFillColor(0)
    # leg.SetBorderSize(0)
    # leg.Draw()
    
    c.SaveAs(plotname+".pdf")
