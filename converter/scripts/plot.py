import ROOT
from collections import OrderedDict
COLZ=[1,2,4,6,3]

cfg = {
    "tanlambda_to_eta":{
        "p_eta_diff": {
            "norm": 0, "logy":0, "ymax":0.4, "ymin":-0.3, "rms":0,
            "tab_sizes" : [6,8,10,12], "leg" : 1
        },
    },
    "pt_inversion":{
        "p_pt_reldiff": {
            "norm": 0, "logy":0, "ymax":0.8, "ymin":-0.6, "rms":0,
            "tab_sizes" : [6,7,8], "leg" : 1
        },
    },
    "prop_phi":{
        "p_dphi_diff": {
            "norm": 0, "logy":0, "ymax":0.01, "ymin":-0.01, "rms":0,
            "tab_sizes" : [8], "leg" : 0
        },
    },
    "prop_tanlambda":{
        "p_tl_tl_calo_diff": {
            "norm": 0, "logy":0, "ymax":0.3, "ymin":-0.3, "rms":0,
            "tab_sizes" : [], "leg" : 0
        },
        "p_z0_tl_calo_diff": {
            "norm": 0, "logy":0, "ymax":0.2, "ymin":-0.2, "rms":0,
            "tab_sizes" : [], "leg" : 0
        },
    },
    "resolution":{
        "p_pt_pt_err_diff": {
            "norm": 0, "logy":0, "ymax":0.1, "ymin":-0.4, "rms":0,
            "tab_sizes" : [], "leg" : 0
        },
        "p_eta_pt_err_diff": {
            "norm": 0, "logy":0, "ymax":0.1, "ymin":-0.6, "rms":0,
            "tab_sizes" : [], "leg" : 0
        },
    },
}
def main():
    for algo in cfg:
        fname="hists/"+algo+".root"
        f = ROOT.TFile(fname,"read")
        for hname in cfg[algo]:
            opts = cfg[algo][hname]
            plot(f, algo, hname, opts)
        f.Close()

    
def plot(f, algo, hname, opts):
    # get hists and process if needed
    tab_sizes = opts["tab_sizes"]
    hists = OrderedDict()
    if len(tab_sizes):
        for n in tab_sizes:
            hists[n] = f.Get("{}_tab{}".format(hname,n))
    else:
        hists[1] = f.Get(hname)
        tab_sizes=[1]

    if opts["norm"]:
        for h in hists.values():
            if h.Integral(): h.Scale(1./h.Integral())

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
    if opts["logy"]:
        pad.SetLogy()

    b = hists[tab_sizes[0]].Clone("tmp")
    b.Reset()
    b.SetLineColor(1)
    b.Draw()
    if "ymax" in opts:
        b.SetMaximum( opts["ymax"] )
    if "ymin" in opts:
        b.SetMinimum( opts["ymin"] )
    
    b.GetXaxis().SetTitleOffset(1.2)
    b.GetYaxis().SetTitleOffset(1.5)
    b.GetYaxis().SetLabelSize(0.04)
    b.GetYaxis().SetTitleSize(0.04)
    b.Draw("")

    for ih,n in enumerate(tab_sizes):
        h = hists[n]
        col = COLZ[ih]
        h.SetLineWidth(2)
        h.SetMarkerSize(0.8)
        h.SetMarkerStyle(20)
        h.SetLineColor(col)
        h.SetMarkerColor(col)        
        #h.Draw("hist same")
        h.Draw("pe same")
        # h.Draw("same plc pmc")

    
    ll = ROOT.TLatex()
    ll.SetNDC()
    ll.SetTextFont(72)
    ll.DrawLatex(.20,.89,"#scale[0.9]{CMS}")
    ll.SetTextFont(42)
    ll.DrawLatex(.30,.89,"#scale[0.9]{Internal}")

    #create and fill legend
    if opts["leg"]:
        leg = ROOT.TLegend(.60,.70, 0.95,.95)
        leg.SetTextFont(42)
        leg.SetHeader("")
        leg.SetNColumns(1)
        for n in tab_sizes:
            suf=""
            if opts["rms"]:
                suf = " (RMS {:.3f})".format(hists[n].GetRMS())
            #leg.AddEntry(hists[n], "{} input{}{}".format(n,"s"*(n>1),suf),"pel")
            leg.AddEntry(hists[n], "LUT bits: {}".format(n),"pel")
        
        #draw the legend
        leg.SetFillStyle(0)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.Draw()
    
    c.SaveAs("output/test_"+algo+"_"+hname+".pdf")
    del c
        
# opts=[]
# for hname in opts:
#     # get hists and process if needed
#     hists = {n: f.Get("{}_tab{}".format(hname,n)) for n in tab_sizes}
    
#     if opts[hname]["norm"]:
#         for h in hists.values():
#             if h.Integral(): h.Scale(1./h.Integral())

#     #configure the canvas, pad
#     ROOT.gStyle.SetOptStat(0)
#     c = ROOT.TCanvas("canv","",750,750)
#     pad = ROOT.TPad("pad", "pad", .005, .01, .995, .995)
#     pad.Draw()
#     pad.cd()

#     ROOT.gPad.SetLeftMargin(0.17)
#     ROOT.gPad.SetRightMargin(0.05)
#     ROOT.gPad.SetBottomMargin(0.1)
#     ROOT.gPad.SetTopMargin(0.05)
#     if opts[hname]["logy"]:
#         pad.SetLogy()

#     b = hists[tab_sizes[0]].Clone("tmp")
#     b.Reset()
#     b.Draw()
#     if "ymax" in opts[hname]:
#         b.SetMaximum( opts[hname]["ymax"] )
#     if "ymin" in opts[hname]:
#         b.SetMinimum( opts[hname]["ymin"] )
    
#     b.GetXaxis().SetTitleOffset(1.2)
#     b.GetYaxis().SetTitleOffset(1.5)
#     b.GetYaxis().SetLabelSize(0.04)
#     b.GetYaxis().SetTitleSize(0.04)
#     b.Draw("")

#     for ih,n in enumerate(tab_sizes):
#         h = hists[n]
#         col = COLZ[ih]
#         h.SetLineWidth(2)
#         h.SetMarkerSize(0.8)
#         h.SetMarkerStyle(20)
#         h.SetLineColor(col)
#         h.SetMarkerColor(col)        
#         #h.Draw("hist same")
#         h.Draw("pe same")
#         # h.Draw("same plc pmc")

    
#     ll = ROOT.TLatex()
#     ll.SetNDC()
#     ll.SetTextFont(72)
#     ll.DrawLatex(.20,.89,"#scale[0.9]{CMS}")
#     ll.SetTextFont(42)
#     ll.DrawLatex(.30,.89,"#scale[0.9]{Internal}")

#     #create and fill legend
#     leg = ROOT.TLegend(.60,.70, 0.95,.95)
#     leg.SetTextFont(42)
#     leg.SetHeader("")
#     leg.SetNColumns(1)
#     for n in tab_sizes:
#         suf=""
#         if opts[hname]["rms"]:
#             suf = " (RMS {:.3f})".format(hists[n].GetRMS())
#         #leg.AddEntry(hists[n], "{} input{}{}".format(n,"s"*(n>1),suf),"pel")
#         leg.AddEntry(hists[n], "LUT bits: {}".format(n),"pel")

#     #draw the legend
#     leg.SetFillStyle(0)
#     leg.SetFillColor(0)
#     leg.SetBorderSize(0)
#     leg.Draw()
    
#     c.SaveAs("output/test_"+hname+".pdf")
#     del c
        
if __name__ == "__main__":
    # execute only if run as a script
    main()
