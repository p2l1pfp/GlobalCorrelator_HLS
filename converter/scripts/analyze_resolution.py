from helper import *

# configurables
testname="resolution"
var_string = "pt pt_err eta pt_hw pt_hwf pt_err_hw pt_err_hwf eta_hw eta_hwf"
tab_sizes = [8]

# tree
fname="../tests/results/test_"+testname+".txt"
t = ROOT.TTree("t","")
var_string=var_string.replace(' ',':')
t.ReadFile(fname,var_string)
#hists
h=OrderedDict()    

for ts in tab_sizes:
    # basic
    pf = ""
    cut= ""
    book(h,"pt_ref" ,60,0,120,";pt(ref) [GeV]", pf=pf)
    book(h,"pt_hw"  ,60,0,120,";pt(HW) [GeV]", pf=pf)
    book(h,"pt_err_ref" ,60,0,120,";pt err(ref)", pf=pf)
    book(h,"pt_err_hw"  ,60,0,120,";pt err(HW)", pf=pf)
    book(h,"eta_ref" ,60,-3,3,";eta(ref)", pf=pf)
    book(h,"eta_hw"  ,60,-3,3,";eta(HW)", pf=pf)
    # diff
    book(h,"pt_err_diff",80,-20,20,";pt err(HW) - pt err(ref)", pf=pf)
    # profile
    bookp(h,"p_pt_pt_err_diff" ,60,0,120,-10,10,";pT [GeV];pt err(HW)-pt err(ref) [GeV]", pf=pf)
    bookp(h,"p_eta_pt_err_diff" ,60,-3,3,-10,10,";eta;pt err(HW)-pt err(ref) [GeV]", pf=pf)

    # fill
    draw(t,"pt","pt_ref", pf=pf, cut=cut)
    draw(t,"pt_hwf","pt_hw", pf=pf, cut=cut)
    draw(t,"pt_err","pt_err_ref", pf=pf, cut=cut)
    draw(t,"pt_err_hwf","pt_err_hw", pf=pf, cut=cut)
    draw(t,"eta","eta_ref", pf=pf, cut=cut)
    draw(t,"eta_hwf","eta_hw", pf=pf, cut=cut)

    draw(t,"pt_err_hwf-pt_err","pt_err_diff", pf=pf, cut=cut)
    draw(t,"pt_err_hwf-pt_err:pt","p_pt_pt_err_diff", pf=pf, cut=cut)
    draw(t,"pt_err_hwf-pt_err:eta","p_eta_pt_err_diff", pf=pf, cut=cut)



outname="hists/"+testname+".root"
f = ROOT.TFile(outname,"recreate")
for x in h: h[x].Write()
f.Close()
