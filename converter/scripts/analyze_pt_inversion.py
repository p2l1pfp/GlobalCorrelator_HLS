from helper import *

# configurables
testname="pt_inversion"
var_string = "pt pt_inv pt_hw pt_hwf pt_inv_hw tab_size max_bits"
tab_sizes = [6,7,8]

# tree
fname="../tests/results/test_"+testname+".txt"
t = ROOT.TTree("t","")
var_string=var_string.replace(' ',':')
t.ReadFile(fname,var_string)
#hists
h=OrderedDict()    

for ts in tab_sizes:
    # basic
    pf = "tab{}".format(ts)
    cut= "tab_size=={}".format(ts)
    book(h,"pt_ref" ,80,0,200,";pt(ref)", pf=pf)
    book(h,"pt_hw"  ,80,0,200,";pt(HW)", pf=pf)
    book(h,"pt_inv_ref" ,55,0,0.55,";1/pt(ref)[1/GeV]", pf=pf)
    book(h,"pt_inv_hw"  ,55,0,0.55,";1/pt(HW)[1/GeV]", pf=pf)
    # diff
    book(h,"pt_diff",80,-20,20,";pt(HW)-pt(ref)", pf=pf)
    book(h,"pt_inv_diff",80,-0.1,0.1,";pt(HW)-pt(ref)", pf=pf)
    # profile
    bookp(h,"p_pt_diff" ,70,0,140,-50,50,";pT[1/GeV];pt(HW)-pt(ref)", pf=pf)
    bookp(h,"p_pt_reldiff" ,70,0,140,-50,50,";pT[GeV];pt(HW)/pt(ref)-1", pf=pf)
    bookp(h,"p_pt_inv_diff",70,0,0.55,-50,50,";1/pT;pt(HW)-pt(ref)", pf=pf)

    # fill
    draw(t,"pt","pt_ref", pf=pf, cut=cut)
    draw(t,"pt_hwf","pt_hw", pf=pf, cut=cut)
    draw(t,"pt_inv","pt_inv_ref", pf=pf, cut=cut)
    draw(t,"pt_inv_hw","pt_inv_hw", pf=pf, cut=cut)
    draw(t,"pt_hwf-pt","pt_diff", pf=pf, cut=cut)
    draw(t,"pt_inv_hw-pt_inv","pt_inv_diff", pf=pf, cut=cut)

    draw(t,"pt_hwf-pt:pt","p_pt_diff", pf=pf, cut=cut)
    draw(t,"pt_hwf/pt-1:pt","p_pt_reldiff", pf=pf, cut=cut)
    draw(t,"pt_hwf-pt:pt_inv","p_pt_inv_diff", pf=pf, cut=cut)

            
outname="hists/"+testname+".root"
f = ROOT.TFile(outname,"recreate")
for x in h: h[x].Write()
f.Close()
