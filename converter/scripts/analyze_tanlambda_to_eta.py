from helper import *

# configurables
testname="tanlambda_to_eta"
var_string = "tl tl_calo z0 tl_hw tl_calo_hw z0_hw"
var_string = "tanlambda eta tanlambda_hw eta_hw eta_hwf tab_size"
tab_sizes = [6,8,10,12]

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
    book(h,"eta_ref" ,60,-3,3,";eta(ref)", pf=pf)
    book(h,"eta_hw"  ,60,-3,3,";eta(HW)", pf=pf)
    # diff
    book(h,"eta_diff",80,-1,1,";eta(HW)-eta(ref)", pf=pf)
    # profile
    bookp(h,"p_eta_diff" ,60,-3,3,-1,1,";ETA(ref);ETA(HW)-ETA(ref)", pf=pf)

    # fill
    draw(t,"eta","eta_ref", pf=pf, cut=cut)
    draw(t,"eta_hwf","eta_hw", pf=pf, cut=cut)
    draw(t,"eta_hwf-eta","eta_diff", pf=pf, cut=cut)
    draw(t,"eta_hwf-eta:eta","p_eta_diff", pf=pf, cut=cut)
    
    # checks
    bookp(h,"p_tanlambda_eta" ,80,-8,8,-3,3,";tan(lambda)(ref);ETA(ref)", pf=pf)
    bookp(h,"p_tanlambda_eta_hw" ,80,-8,8,-3,3,";tan(lambda)(HW);ETA(HW)", pf=pf)
    bookp(h,"p_tanlambda_eta_hw_int" ,80,-8,8,-550,550,";tan(lambda)(HW);ETA(HW int)", pf=pf)
    draw(t,"eta:tanlambda","p_tanlambda_eta", pf=pf, cut=cut)
    draw(t,"eta_hwf:tanlambda_hw","p_tanlambda_eta_hw", pf=pf, cut=cut)
    draw(t,"eta_hw:tanlambda_hw","p_tanlambda_eta_hw_int", pf=pf, cut=cut)
    
            
outname="hists/"+testname+".root"
f = ROOT.TFile(outname,"recreate")
for x in h: h[x].Write()
f.Close()
