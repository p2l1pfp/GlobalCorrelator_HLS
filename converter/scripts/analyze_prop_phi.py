from helper import *

# configurables
testname="prop_phi"
var_string = "pt_inv dphi pt_inv_hw pt_inv_hwf dphi_hw dphi_hwf tab_size"
tab_sizes = [8,10]

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
    book(h,"pt_inv_ref" ,55,0,0.55,";1/pt(ref)", pf=pf)
    book(h,"pt_inv_hw"  ,55,0,0.55,";1/pt(HW)", pf=pf)
    book(h,"dphi_ref" ,50,0,0.5,";dphi(ref)", pf=pf)
    book(h,"dphi_hw"  ,50,0,0.5,";dphi(HW)", pf=pf)
    # diff
    book(h,"dphi_diff",80,-0.5,0.5,";dphi(HW)-dphi(ref)", pf=pf)
    # profile
    bookp(h,"p_dphi_diff" ,55,0,0.55,-1,1,";1/pt;dphi(HW)-dphi(ref)", pf=pf)

    # fill
    draw(t,"pt_inv","pt_inv_ref", pf=pf, cut=cut)
    draw(t,"pt_inv_hw","pt_inv_hw", pf=pf, cut=cut)
    draw(t,"dphi","dphi_ref", pf=pf, cut=cut)
    draw(t,"dphi_hwf","dphi_hw", pf=pf, cut=cut)
    draw(t,"dphi_hwf-dphi","dphi_diff", pf=pf, cut=cut)
    draw(t,"dphi_hwf-dphi:pt_inv","p_dphi_diff", pf=pf, cut=cut)
    
            
outname="hists/"+testname+".root"
f = ROOT.TFile(outname,"recreate")
for x in h: h[x].Write()
f.Close()
