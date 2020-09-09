#
import ROOT
from collections import OrderedDict
PI=ROOT.TMath.Pi() #3.141593

# helpers 
def book(h,name,n,a,b,title="", pf=""):
    if pf: name = name+"_"+pf
    h[name]=ROOT.TH1F(name,title,n,a,b)
def book2(h,name,nx,ax,bx,ny,ay,by,title="", pf=""):
    if pf: name = name+"_"+pf
    h[name]=ROOT.TH2F(name,title,nx,ax,bx,ny,ay,by)
def bookp(h,name,nx,ax,bx,ay,by,title="",err="s", pf=""):
    if pf: name = name+"_"+pf
    h[name]=ROOT.TProfile(name,title,nx,ax,bx,ay,by,err)
def bookp2(h,name,nx,ax,bx,ny,ay,by,az,bz,title="",err="s", pf=""):
    if pf: name = name+"_"+pf
    h[name]=ROOT.TProfile2D(name,title,nx,ax,bx,ny,ay,by,az,bz,err)
def draw(t,var,hname, pf="", cut=""):
    if pf: hname = hname+"_"+pf
    t.Draw(var+">>"+hname,cut)

