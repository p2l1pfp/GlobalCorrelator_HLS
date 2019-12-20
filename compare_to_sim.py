import sys
import os
import argparse
import subprocess
import time

NEVENTS=10
#NEVENTS=1
NLARGE=1
NSMALL=18

#CLKMAX=150
#CLKMAX=108
TMUX_OUT=18
CLKMAX=TMUX_OUT*NEVENTS
CLKMAX=140 # 134 lines in sim out

NTRACK=22
NEM=13
NCALO=15
NMU=2
#
NPHOTON=NEM
NSELCALO=10

#NTRACK+NPHOTON+NSELCALO+NMU=47

ORDERING = [(0,0), (0,2), (0,1), (0,3), (1,2),
            (0,4), (1,3), (0,5), (1,4), (0,6),
            (1,5), (0,7), (1,6), (0,8), (1,7), 
            (1,0), (1,8), (1,1)]

class Object:
    def __init__(self, 
                 pt   = 0,
                 pte  = 0,
                 eta  = 0,
                 phi  = 0,
                 z0   = 0,
                 qual = 0):
        self.pt   = pt  
        self.pte  = pte 
        self.eta  = eta 
        self.phi  = phi 
        self.z0   = z0  
        self.qual = qual

class Region:
    def __init__(self, n_sub=0):
        #add subregions if applicable
        self.n_subregions = n_sub
        self.subregions = []
        for i in range(self.n_subregions):
            smallR = Region()
            smallR.reg_type = "small"
            self.subregions.append(smallR)

        self.reg_type  = None
        self.eta_index = -1
        self.phi_index = -1
        self.tracks    = []
        self.ems       = []
        self.calos     = []
        self.mus       = []
        self.ID        = -1 #unique input ID
        self.track_clks= {} # record clock assoc w each track
        self.ntracks = 0
        self.nems    = 0
        self.ncalos  = 0
        self.nmus    = 0
        self.pf_ch=[]
        self.pf_ph=[]
        self.pf_ne=[]
        self.pf_mu=[]
        self.npf_ch=0
        self.npf_ph=0
        self.npf_ne=0
        self.npf_mu=0
    def counts(self): return (
            self.ntracks,
            self.nems   ,
            self.ncalos ,
            self.nmus   ,)
    def countsPFPUPPI(self): return (
            self.npf_ch,
            self.npf_ph,
            self.npf_ne,
            self.npf_mu,)
    def resize_subregions(self, n_resize): 
        if n_resize <= self.n_subregions: return
        for i in range(n_resize-self.n_subregions):
            smallR = Region()
            smallR.reg_type = "small"
            self.subregions.append(smallR)        
        self.n_subregions = n_resize
        return

class Event:
    def __init__(self, n_large = 1, n_small = 1):
        self.n_small_regions = n_small
        self.n_big_regions   = n_large
        self.large_regions   = []
        for i in range(self.n_big_regions):
            bigR = Region(n_sub = self.n_small_regions)
            bigR.reg_type = "large"
            self.large_regions.append(bigR)
        self.ntracks = 0
        self.nems    = 0
        self.ncalos  = 0
        self.nmus    = 0
        self.npf_ch=0
        self.npf_ph=0
        self.npf_ne=0
        self.npf_mu=0
        self.obj_to_regions={}
        self.obj_to_regions_filled=False
        self.obj_to_regions_pf={}
        self.obj_to_regions_filled_pf=False
        self.index   = -1
    def counts(self): return (
            self.ntracks,
            self.nems   ,
            self.ncalos ,
            self.nmus   ,)
    def countsPFPUPPI(self): return (
            self.npf_ch,
            self.npf_ph,
            self.npf_ne,
            self.npf_mu,)
    def allPF(self,tag):
        x=[]
        for l in self.large_regions: 
            for s in l.subregions: 
                x += getattr(s,"pf_"+tag)
                # if tag=="ch": x += s.pf_ch
                # if tag=="ph": x += s.pf_ph
                # if tag=="ne": x += s.pf_ne
                # if tag=="mu": x += s.pf_mu
        return x
    def allTracks(self):
        x=[]
        for l in self.large_regions: 
            for s in l.subregions: x += s.tracks
        return x
    def allEMs(self):
        x=[]
        for l in self.large_regions: 
            for s in l.subregions: x += s.ems
        return x
    def allCalos(self):
        x=[]
        for l in self.large_regions: 
            for s in l.subregions: x += s.calos
        return x
    def allMuons(self):
        x=[]
        for l in self.large_regions: 
            for s in l.subregions: x += s.mus
        return x
    def findRegions(self,x):
        if(self.obj_to_regions_filled==False):
            for li,l in enumerate(self.large_regions): 
                for si,s in enumerate(l.subregions): 
                    for obj in s.tracks+s.calos+s.ems+s.mus: 
                        # if obj == 0x793fdfa5000c000e:
                        #     print(li,si)
                        if obj in self.obj_to_regions: self.obj_to_regions[obj] += [ (li,si) ]
                        else: self.obj_to_regions[obj] = [ (li,si) ]
            self.obj_to_regions_filled=True
            return self.obj_to_regions[x]
        # if x == 0x793fdfa5000c000e:
        #     if x in self.obj_to_regions: 
        #         print (self.obj_to_regions[x])
        #     else: print ((-1,-1))
        if x in self.obj_to_regions: return self.obj_to_regions[x]
        else: return (-1,-1)
    def findRegionsPF(self,x):
        if(self.obj_to_regions_filled_pf==False):
            for li,l in enumerate(self.large_regions): 
                for si,s in enumerate(l.subregions): 
                   for obj in s.pf_ch+s.pf_ph+s.pf_ne+s.pf_mu:
                        if obj in self.obj_to_regions_pf: self.obj_to_regions_pf[obj] += [ (li,si) ]
                        else: self.obj_to_regions_pf[obj] = [ (li,si) ]
            self.obj_to_regions_filled_pf=True
            return self.obj_to_regions_pf[x]
        if x in self.obj_to_regions_pf: return self.obj_to_regions_pf[x]
        else: return (-1,-1)

    # def AddLargeRegion(self):
    #     pass #self.small_regions = []

def SelectBits(x, nbits, shift):
    return ((2**nbits-1 << shift) & x) >> shift
def BitsToInt(x,nbits):
    if x < 2**(nbits-1): return x
    else: return x - 2**nbits

def SelectBitsInt(x, nbits, shift):
    return BitsToInt(SelectBits(x, nbits, shift),nbits)

def GetTrackParams(x):
    pt   =  SelectBitsInt(x,16, 0)
    pte  =  SelectBitsInt(x,16,16)
    eta  =  SelectBitsInt(x,10,32)#-2**10 # allow negative
    phi  =  SelectBitsInt(x,10,42)
    z0   =  SelectBitsInt(x,10,52)
    qual =  SelectBitsInt(x, 1,62)
    return pt, pte, eta, phi, z0, qual
    #print("track: pt={}, pte={}, eta={}, phi={}, z0={}, qual={} ({})".format(pt, pte, eta, phi, z0, qual, x))

def GetCaloParams(x):
    pt   =  SelectBitsInt(x,16, 0)
    empt =  SelectBitsInt(x,16,16)
    eta  =  SelectBitsInt(x,10,32)#-2**10 # allow negative
    phi  =  SelectBitsInt(x,10,42)
    isEM =  SelectBitsInt(x, 1,52)
    return pt, empt, eta, phi, isEM
    #print("calo: pt={}, empt={}, eta={}, phi={}, isEM={} ({})".format(pt, empt, eta, phi, isEM, x))

def GetEMParams(x):
    pt   =  SelectBitsInt(x,16, 0)
    pte  =  SelectBitsInt(x,16,16)
    eta  =  SelectBitsInt(x,10,32)#-2**10 # allow negative
    phi  =  SelectBitsInt(x,10,42)
    return pt, pte, eta, phi
    #print("EM: pt={}, pte={}, eta={}, phi={} ({})".format(pt, pte, eta, phi, x))
def GetPFChParams(x):
    # (2 null)+(10 z0)+(10 phi)+(10 eta)+(13 NULL)+(3 ID)+(16 pt)
    pt   =  SelectBitsInt(x,16, 0)
    ID   =  SelectBits   (x, 3,16)
    eta  =  SelectBitsInt(x,10,32)
    phi  =  SelectBitsInt(x,10,42)
    z0   =  SelectBitsInt(x,10,52)
    return ID, pt, z0, eta, phi
def GetPFNeParams(x):
    # (9 null)+(3 ID)++(10 phi)+(10 eta)+(16 Puppi)+(16 pt)
    pt   =  SelectBitsInt(x,16, 0)
    puppi=  SelectBitsInt(x,16,16)
    eta  =  SelectBitsInt(x,10,32)
    phi  =  SelectBitsInt(x,10,42)
    ID   =  SelectBits   (x, 3,52)
    return puppi, pt, ID, eta, phi

def GetPtEtaPhi(x,tag):
    if tag=="track" or tag=="tk" or tag=="mu":
        pt, pte, eta, phi, z0, qual = GetTrackParams(x)
        return pt, eta, phi
    if tag=="calo":
        pt, empt, eta, phi, isEM = GetCaloParams(x)
        return pt, eta, phi
    if tag=="em":
        pt, pte, eta, phi = GetEMParams(x)
        return pt, eta, phi
    if tag=="PFch":
        ID, pt, z0, eta, phi = GetPFChParams(x)
        return pt, eta, phi
    if tag=="PFne":
        puppi, pt, ID, eta, phi = GetPFNeParams(x)
        return pt, eta, phi
    return None

def GetEmulationData(parser):
    '''
    Emulation data is in a single file.
    Each line has 48=(15+15+15+2+1) entries
    Corresponding to track,em,calo,mu,pv objects.
    The number of lines is equal to TMUX_OUT*NTEST
    This may be 18*6=108
    '''
    
    # from ryan. this is the order of (eta,phi) regions output by the emulator
    # 0,18,36,.. are (0,0), 1,19,37,... are (0,2) and so on...

    tracks=[]
    ems=[]
    calos=[]
    mus=[]

    # local vars
    EM_NEVENTS = NEVENTS if (NEVENTS>0) else 6
    #EM_NEVENTS=21
    EM_NLARGE=NLARGE if (NLARGE>0) else 1
    EM_NSMALL=NSMALL if (NSMALL>0) else 18

    evts = []

    # setup for the first region
    evt = Event(EM_NLARGE, EM_NSMALL)
    evt.index=0

    with open(parser.emulator_output,'r') as f:
        ctr = 0
        sreg_ctr = 0
        lreg_ctr = 0
        evt_ctr  = 0
        for l in f:

            # retrieve the small region
            sr = evt.large_regions[lreg_ctr].subregions[sreg_ctr]

            # extract data and fill
            arr = list(map(lambda x : int(x,16), l.split()))
            step = arr[0] # unused
            sr.tracks = list(filter(lambda x: x, arr[1:NTRACK+1]                               ))
            sr.ems    = list(filter(lambda x: x, arr[1+NTRACK:NTRACK+NEM+1]                    ))
            sr.calos  = list(filter(lambda x: x, arr[1+NTRACK+NEM:NTRACK+NEM+NCALO+1]          ))
            sr.mus    = list(filter(lambda x: x, arr[1+NTRACK+NEM+NCALO:NTRACK+NEM+NCALO+NMU+1]))
            #remove 0 pt
            sr.tracks = list(filter(lambda x: GetPtEtaPhi(x,"track")[0], sr.tracks))
            sr.ems    = list(filter(lambda x: GetPtEtaPhi(x,"em")   [0], sr.ems   ))
            sr.calos  = list(filter(lambda x: GetPtEtaPhi(x,"calo") [0], sr.calos ))
            sr.mus    = list(filter(lambda x: GetPtEtaPhi(x,"mu")   [0], sr.mus   ))
            #remove eta=32 TEMPORARY TEST!!
            # sr.tracks = list(filter(lambda x: GetPtEtaPhi(x,"track")[1]!=32, sr.tracks))
            # sr.ems    = list(filter(lambda x: GetPtEtaPhi(x,"em")   [1]!=32, sr.ems   ))
            # sr.calos  = list(filter(lambda x: GetPtEtaPhi(x,"calo") [1]!=32, sr.calos ))
            # sr.mus    = list(filter(lambda x: GetPtEtaPhi(x,"mu")   [1]!=32, sr.mus   ))
            sr.ID     = step

            #if 0x0463c00d000a000a in sr.tracks: print("We found it",ctr,sreg_ctr)
            #if 0x793fdfa5000c000e in sr.tracks: print("We found it",ctr,sreg_ctr)

            sr.ntracks  += len(sr.tracks)
            sr.nems     += len(sr.ems   )
            sr.ncalos   += len(sr.calos )
            sr.nmus     += len(sr.mus   )
            evt.large_regions[lreg_ctr].ntracks  += len(sr.tracks)
            evt.large_regions[lreg_ctr].nems     += len(sr.ems   )
            evt.large_regions[lreg_ctr].ncalos   += len(sr.calos )
            evt.large_regions[lreg_ctr].nmus     += len(sr.mus   )
            evt.ntracks += len(sr.tracks)
            evt.nems    += len(sr.ems   )
            evt.ncalos  += len(sr.calos )
            evt.nmus    += len(sr.mus   )


            # increment counters
            # (add support for large region ctr later)
            ctr = ctr+1
            sreg_ctr = ctr % EM_NSMALL
            if ctr and (ctr % EM_NSMALL)==0: 
                evt_ctr = evt_ctr+1
                # create new event
                evts.append( evt )
                evt = Event(EM_NLARGE, EM_NSMALL)
                evt.index=evt_ctr
                # large_reg = Region(n_sub=EM_NSMALL)
                # evt.large_regions.append( large_reg )


    # add empty events for now
    while len(evts) < EM_NEVENTS:
        print("adding an empty emulator event")
        evt = Event(EM_NLARGE, EM_NSMALL)
        evts.append( evt )
        
    print( "EM_NEVENTS",EM_NEVENTS,len(evts))
    #'serialized' list of 18 tmux regions * however many events...
    #return tracks, ems, calos, mus
    return evts

def GetSimulationData(parser):
    '''
    Simulation data is stored across multiple files,
    each corresponding to a time series of single-object outputs.
    #files is #em+#calo+#track, which may equal 15+20+25=60 --> 22+13+15+2=52?
    Each line corresonds to a 'step' of the output and...
    TODO add comment on length of each file.
    '''

    # local vars
    SIM_NEVENTS = NEVENTS if (NEVENTS>0) else 6
    SIM_NLARGE  = NLARGE if (NLARGE>0) else 1
    SIM_NSMALL  = NSMALL if (NSMALL>0) else 18
    SIM_NSMALL  = 18
    #SIM_NSMALL  = 80

    evts = []
    for i in range(SIM_NEVENTS):
        e = Event(SIM_NLARGE, SIM_NSMALL) 
        e.index = i
        evts.append( e )

    sreg_ctr = 0
    lreg_ctr = 0
    evt_ctr  = 0

    for ii in range(NEM+NCALO+NTRACK+NMU):
        with open("{}/sim_HLS_input_object_{}.dat".format(parser.sim_output_dir,ii),'r') as f:
            ctr=0
            for l in f:
                #print('ctr = ',ctr)
                val = int(l,base=16)
                # ignore empty lines
                sreg_ctr = ctr % SIM_NSMALL
                evt_ctr  = int(ctr / SIM_NSMALL)
                if evt_ctr>=SIM_NEVENTS: continue
                # (assume only one large region for now)
                #print(evt_ctr,lreg_ctr,sreg_ctr)
                sr = evts[evt_ctr].large_regions[lreg_ctr].subregions[sreg_ctr]
                if(ii==0): sr.ID=ctr
                if val:
                    if ii<NEM:
                        #if GetPtEtaPhi(val,"em")[1]==32: continue # TODO TEMP TEST
                        sr.ems.append(val)
                        sr.nems  += 1
                        evts[evt_ctr].nems  += 1
                    elif ii<NEM+NCALO: 
                        #if GetPtEtaPhi(val,"calo")[1]==32: continue # TODO TEMP TEST
                        if GetPtEtaPhi(val,"calo")[1]!=32: 
                            sr.calos.append(val)
                            sr.ncalos  += 1
                            evts[evt_ctr].ncalos  += 1
                    elif ii<NEM+NCALO+NTRACK:
                        #if GetPtEtaPhi(val,"tk")[1]==32: continue # TODO TEMP TEST
                        #if 0x793fdfa5000c000e ==val: print("We found it",ctr,sreg_ctr)
                        sr.tracks.append(val)
                        sr.track_clks[val]=ctr
                        sr.ntracks  += 1
                        evts[evt_ctr].ntracks  += 1
                    elif ii<NEM+NCALO+NTRACK+NMU:
                        #if GetPtEtaPhi(val,"mu")[1]==32: continue # TODO TEMP TEST
                        sr.mus.append(val)
                        sr.nmus  += 1
                        evts[evt_ctr].nmus  += 1
                ctr += 1

    return evts

def GetEmulationDataPFPUPPI(evts, parser):
    '''
    Get sequence of regional info out
    Sequence should be same as for regionizer
    (see 'ORDERING' list above) = [(0,0), (0,2), ...]

    Emulation data (after layer-1 PF+PUPPI alg) is in a single file.
    Each line has 47 = NTRACK + NNEUTRALS(=NPHOTON+NSELCALO) + NMU + ZPV = 22+(13+10)+2+1 entries
    The number of lines is equal to TMUX_OUT*NTEST
    This may be 18*6=108 if running 6 events
    '''
    
    # local vars
    EM_NEVENTS = NEVENTS if (NEVENTS>0) else 6
    #EM_NEVENTS=21
    EM_NLARGE=NLARGE if (NLARGE>0) else 1
    EM_NSMALL=NSMALL if (NSMALL>0) else 18


    with open(parser.emulator_output_pfpuppi,'r') as f:
        ctr = 0
        sreg_ctr = 0
        lreg_ctr = 0
        evt_ctr  = 0
        for l in f:

            # retrieve the small region
            if evt_ctr >= len(evts): 
                print("Trying to fill an emulation event that doesn't exist!")
                return
            evt=evts[evt_ctr]
            sr = evt.large_regions[lreg_ctr].subregions[sreg_ctr]

            # extract data and fill
            arr = list(map(lambda x : int(x,16), l.split()))
            step = arr[0] # unused
            sr.pf_ch = list(filter(lambda x: x, arr[1:NTRACK+1]                                             ))
            sr.pf_ph = list(filter(lambda x: x, arr[1+NTRACK:NTRACK+NPHOTON+1]                              ))
            sr.pf_ne = list(filter(lambda x: x, arr[1+NTRACK+NPHOTON:NTRACK+NPHOTON+NSELCALO+1]             ))
            sr.pf_mu = list(filter(lambda x: x, arr[1+NTRACK+NPHOTON+NSELCALO:NTRACK+NPHOTON+NSELCALO+NMU+1]))
            #remove 0 pt
            sr.pf_ch = list(filter(lambda x: GetPtEtaPhi(x,"PFch")[0], sr.pf_ch))
            sr.pf_ph = list(filter(lambda x: GetPtEtaPhi(x,"PFne")[0], sr.pf_ph))
            sr.pf_ne = list(filter(lambda x: GetPtEtaPhi(x,"PFne")[0], sr.pf_ne))
            sr.pf_mu = list(filter(lambda x: GetPtEtaPhi(x,"PFch")[0], sr.pf_mu))
            sr.ID     = step

            sr.npf_ch += len(sr.pf_ch)
            sr.npf_ph += len(sr.pf_ph)
            sr.npf_ne += len(sr.pf_ne)
            sr.npf_mu += len(sr.pf_mu)
            evt.large_regions[lreg_ctr].npf_ch += len(sr.pf_ch)
            evt.large_regions[lreg_ctr].npf_ph += len(sr.pf_ph)
            evt.large_regions[lreg_ctr].npf_ne += len(sr.pf_ne)
            evt.large_regions[lreg_ctr].npf_mu += len(sr.pf_mu)
            evt.npf_ch += len(sr.pf_ch)
            evt.npf_ph += len(sr.pf_ph)
            evt.npf_ne += len(sr.pf_ne)
            evt.npf_mu += len(sr.pf_mu)

            # increment counters
            # (add support for large region ctr later)
            ctr = ctr+1
            sreg_ctr = ctr % EM_NSMALL
            if ctr and (ctr % EM_NSMALL)==0: 
                evt_ctr = evt_ctr+1

    return evts

def GetSimulationDataPFPUPPI(evts, parser):
    '''
    Simulation data is stored across multiple files,
    each corresponding to a time series of single-object outputs.
    #files is #em+#calo+#track, which may equal 15+20+25=60 --> 22+13+15+2=52?
    Each line corresonds to a 'step' of the output and...
    TODO add comment on length of each file.

    LOL
    1
    2
    2
    3
    4
    4
    5
    6
    6
    ...

    '''

    # local vars
    SIM_NEVENTS = NEVENTS if (NEVENTS>0) else 6
    SIM_NLARGE  = NLARGE if (NLARGE>0) else 1
    SIM_NSMALL  = NSMALL if (NSMALL>0) else 18
    #SIM_NSMALL  = 18
    # SIM_NSMALL  = 27

    sreg_ctr = 0
    lreg_ctr = 0
    evt_ctr  = 0

    # mask = [[0]*SIM_NSMALL for i in range(SIM_NEVENTS)]
    # summask = [0]*SIM_NSMALL
    # for i in range(SIM_NEVENTS): print("".join(str(mask[i])))
    # print('-')
    # mask[3][7]=1

    for ii in range(NTRACK+NPHOTON+NSELCALO+NMU):
        with open("{}/sim_output_fiber_{}.dat".format(parser.sim_output_dir,ii),'r') as f:
            ctr=0
            line_ctr=-1 # so syncs to 0 after initial increment
            for l in f:
                #print('ctr = ',ctr)
                val = int(l,base=16)
                #ignore repeated lines (18->27 issue)
                line_ctr+=1 #_always_ increment
                #if line_ctr%3==2: continue
                # ignore empty lines
                sreg_ctr = ctr % SIM_NSMALL
                evt_ctr  = int(ctr / SIM_NSMALL)
                if evt_ctr >= NEVENTS:
                    break
                if evt_ctr >= len(evts): 
                    print("Trying to fill a simulation event that doesn't exist! evt_ctr={}".format(evt_ctr))
                    return
                # (assume only one large region for now)
                if sreg_ctr >= len(evts[evt_ctr].large_regions[lreg_ctr].subregions):
                    evts[evt_ctr].large_regions[lreg_ctr].resize_subregions(sreg_ctr+1)
                    # todo todo
                    # print("adding subregion {} in PF/PUPPI step (event {})".format(sreg_ctr+1,evt_ctr))

                sr = evts[evt_ctr].large_regions[lreg_ctr].subregions[sreg_ctr]
                if GetPFChParams(val)[1]: #pt is stored in same bits for charged and neutral PF objs
                    # print("Event {} SR {} file {}".format(evt_ctr,sreg_ctr,ii))
                    # print("  Charged pt={1}, eta={3}, phi={4}, puppi={0}, ID={2}".format(*GetPFChParams(val)) )
                    # print("  Neutral pt={1}, eta={3}, phi={4}, z0={2}, ID={0}".format(*GetPFNeParams(val)) )
                    # mask[evt_ctr][sreg_ctr]+=1
                    if ii<NTRACK:
                        #if GetPtEtaPhi(val,"em")[1]==32: continue # TODO TEMP TEST
                        sr.pf_ch.append(val)
                        sr.npf_ch  += 1
                        evts[evt_ctr].npf_ch  += 1
                        evts[evt_ctr].large_regions[lreg_ctr].npf_ch += 1
                    elif ii<NTRACK+NPHOTON: 
                        #if GetPtEtaPhi(val,"calo")[1]==32: continue # TODO TEMP TEST
                        #if GetPtEtaPhi(val,"calo")[1]!=32: 
                        sr.pf_ph.append(val)
                        sr.npf_ph  += 1
                        evts[evt_ctr].npf_ph  += 1
                        evts[evt_ctr].large_regions[lreg_ctr].npf_ph += 1
                    elif ii<NTRACK+NPHOTON+NSELCALO:
                        #if GetPtEtaPhi(val,"tk")[1]==32: continue # TODO TEMP TEST
                        sr.pf_ne.append(val)
                        sr.npf_ne  += 1
                        evts[evt_ctr].npf_ne  += 1
                        evts[evt_ctr].large_regions[lreg_ctr].npf_ne += 1
                    elif ii<NTRACK+NPHOTON+NSELCALO+NMU:
                        #if GetPtEtaPhi(val,"mu")[1]==32: continue # TODO TEMP TEST
                        sr.pf_mu.append(val)
                        sr.npf_mu  += 1
                        evts[evt_ctr].npf_mu  += 1
                        evts[evt_ctr].large_regions[lreg_ctr].npf_mu += 1
                ctr += 1

    # #for i in range(SIM_NEVENTS): print("".join(str(mask[i])))
    # for i in range(SIM_NEVENTS): print((27*"{:2} ").format(*mask[i]))
    # mask2 = [[mask[i][j] if j>=18 else (mask[i][j]+mask[i][j-9] if j>=9 else 0) for j in range(SIM_NSMALL)] for i in range(SIM_NEVENTS)]
    # print('-')
    # for i in range(SIM_NEVENTS): print((27*"{:2} ").format(*mask2[i]))
    # print('-')
    # for j in range(SIM_NSMALL):
    #     for i in range(SIM_NEVENTS):
    #         summask[j] += mask[i][j]
    # print((27*"{:2} ").format(*summask))
    # for i in range(SIM_NEVENTS): print("".join(mask[i]))
    #return evts

def GetInputLink(xin, parser):
    #could definitely make this function fancier
    x = "{:0>16x}".format(xin)
    x = x[0:8] # often spread across two link check half
    infiles="{}/sim_input_fiber_*.dat".format(parser.sim_output_dir)
    cmd = "grep -i -n {} {}".format(x,infiles)
    output = subprocess.check_output([cmd],shell=True)
    if len(output) > 10000: return -1,-1 # e.g. x=0
    output = str(output)
    ind = output.find("sim_input_fiber")
    link = int((output[ind+len("sim_input_fiber_"):]).split(".")[0])
    clk = int(output.split(":")[1]) - 1 #0 index
    return link,clk



def GetInputs(parser, dump=True):
    NLINKS_TRACK=10
    NLINKS_EM=10
    NLINKS_CALO=10
    NLINKS_MUON=2
    NLINKS_REG=NLINKS_TRACK+NLINKS_EM+NLINKS_CALO+NLINKS_MUON
    BEG_TRACK=0
    BEG_EM=BEG_TRACK+NLINKS_TRACK
    BEG_CALO=BEG_EM+NLINKS_EM
    BEG_MUON=BEG_CALO+NLINKS_CALO
    # NTRACK=15
    # NEM=15
    # NCALO=15

    tracks=[]
    ems=[]
    calos=[]
    """
    Each line of the file has 96 entries, 1 per link
    There are a total of 378 lines in the file

    Link structure is 10+10+10+2=32 with dedicated links for each object
    (track+EM+calo+mu)
    This is duplicated three times to get 96. 
    (There is nothing special about this factor of 3, just fits nicely.)

    Each link is 'saturated' before moving to write to the next one.
    
    Links are rotated through 
    

    Numerology here is more complicated, but this is just the depth of the buffer 
    needed to accomplish the time multiplexed division into 96 links.
    Many (6) events each gets split into 18 TMUX_OUT slices.
    Data in arrives with TMUX_IN 36.    
    ... I still don't fully understand this
    """
    with open(parser.input_file,'r') as f:
        ctr=0
        for l in f:
            arr = list(map(lambda x : int(x,16), l.split()))
            # print (l)
            # print(len(arr))
            # print(arr)
            # exit(0) # 32*3+1
            step = arr[0]
            for off in [NLINKS_REG*i for i in range(3)]:
                # print( arr[1+BEG_TRACK+off:BEG_TRACK+NLINKS_TRACK+1+off] )
                # exit(0)
                tracks += arr[1+BEG_TRACK+off:BEG_TRACK+NLINKS_TRACK+1+off]
                ems    += arr[1+BEG_EM+off:BEG_EM+NLINKS_EM+1+off]
                calos  += arr[1+BEG_CALO+off:BEG_CALO+NLINKS_CALO+1+off]
                # these are 32-object arrays, holding all info


    return tracks, ems, calos

# print in red
def Warn(x): print('\033[91m'+x+'\033[0m')

def ComparePerEvent(em_events,sim_events,nevts=1):
    '''
    Loop over each event / tmuxed region and compare the outputs for each
    '''
    print("\n"+"="*80)
    print("""Test: Compare output tracks for all events
    The comparison is shown in 'steps', out of #EVENTS*(TMUX_OUT=18) total.
    Format is 'pt eta phi bits'. """)
    print("="*80)

    
    for ievt in range(nevts):
        print("\nEvent %u"%ievt)
        em_evt = em_events[ievt]
        sim_evt = sim_events[ievt]
        nsub = len(em_evt.large_regions[0].subregions)
        if nsub != len(sim_evt.large_regions[0].subregions):
            print("nsubreg mismatch!")
        for isr in range(nsub):
            em_sr = em_evt.large_regions[0].subregions[isr]
            sim_sr = sim_evt.large_regions[0].subregions[isr]
            # print(">em",isr,em_sr.tracks)
            # print("sim",isr,sim_sr.tracks)
            print(">em",isr,[GetTrackParams(t) for t in em_sr.tracks])
            print("sim",isr,[GetTrackParams(t) for t in sim_sr.tracks])
            print(">em",isr,["{:0>16x}".format(t) for t in em_sr.tracks])
            print("sim",isr,["{:0>16x}".format(t) for t in sim_sr.tracks])
    return
    for ii in range(CLKMAX):
        print("Step {}".format(ii))
        # remove zero entries from the track lists
        em_tracks[ii] = [x for x in em_tracks[ii] if x]
        sim_tracks[ii] = [x for x in sim_tracks[ii] if x]

        # Compare track sets
        common_tks = set(em_tracks[ii]).intersection(set(sim_tracks[ii]))
        em_only = set(em_tracks[ii]).difference(common_tks)
        sim_only = set(sim_tracks[ii]).difference(common_tks)
        if(len(set(em_tracks[ii])) != len(em_tracks[ii])): Warn("Warning! duplicate emulation track in this event!!")
        if(len(set(sim_tracks[ii])) != len(sim_tracks[ii])): Warn("Warning! duplicate simulation track in this event!!")
        print("  {} common tracks, {} only emulation, {} only simulation".format(len(common_tks),len(em_only),len(sim_only)))

        # dump the information for each of the tracks in this event
        for tk in em_tracks[ii]:
            pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
            if tk: 
                if tk in em_only: print("    EM only  tk {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
                else: print("    EM       tk {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
        for tk in sim_tracks[ii]:
            pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
            if tk:
                if tk in sim_only: print("    SIM only tk {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
                else: print("    SIM      tk {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
    print("="*80+"\n")
    return


def CheckAllObjects(em_events,sim_events):
    '''
    Compare the links that the common tracks arrive on
    '''
    all_em_tks=set()
    all_sim_tks=set()
    nevts = max(len(em_events),len(sim_events))
    if len(em_events) != len(sim_events):
        print("mismatch event length")
        #return
    else: print("checking",nevts,"events")
    #for ievt in range(1):
    for ievt in range(nevts):
        em_evt = em_events[ievt]
        sim_evt = sim_events[ievt]
        nsub = len(em_evt.large_regions[0].subregions)
        if nsub != len(sim_evt.large_regions[0].subregions):
            print("nsubreg mismatch!")
        for isr in range(nsub):
            em_sr = em_evt.large_regions[0].subregions[isr]
            sim_sr = sim_evt.large_regions[0].subregions[isr]
            for t in em_sr.tracks:
                all_em_tks.add(t)
            for t in sim_sr.tracks:
                all_sim_tks.add(t)
            
            # print(">em",isr,em_sr.tracks)
            # print("sim",isr,sim_sr.tracks)
            # print(">em",isr,[GetTrackParams(t) for t in em_sr.tracks])
            # print("sim",isr,[GetTrackParams(t) for t in sim_sr.tracks])
            # print(">em",isr,["{:0>16x}".format(t) for t in em_sr.tracks])
            # print("sim",isr,["{:0>16x}".format(t) for t in sim_sr.tracks])

    common_tks = all_em_tks.intersection(all_sim_tks)
    em_only = all_em_tks.difference(common_tks)
    sim_only = all_sim_tks.difference(common_tks)
    print()
    print("common_tks",["{:0>16X}".format(t) for t in common_tks])
    print()
    print("em_only   ",["{:0>16X}".format(t) for t in em_only   ])
    print()
    print("sim_only  ",["{:0>16X}".format(t) for t in sim_only  ])
    print()
    #print("em_only   ",[GetTrackParams(t) for t in em_only   ])
    print("sim_only   ",[GetTrackParams(t) for t in sim_only   ])

    return

    return
    all_em=[]
    all_sim=[]
    for ii in range(CLKMAX):
        for tk in [x for x in em_tracks[ii] if x]: all_em.append(tk)
        for tk in [x for x in sim_tracks[ii] if x]: all_sim.append(tk)
    common_tks = set(all_em).intersection(set(all_sim))
    
    em_links={}
    sim_links={}
    for tk in common_tks:
        em_links[tk]=[]
        sim_links[tk]=[]

    for ii in range(CLKMAX):
        for tk in [x for x in em_tracks[ii] if x and x in common_tks]:
            em_links[tk].append(ii % TMUX_OUT)
        for tk in [x for x in sim_tracks[ii] if x and x in common_tks]:
            sim_links[tk].append(ii % TMUX_OUT)

    # print the track characteristics
    print("\n"+"="*80)
    print("""Test: Compare output links over which common tracks are sent (disagreeing only!)
    Format is 'pt eta phi bits'. """)
    print("="*80)
    for tk in common_tks:
        if tuple(em_links[tk]) != tuple(sim_links[tk]):
            pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
            print("For track {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
            print("  Emulator  links: "+", ".join(map(str,em_links[tk])) )
            print("  Simulator links: "+", ".join(map(str,sim_links[tk])) )
    print("="*80+"\n")

    return

def TrackEmInSim(em_events,sim_events):
    '''
    Compare the links that the common tracks arrive on
    '''
    all_em_tks=set()
    sim_tks={} # maps tracks to event
    nevts = max(len(em_events),len(sim_events))
    if len(em_events) != len(sim_events):
        print("mismatch event length")
        #return
    else: print("checking",nevts,"events")
    #for ievt in range(1):

    import ROOT
    from collections import OrderedDict
    hists=OrderedDict()
    for i in range(6):
        n = "sim_clk_em_evt"+str(i)
        hists[n] = ROOT.TH1F(n,"",400,0,400)

    for evt in sim_events:
        for sr in evt.large_regions[0].subregions:
            for t in sr.tracks:
                # fill map with clock, maybe other things as well (ckl->tuple)
                clk = sr.track_clks[t]
                sim_tks[t] = clk

    for evt in em_events:
        for sr in evt.large_regions[0].subregions:
            for t in sr.tracks:
                if t in sim_tks:
                    #check where the tracks go...
                    l = sim_tks[t]
                    n = "sim_clk_em_evt"+str(i)
                    hists[n].Fill(float(l))
                else:
                    #print("emulation is missing a track")
                    pass
            
    f = ROOT.TFile("out.root","recreate")
    for i in range(len(em_events)):
        n = "sim_clk_em_evt"+str(i)
        hists[n].Write()
    f.Close()
                    
    return
    for ievt in range(nevts):
        em_evt = em_events[ievt]
        sim_evt = sim_events[ievt]
        nsub = len(em_evt.large_regions[0].subregions)
        if nsub != len(sim_evt.large_regions[0].subregions):
            print("nsubreg mismatch!")
        for isr in range(nsub):
            em_sr = em_evt.large_regions[0].subregions[isr]
            sim_sr = sim_evt.large_regions[0].subregions[isr]
            for t in em_sr.tracks:
                all_em_tks.add(t)
            for t in sim_sr.tracks:
                all_sim_tks.add(t)
            
            # print(">em",isr,em_sr.tracks)
            # print("sim",isr,sim_sr.tracks)
            # print(">em",isr,[GetTrackParams(t) for t in em_sr.tracks])
            # print("sim",isr,[GetTrackParams(t) for t in sim_sr.tracks])
            # print(">em",isr,["{:0>16x}".format(t) for t in em_sr.tracks])
            # print("sim",isr,["{:0>16x}".format(t) for t in sim_sr.tracks])

    common_tks = all_em_tks.intersection(all_sim_tks)
    em_only = all_em_tks.difference(common_tks)
    sim_only = all_sim_tks.difference(common_tks)
    print()
    print("common_tks",["{:0>16X}".format(t) for t in common_tks])
    print()
    print("em_only   ",["{:0>16X}".format(t) for t in em_only   ])
    print()
    print("sim_only  ",["{:0>16X}".format(t) for t in sim_only  ])
    print()
    #print("em_only   ",[GetTrackParams(t) for t in em_only   ])
    print("sim_only   ",[GetTrackParams(t) for t in sim_only   ])

def DumpHistograms(em_tracks,sim_tracks, outfile="out.root"):
    '''
    Loop over each event / tmuxed region and dump a number of outputs into a histogram
    '''
    import ROOT
    f = ROOT.TFile(outfile,"RECREATE")
    hPhi_EM = ROOT.TH1D("hPhi_EM","",40,0,2**10)
    hEta_EM = ROOT.TH1D("hEta_EM","",40,0,2**10)
    hPhi_SIM = ROOT.TH1D("hPhi_SIM","",40,0,2**10)
    hEta_SIM = ROOT.TH1D("hEta_SIM","",40,0,2**10)
    hPhi_COMMON = ROOT.TH1D("hPhi_COMMON","",40,0,2**10)
    hEta_COMMON = ROOT.TH1D("hEta_COMMON","",40,0,2**10)

    # add all non-zero entries to the total track lists
    all_em=[]
    all_sim=[]
    for ii in range(CLKMAX):
        for tk in [x for x in em_tracks[ii] if x]: all_em.append(tk)
        for tk in [x for x in sim_tracks[ii] if x]: all_sim.append(tk)

    # categorize
    common_tks = set(all_em).intersection(set(all_sim))
    em_only = set(all_em).difference(common_tks)
    sim_only = set(all_sim).difference(common_tks)
    
    for tk in em_only:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        hPhi_EMonly.Fill(phi)
        hEta_EMonly.Fill(eta)
    for tk in sim_only:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        hPhi_SIMonly.Fill(phi)
        hEta_SIMonly.Fill(eta)
    for tk in common_tks:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        hPhi_COMMON.Fill(phi)
        hEta_COMMON.Fill(eta)

    hPhi_EM.Write()
    hEta_EM.Write()
    hPhi_SIM.Write()
    hEta_SIM.Write()
    hPhi_COMMON.Write()
    hEta_COMMON.Write()
    f.Close()                
    return

def CompareOverlap(em_evts,sim_evts):
    '''
    Loop over each event / tmuxed region and dump a number of outputs into a histogram
    '''

    # add all non-zero entries to the total track lists
    all_em=[]
    all_sim=[]
    # for ii in range(CLKMAX):
    #     for tk in [x for x in em_tracks[ii] if x]: all_em.append(tk)
    #     for tk in [x for x in sim_tracks[ii] if x]: all_sim.append(tk)

    for e in em_evts:
        for l in e.large_regions:
            for s in l.subregions:
                all_em += s.tracks
    for e in sim_evts:
        for l in e.large_regions:
            for s in l.subregions:
                all_sim += s.tracks

    # categorize
    common_tks = set(all_em).intersection(set(all_sim))
    em_only = set(all_em).difference(common_tks)
    sim_only = set(all_sim).difference(common_tks)
    
    print("\n"+"="*80)
    print("""Test: Compare output tracks for all events
    Format is 'pt eta phi bits'. """)
    print("="*80)
    print("Emulation-only tracks")
    for tk in em_only:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        print("  tk {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
    if not len(em_only): print("  (...none found...)")

    print("Simulation-only tracks")
    for tk in sim_only:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        print("  tk {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
    if not len(sim_only): print("  (...none found...)")

    print("Commonly found tracks")
    for tk in common_tks:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        print("  tk {:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,tk))
    if not len(common_tks): print("  (...none found...)")
    print("="*80+"\n")

    return

def DumpInputs(in_tracks, in_ems, in_calos):
                        
    # Process inputs
    # ii=0
    # while ii < len(in_tracks):
    #     print (in_tracks[ii:ii+3])
    #     ii += 3
    # exit(0)

    #first get rid of all of the zeros
    in_tracks = [x for x in in_tracks if x]
    in_ems    = [x for x in in_ems    if x]
    in_calos  = [x for x in in_calos  if x]
    for tk in in_tracks:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        print("Input tracks", pt, eta, phi, "{:0>16x}".format(tk))

    exit(0)
    # do we find EM only in the inputs?
    all_in = set(in_tracks)
    em_no_input = em_only.difference(all_in)
    print("We find a total of {} input tracks ({} unique),".format(len(in_tracks),len(all_in)))
    print(" and we find that {} of the {} emulated-but-not-simulated tracks".format(len(em_no_input),len(em_only)))
    print(" do not match any input track!!")

def DumpEventsToText(evts, fname, fine=False):
    with open(fname,'w') as f:
        f.write("Dumping events to file: {}\n".format(fname))

        f.write("Found {} events \n".format(len(evts)))
        ei=0
        for e in evts:
            lrs = e.large_regions
            f.write("Event {} has {} large regions, ".format(ei,len(lrs)))
            f.write("{} tracks, ems, calos, mus \n".format(e.counts()))
            li=0
            for lr in lrs:
                srs = lr.subregions
                f.write("  Large region {} has {} small regions, ".format(li,len(srs)))
                f.write("{} tracks, ems, calos, mus \n".format(lr.counts()))
                si=0
                for sr in srs:
                    f.write("    Small region {} has ".format(si))
                    f.write("{} tracks, ems, calos, mus \n".format(sr.counts()))
                    if fine:
                        f.write("      Tracks ({}):\n".format(len(sr.tracks)))
                        for x in sr.tracks:
                            pt, pte, eta, phi, z0, qual = GetTrackParams(x)
                            f.write("         {:5} {:5} {:5} {:0>16x}\n".format(pt,eta,phi,x))
                        f.write("      EMs ({}):\n".format(len(sr.ems)))
                        for x in sr.ems:
                            pt, pte, eta, phi = GetEMParams(x)
                            f.write("         {:5} {:5} {:5} {:0>16x}\n".format(pt,eta,phi,x))
                        f.write("      Calos ({}):\n".format(len(sr.calos)))
                        for x in sr.calos:
                            pt, empt, eta, phi, isEM = GetCaloParams(x)
                            f.write("         {:5} {:5} {:5} {:0>16x}\n".format(pt,eta,phi,x))
                        f.write("      Muons ({}):\n".format(len(sr.mus)))
                        for x in sr.mus:
                            pt, pte, eta, phi, z0, qual = GetTrackParams(x)
                            f.write("         {:5} {:5} {:5} {:0>16x}\n".format(pt,eta,phi,x))


                    si += 1
                li += 1
            ei += 1

                    # print("  ntracks",s.ntracks)
                    # print("  nems",s.nems)
                    # print("  ncalos",s.ncalos)
                    # print("  nmus",s.nmus)

        f.write("\n")
    print("Wrote "+fname)

def DumpEventsToTextFine(evts, fname):
    DumpEventsToText(evts,fname,True)

def GetCommonEmSim(emlist, simlist):
    em = set(emlist)
    sim = set(simlist)
    common= em.intersection(sim)
    return common, em.difference(common), sim.difference(common)

def DumpCollection(parser, xlist, title, tag, mult=True, inLink=True,
                   printReg=False, printRegPF=False, sim_evt=None, em_evt=None):
    record=title
    suff=""
    if mult: record += " ({})".format(len(xlist))
    record += ":\n"
    nblank = 0
    match = True
    while title[nblank]==" ": nblank+=1
    for x in xlist:
        suff=""
        pt, eta, phi = GetPtEtaPhi(x, tag)
        record += " "*nblank + "{:5} {:5} {:5} {:0>16x}".format(pt,eta,phi,x)
        if sim_evt and em_evt:
            if printReg:
                sim_srs = tuple(sorted([ xx[1] for xx in sim_evt.findRegions(x)]))
                em_srs =  tuple(sorted([ xx[1] for xx in em_evt.findRegions(x) ]))
                sim_srs2 = [ORDERING[x] for x in sim_srs]
                em_srs2 =  [ORDERING[x] for x in  em_srs]
                match = (sim_srs==em_srs)
                suff = " match?={} (SRs sim vs emulator are  {} vs {} OR IN ETA/PHI {} vs {})".format(int(match),sim_srs,em_srs,sim_srs2,em_srs2)
            if printRegPF:
                sim_srs = tuple(sorted([ xx[1] for xx in sim_evt.findRegionsPF(x)]))
                em_srs =  tuple(sorted([ xx[1] for xx in  em_evt.findRegionsPF(x) ]))
                sim_srs2 = [ORDERING[x] for x in sim_srs]
                em_srs2 =  [ORDERING[x] for x in  em_srs]
                match = (sim_srs==em_srs)
                suff = " match?={} (SRs sim vs emulator are  {} vs {} OR IN ETA/PHI {} vs {})".format(int(match),sim_srs,em_srs,sim_srs2,em_srs2)
        if inLink and ((("ommon" in title) and not match) or ("mulation" in title)): 
        # adding this criterion so that the grepping is only run if needed
            suff = " (link {:2}, clock {:4})".format(*GetInputLink(x, parser)) + suff
        record += suff+"\n"
    return record

def DumpCollectionReg(parser, xlist, title, tag, sim_evt, em_evt, mult=True, inLink=True):
    return DumpCollection(parser, xlist, title, tag, mult=mult, inLink=inLink,
                          printReg=True, sim_evt=sim_evt, em_evt=em_evt)

def DumpEventComparison(parser, em_evts, sim_evts, fname, 
                        nevts=1, compareSmallRegion=False, compareEvent=True, tracksOnly=True, summary=True):
    
    ts = time.time()
    with open(fname,'w') as f:
        f.write("Dumping comparison to file: {}\n".format(fname))
        f.write("    (object outputs are all: pt  eta  phi  64bID)\n")

        f.write("Found {} events in sim, {} in emulation \n".format(len(sim_evts),len(em_evts)))
        if len(sim_evts) != len(em_evts): 
            f.write("Event mismatch!\n")
            return

        for ei in range(len(sim_evts)):
            sim_e = sim_evts[ei]
            em_e = em_evts[ei]
            sim_lrs = sim_e.large_regions
            em_lrs = em_e.large_regions
            if len(sim_lrs) != len(em_lrs): 
                f.write("Large region mismatch!\n")
                return
            f.write("Event {} has (#track,em,calo,mu) = EM {} vs SIM {}\n".format(ei,em_evts[ei].counts(),sim_evts[ei].counts()))

            if summary:
                com_tracks, em_tracks, sim_tracks = GetCommonEmSim(em_e.allTracks(), sim_e.allTracks())
                f.write("     {} Common, {} Emulation-only, {} Simulation-only Tracks".format(
                    len(com_tracks),len(em_tracks),len(sim_tracks)))
                com_ems, em_ems, sim_ems = GetCommonEmSim(em_e.allEMs(), sim_e.allEMs())
                f.write("     {} Common, {} Emulation-only, {} Simulation-only EM Calos".format(
                    len(com_ems),len(em_ems),len(sim_ems)))
                com_calos, em_calos, sim_calos = GetCommonEmSim(em_e.allCalos(), sim_e.allCalos())
                f.write("     {} Common, {} Emulation-only, {} Simulation-only Calos".format(
                    len(com_calos),len(em_calos),len(sim_calos)))
                com_muons, em_muons, sim_muons = GetCommonEmSim(em_e.allMuons(), sim_e.allMuons())
                f.write("     {} Common, {} Emulation-only, {} Simulation-only Muons".format(
                    len(com_muons),len(em_muons),len(sim_muons)))
                continue
        
            if compareEvent:
                tag="tk"
                com_tracks, em_tracks, sim_tracks = GetCommonEmSim(em_e.allTracks(), sim_e.allTracks())
                f.write( DumpCollectionReg(parser, com_tracks,"      Common Tracks",tag,sim_e, em_e) )
                f.write( DumpCollection(parser,  em_tracks,"      Emulation-only Tracks",tag) )
                f.write( DumpCollection(parser, sim_tracks,"      Simulation-only Tracks",tag) )
                #findRegions()

                if not tracksOnly:
                    tag="em"
                    com_ems, em_ems, sim_ems = GetCommonEmSim(em_e.allEMs(), sim_e.allEMs())
                    f.write( DumpCollection(parser, com_ems,"      Common EMs",tag) )
                    f.write( DumpCollection(parser,  em_ems,"      Emulation-only EMs",tag) )
                    f.write( DumpCollection(parser, sim_ems,"      Simulation-only EMs",tag) )
                    
                    tag="calo"
                    com_calos, em_calos, sim_calos = GetCommonEmSim(em_e.allCalos(), sim_e.allCalos())
                    f.write( DumpCollection(parser, com_calos,"      Common Calos",tag) )
                    f.write( DumpCollection(parser,  em_calos,"      Emulation-only Calos",tag) )
                    f.write( DumpCollection(parser, sim_calos,"      Simulation-only Calos",tag) )
                    
                    tag="mu"
                    com_mus, em_mus, sim_mus = GetCommonEmSim(em_e.allMuons(), sim_e.allMuons())
                    f.write( DumpCollection(parser, com_mus,"      Common Muons",tag) )
                    f.write( DumpCollection(parser,  em_mus,"      Emulation-only Muons",tag) )
                    f.write( DumpCollection(parser, sim_mus,"      Simulation-only Muons",tag) )


                   
            for li in range(len(sim_lrs)):
                sim_srs = sim_lrs[li].subregions
                em_srs  = em_lrs[li].subregions
                if len(sim_srs) != len(em_srs): 
                    f.write("Small region mismatch!\n")
                    return
                f.write("  Large region {} has (#track,em,calo,mu) = EM {} vs SIM {}\n".format(
                    li,em_lrs[li].counts(),sim_lrs[li].counts()))

                for si in range(len(sim_srs)):
                    sim_sr = sim_srs[si]
                    em_sr  = em_srs[si]
                    f.write("    Small region {} has (#track,em,calo,mu) = EM {} vs SIM {}, (matchSR={})\n".format(
                        si,em_sr.counts(),sim_sr.counts(),em_sr.counts()==sim_sr.counts()))
                    if not (em_sr.counts()==sim_sr.counts()):
                        f.write( DumpCollection(parser,  em_sr.calos, "DEBUG  EM ","calo") )
                        f.write( DumpCollection(parser, sim_sr.calos, "DEBUG SIM ","calo") )

                    if compareSmallRegion:
                        tag="tk"
                        com_tracks, em_tracks, sim_tracks = GetCommonEmSim(em_sr.tracks, sim_sr.tracks)
                        f.write( DumpCollection(parser, com_tracks,"      Common Tracks",tag) )
                        f.write( DumpCollection(parser,  em_tracks,"      Emulation-only Tracks",tag) )
                        f.write( DumpCollection(parser, sim_tracks,"      Simulation-only Tracks",tag) )

                        if not tracksOnly:
                            tag="em"
                            com_ems, em_ems, sim_ems = GetCommonEmSim(em_sr.ems, sim_sr.ems)
                            f.write( DumpCollection(parser, com_ems,"      Common EMs",tag) )
                            f.write( DumpCollection(parser,  em_ems,"      Emulation-only EMs",tag) )
                            f.write( DumpCollection(parser, sim_ems,"      Simulation-only EMs",tag) )
                            
                            tag="calo"
                            com_calos, em_calos, sim_calos = GetCommonEmSim(em_sr.calos, sim_sr.calos)
                            f.write( DumpCollection(parser, com_calos,"      Common Calos",tag) )
                            f.write( DumpCollection(parser,  em_calos,"      Emulation-only Calos",tag) )
                            f.write( DumpCollection(parser, sim_calos,"      Simulation-only Calos",tag) )
                            
                            tag="mu"
                            com_mus, em_mus, sim_mus = GetCommonEmSim(em_sr.mus, sim_sr.mus)
                            f.write( DumpCollection(parser, com_mus,"      Common Muons",tag) )
                            f.write( DumpCollection(parser,  em_mus,"      Emulation-only Muons",tag) )
                            f.write( DumpCollection(parser, sim_mus,"      Simulation-only Muons",tag) )


        f.write("\n")
    te = time.time()
    print("Wrote {} (in {:.1f} sec)".format(fname, te-ts) )


def DumpEventComparisonPFPUPPI(parser, em_evts, sim_evts, fname, 
                               nevts=1, compareSmallRegion=False, compareEvent=True, chargedOnly=True, summary=True, skipCommon=False):
    
    with open(fname,'w') as f:
        f.write("Dumping comparison to file: {}\n".format(fname))
        f.write("    (object outputs are all: pt eta  phi  64bID)\n")

        f.write("Found {} events in sim, {} in emulation \n".format(len(sim_evts),len(em_evts)))
        if len(sim_evts) != len(em_evts): 
            f.write("Event mismatch!\n")
            return

        for ei in range(len(sim_evts)):
            sim_e = sim_evts[ei]
            em_e = em_evts[ei]
            sim_lrs = sim_e.large_regions
            em_lrs = em_e.large_regions
            if len(sim_lrs) != len(em_lrs): 
                f.write("Large region mismatch!\n")
                return
            f.write("Event {} has (#ch,ph,neu,mu) = EM {} vs SIM {}\n".format(
                ei,em_evts[ei].countsPFPUPPI(),sim_evts[ei].countsPFPUPPI()))
            # f.write( DumpCollection(parser,  em_e.allPF("ch"), "DEBUG  EM ","PFch") )
            # f.write( DumpCollection(parser, sim_e.allPF("ch"), "DEBUG SIM ","PFch") )

            for tag in ["ch","ph","ne","mu"]:
                tagpf="PFne"
                if tag=="ch" or tag=="mu": tagpf="PFch" 
                #evt.findRegions(x)
                common, em_only, sim_only = GetCommonEmSim(em_e.allPF(tag), sim_e.allPF(tag))
                if not skipCommon: f.write( DumpCollection(parser,   common,"      Common "+tag,         tagpf,
                                        inLink = False,printRegPF=True,sim_evt=sim_e,em_evt=em_e) )
                f.write( DumpCollection(parser,  em_only,"      Emulation-only "+tag, tagpf,inLink = False) )
                f.write( DumpCollection(parser, sim_only,"      Simulation-only "+tag,tagpf,inLink = False) )
                    
            continue
            for li in range(len(sim_lrs)):
                sim_srs = sim_lrs[li].subregions
                em_srs  = em_lrs[li].subregions
                # if len(sim_srs) != len(em_srs): 
                #     f.write("Small region mismatch!\n")
                #     return
                f.write("  Large region {} has (#ch,ph,neu,mu) = EM {} vs SIM {}\n".format(
                    li,em_lrs[li].countsPFPUPPI(),sim_lrs[li].countsPFPUPPI()))


                for si in range(len(sim_srs)):
                    if si >= len(em_srs):
                        sim_sr = sim_srs[si]
                        f.write("    Small region {} has (#ch,ph,neu,mu) = SIM {} \n".format(si, sim_sr.countsPFPUPPI()))
                    else:
                        sim_sr = sim_srs[si]
                        em_sr  = em_srs[si]
                        f.write("    Small region {} has (#ch,ph,neu,mu) = EM {} vs SIM {}, (matchSR={})\n".format(
                            si,em_sr.countsPFPUPPI(),sim_sr.countsPFPUPPI(),em_sr.countsPFPUPPI()==sim_sr.countsPFPUPPI()))


            #     for si in range(len(sim_srs)):
            #         sim_sr = sim_srs[si]
            #         em_sr  = em_srs[si]
            #         f.write("    Small region {} has (#ch,ph,neu,mu) = EM {} vs SIM {}, (matchSR={})\n".format(
            #             si,em_sr.countsPFPUPPI(),sim_sr.countsPFPUPPI(),em_sr.countsPFPUPPI()==sim_sr.countsPFPUPPI()))
            #         if not (em_sr.countsPFPUPPI()==sim_sr.countsPFPUPPI()):
            #             f.write( DumpCollection(parser,  em_sr.calos, "DEBUG  EM ","calo") )
            #             f.write( DumpCollection(parser, sim_sr.calos, "DEBUG SIM ","calo") )

            #         if compareSmallRegion:
            #             tag="tk"
            #             com_tracks, em_tracks, sim_tracks = GetCommonEmSim(em_sr.tracks, sim_sr.tracks)
            #             f.write( DumpCollection(parser, com_tracks,"      Common Tracks",tag) )
            #             f.write( DumpCollection(parser,  em_tracks,"      Emulation-only Tracks",tag) )
            #             f.write( DumpCollection(parser, sim_tracks,"      Simulation-only Tracks",tag) )

            #             if not tracksOnly:
            #                 tag="em"
            #                 com_ems, em_ems, sim_ems = GetCommonEmSim(em_sr.ems, sim_sr.ems)
            #                 f.write( DumpCollection(parser, com_ems,"      Common EMs",tag) )
            #                 f.write( DumpCollection(parser,  em_ems,"      Emulation-only EMs",tag) )
            #                 f.write( DumpCollection(parser, sim_ems,"      Simulation-only EMs",tag) )
                            
            #                 tag="calo"
            #                 com_calos, em_calos, sim_calos = GetCommonEmSim(em_sr.calos, sim_sr.calos)
            #                 f.write( DumpCollection(parser, com_calos,"      Common Calos",tag) )
            #                 f.write( DumpCollection(parser,  em_calos,"      Emulation-only Calos",tag) )
            #                 f.write( DumpCollection(parser, sim_calos,"      Simulation-only Calos",tag) )
                            
            #                 tag="mu"
            #                 com_mus, em_mus, sim_mus = GetCommonEmSim(em_sr.mus, sim_sr.mus)
            #                 f.write( DumpCollection(parser, com_mus,"      Common Muons",tag) )
            #                 f.write( DumpCollection(parser,  em_mus,"      Emulation-only Muons",tag) )
            #                 f.write( DumpCollection(parser, sim_mus,"      Simulation-only Muons",tag) )


        f.write("\n")
    print("Wrote "+fname)


def HistogramEventComparisons(parser, em_evts, sim_evts, fname):
    
    import ROOT
    from collections import OrderedDict
    hists=OrderedDict()

    # helper func (no errors)
    def book(h,name,n,a,b,title=""):
        h[name]=ROOT.TH1F(name,title,n,a,b)
    def book2(h,name,nx,ax,bx,ny,ay,by,title=""):
        h[name]=ROOT.TH2F(name,title,nx,ax,bx,ny,ay,by)
    def bookp(h,name,nx,ax,bx,ay,by,title="",err="s"):
        h[name]=ROOT.TProfile(name,title,nx,ax,bx,ay,by,err)
    def IncrementBinContent(h, b, x):
        h.SetBinContent(b, h.GetBinContent(b) + x)
    
    #plots
    objs = ["track","em","calo","mu"]
    for obj in objs:
        for x in ["sim","em"]:
            # number of objects found versus event number
            book(hists,"{}_{}_perEvent".format(x,obj),NEVENTS,-0.5,NEVENTS-0.5,";Event #")
            # number of objects found versus small region
            book(hists,"{}_{}_perSmallRegion".format(x,obj),NSMALL,-0.5,NSMALL-0.5,";Region index")
            # number of objects found versus eta 
            book(hists,"{}_{}_eta".format(x,obj),300,-250-0.5,50-0.5,";eta bits")
            # number of objects found versus phi
            book(hists,"{}_{}_phi".format(x,obj),1040,-520-0.5,520-0.5,";phi bits")
        # number of objects found versus input link
        book(hists,"{}_perInputLink".format(obj),96,-0.5,96-0.5,";Input link index")
        # number of objects found versus clock
        book(hists,"{}_perInputClock".format(obj),400,-0.5,399.5,";Input clock") # 384 = 6 events
        # number of objects found versus clock
        book2(hists,"{}_perLink_perInputClock".format(obj),96,-0.5,96-0.5, 400,-0.5,399.5,";Link index;Input clock") # 384 = 6 events
    

    for ei in range(len(sim_evts)):
        sim_e = sim_evts[ei]
        em_e = em_evts[ei]
        sim_lrs = sim_e.large_regions
        em_lrs = em_e.large_regions
        if len(sim_lrs) != len(em_lrs): 
            print("Large region mismatch!\n")
            return

        # note: assumed 'counts()' is indexed same as 'objs'
        for oi, obj in enumerate(objs):
            hists["em_{}_perEvent".format(obj)].SetBinContent(ei+1, em_evts[ei].counts()[oi])
            hists["sim_{}_perEvent".format(obj)].SetBinContent(ei+1, sim_evts[ei].counts()[oi])
            
        # just one LR for now
        for li in range(len(sim_lrs)):
            sim_srs = sim_lrs[li].subregions
            em_srs  = em_lrs[li].subregions
            if len(sim_srs) != len(em_srs): 
                print("Small region mismatch!\n")
                return

            for si in range(len(sim_srs)):
                sim_sr = sim_srs[si]
                em_sr  = em_srs[si]

                for oi, obj in enumerate(objs):
                    # small region counts
                    IncrementBinContent(hists["em_{}_perSmallRegion".format(obj)],si+1, em_sr.counts()[oi])
                    IncrementBinContent(hists["sim_{}_perSmallRegion".format(obj)],si+1, sim_sr.counts()[oi])

                    for x in em_sr.tracks: hists["em_track_eta"].Fill(GetPtEtaPhi(x,"tk")[1])
                    for x in em_sr.calos: hists["em_calo_eta"].Fill(GetPtEtaPhi(x,"calo")[1])
                    for x in em_sr.ems: hists["em_em_eta"].Fill(GetPtEtaPhi(x,"em")[1])
                    for x in em_sr.mus: hists["em_mu_eta"].Fill(GetPtEtaPhi(x,"mu")[1])
                    for x in sim_sr.tracks: hists["sim_track_eta"].Fill(GetPtEtaPhi(x,"tk")[1])
                    for x in sim_sr.calos: hists["sim_calo_eta"].Fill(GetPtEtaPhi(x,"calo")[1])
                    for x in sim_sr.ems: hists["sim_em_eta"].Fill(GetPtEtaPhi(x,"em")[1])
                    for x in sim_sr.mus: hists["sim_mu_eta"].Fill(GetPtEtaPhi(x,"mu")[1])

                    for x in em_sr.tracks: hists["em_track_phi"].Fill(GetPtEtaPhi(x,"tk")[2])
                    for x in em_sr.calos: hists["em_calo_phi"].Fill(GetPtEtaPhi(x,"calo")[2])
                    for x in em_sr.ems: hists["em_em_phi"].Fill(GetPtEtaPhi(x,"em")[2])
                    for x in em_sr.mus: hists["em_mu_phi"].Fill(GetPtEtaPhi(x,"mu")[2])
                    for x in sim_sr.tracks: hists["sim_track_phi"].Fill(GetPtEtaPhi(x,"tk")[2])
                    for x in sim_sr.calos: hists["sim_calo_phi"].Fill(GetPtEtaPhi(x,"calo")[2])
                    for x in sim_sr.ems: hists["sim_em_phi"].Fill(GetPtEtaPhi(x,"em")[2])
                    for x in sim_sr.mus: hists["sim_mu_phi"].Fill(GetPtEtaPhi(x,"mu")[2])
                    
                    # link info for common objects 
                    common_collection=None
                    #try eval here :) ?
                    if obj=="track": common_collection = GetCommonEmSim(em_sr.tracks, sim_sr.tracks)[0]
                    if obj=="em"   : common_collection = GetCommonEmSim(em_sr.ems, sim_sr.ems)[0]
                    if obj=="calo" : common_collection = GetCommonEmSim(em_sr.calos, sim_sr.calos)[0]
                    if obj=="mu"   : common_collection = GetCommonEmSim(em_sr.mus, sim_sr.mus)[0]
                    
                    for link, clk in [GetInputLink(x,parser) for x in common_collection]:
                        hists["{}_perInputLink".format(obj)].Fill(link)
                        hists["{}_perInputClock".format(obj)].Fill(clk)
                        hists["{}_perLink_perInputClock".format(obj)].Fill(link,clk)


    f = ROOT.TFile(fname,"recreate")
    for n in hists:
        hists[n].Write()
    f.Close()
    print("Wrote file: "+fname)

    
def PlotEventComparisons(fin_name, fout_name,outdir):
    import ROOT
    from collections import OrderedDict
    hists=OrderedDict()
    
    #get hists from file
    f = ROOT.TFile(fin_name,"read")
    objs = ["track","em","calo","mu"]
    for obj in objs:
        for x in ["sim","em"]:
            for n in ["perEvent","perSmallRegion","eta","phi"]:
                s = "{}_{}_{}".format(x,obj,n)
                hists[s] = f.Get(s)
        for n in ["perInputLink","perInputClock","perLink_perInputClock"]:
            s = "{}_{}".format(obj,n)
            hists[s] = f.Get(s)

    from histtools import DrawHistsEmSim,DrawHistsLinkClk,DrawHistsLinkClk2D
    em_hists  = [hists[ "em_{}_perEvent".format(obj)] for obj in objs]
    sim_hists = [hists["sim_{}_perEvent".format(obj)] for obj in objs]
    DrawHistsEmSim(em_hists, sim_hists,outdir+"/perEvent")
    em_hists  = [hists[ "em_{}_perSmallRegion".format(obj)] for obj in objs]
    sim_hists = [hists["sim_{}_perSmallRegion".format(obj)] for obj in objs]
    DrawHistsEmSim(em_hists, sim_hists,outdir+"/perSmallRegion")

    em_hists  = [hists[ "em_{}_eta".format(obj)] for obj in objs]
    sim_hists = [hists["sim_{}_eta".format(obj)] for obj in objs]
    DrawHistsEmSim(em_hists, sim_hists,outdir+"/eta")
    em_hists  = [hists[ "em_{}_phi".format(obj)] for obj in objs]
    sim_hists = [hists["sim_{}_phi".format(obj)] for obj in objs]
    DrawHistsEmSim(em_hists, sim_hists,outdir+"/phi")



    DrawHistsLinkClk( [hists[o+"_perInputLink"] for o in objs], outdir+"/perInputLink")
    DrawHistsLinkClk( [hists[o+"_perInputClock"] for o in objs], outdir+"/perInputClock")

    DrawHistsLinkClk2D( [hists[o+"_perLink_perInputClock"] for o in objs], outdir+"/perLinkperInputClock")

    fout = ROOT.TFile(fout_name,"recreate")
    fout.Close()
    f.Close()
    print("Wrote file: "+fout_name)


if __name__ == "__main__":

    # GetTrackParams(0x40411333000B0012)
    # GetCaloParams(0x40411333000B0012)
    # GetEMParams(0x40411333000B0012)
    # exit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="",  dest = "input_file", help = "regionizer input from emulator")
    parser.add_argument("-e", "--emulation-output", default="outputs.txt",  dest = "emulator_output", help = "regionizer output from emulator")
    parser.add_argument("-l", "--emulation-output-pfpuppi", default="",  dest = "emulator_output_pfpuppi", help = "PF+PUPPI output text file from emulator")
    parser.add_argument("-s", "--simulation-data-dir", default="/home/therwig/sandbox/otsdaq-cms-firmware/regionizer_full_small/sim/sim_data/", dest = "sim_output_dir", help = "regionizer output directory from simulation")
    parser.add_argument("-a", "--all", action='store_true', default=False,  dest = "all_dumps", help = "produce object dumps for all object types")
    parser.add_argument("-r", "--root", action='store_true', default=False, help = "produce root histograms from logs")
    parser.add_argument("-p", "--plot", action='store_true', default=False, help = "produce root plots from hists")
    parser = parser.parse_args(sys.argv[1:])
    print()
    if len(parser.input_file): print("Checking inputs from: "+parser.input_file)
    print("Reading from emulation regionizer output:        "+parser.emulator_output)
    print("Reading from emulation PF+PUPPI (layer1) output: "+parser.emulator_output_pfpuppi)
    print("Reading from simulation output: "+parser.sim_output_dir)

    if len(parser.input_file):
        in_tracks, in_ems, in_calos = GetInputs(parser)
        print (in_tracks)

    # information is dumped here
    dname = "dumps"
    os.system("mkdir -p "+dname)

    if parser.plot:
        import sys 
        sys.argv.append( '-b' )
        PlotEventComparisons(dname+"/verification_hists.root",dname+"/verification_plots.root",dname)
        print("Finished plotting, exiting now...")
        exit(0)

    # fill objects from log files
    print("Getting emulation data")
    em_events  = GetEmulationData(parser)
    print("Getting simulation data")
    sim_events = GetSimulationData(parser)

    # get PF+PUPPI outputs
    if len(parser.emulator_output_pfpuppi):
        GetEmulationDataPFPUPPI(em_events, parser)
        GetSimulationDataPFPUPPI(sim_events, parser)
    
    # for e in sim_events:
    #     print( len(e.allPF("ch")) )
    # exit(0)
    # for ei,e in enumerate(em_events):
    #     #        print ( e.large_regions )
    #     for si,sr in enumerate(e.large_regions[0].subregions):
    #         print("Evt {}, SR {:2}: Found {:2} charged, {:2} photons, {:2} neutrals, and {} muons".format(ei,si,sr.npf_ch,sr.npf_ph,sr.npf_ne,sr.npf_mu))
    #         if True:
    #             for x in sr.pf_ch:
    #                 print("Charged pt={1}, eta={3}, phi={4}, puppi={0}, ID={2}".format(*GetPFChParams(x)) )
    #             for x in sr.pf_ph:
    #                 print("Photon  pt={1}, eta={3}, phi={4}, z0={2}, ID={0}".format(*GetPFNeParams(x)) )
    #             for x in sr.pf_ne:
    #                 print("Neutral pt={1}, eta={3}, phi={4}, z0={2}, ID={0}".format(*GetPFNeParams(x)) )
    #             for x in sr.pf_mu:
    #                 print("Muon    pt={1}, eta={3}, phi={4}, puppi={0}, ID={2}".format(*GetPFChParams(x)) )
    #             # for x in sr.pf_ch:
    #             #     print (x)
    #             #     print (GetPFChParams(x))
    #             #     print( "charged track with pt {}, eta {}, phi {}".format(*GetPtEtaPhi(x,"PFch")) )
    # exit(0)


    if len(parser.emulator_output_pfpuppi):
        DumpEventComparisonPFPUPPI(parser, em_events,sim_events,dname+"/comp_pfpuppi.txt",nevts=NEVENTS)
        DumpEventComparisonPFPUPPI(parser, em_events,sim_events,dname+"/comp_pfpuppi_short.txt",nevts=NEVENTS,skipCommon=True)
    #exit(0)

    #basic event information w/ and w/o object details
    DumpEventsToText(em_events,dname+"/events_emu.txt")
    DumpEventsToText(sim_events,dname+"/events_sim.txt")
    DumpEventsToTextFine(em_events,dname+"/events_emu_fine.txt")
    DumpEventsToTextFine(sim_events,dname+"/events_sim_fine.txt")
    
    #check which object overlap, by event and region
    DumpEventComparison(parser, em_events,sim_events,dname+"/summary_comp_event_tracks.txt",nevts=NEVENTS,
                        compareSmallRegion=False, compareEvent=True, tracksOnly=True, summary=True)
    DumpEventComparison(parser, em_events,sim_events,dname+"/comp_event_tracks.txt",nevts=NEVENTS,
                        compareSmallRegion=False, compareEvent=True, tracksOnly=True, summary=False)
    DumpEventComparison(parser, em_events,sim_events,dname+"/comp_smallregion_tracks.txt",nevts=NEVENTS,
                        compareSmallRegion=True, compareEvent=False, tracksOnly=True, summary=False)

    DumpEventComparison(parser, em_events,sim_events,dname+"/summary_comp_event_all.txt",nevts=NEVENTS,
                        compareSmallRegion=False, compareEvent=True, tracksOnly=False, summary=True)
    if parser.all_dumps:
        DumpEventComparison(parser, em_events,sim_events,dname+"/comp_event_all.txt",nevts=NEVENTS,
                            compareSmallRegion=False, compareEvent=True, tracksOnly=False, summary=False)
        DumpEventComparison(parser, em_events,sim_events,dname+"/comp_smallregion_all.txt",nevts=NEVENTS,
                            compareSmallRegion=True, compareEvent=False, tracksOnly=False, summary=False)

    if parser.root:
        HistogramEventComparisons(parser, em_events,sim_events,dname+"/verification_hists.root")


    if False: # quick hack to make plots
        ctre=0
        ctrs=0
        print("Emulation")
        for e in em_events:
            if ctre==1: break
            ctre += 1
            for l in e.large_regions:
                for s in l.subregions:
                    # print("EM SR",ctrs)
                    # print("  ntracks",s.ntracks)
                    # print("  nems",s.nems)
                    # print("  ncalos",s.ncalos)
                    # print("  nmus",s.nmus)
                    ctrs += 1
                    print('h["em_sr_track"].SetBinContent({},{})'.format(ctrs,s.ntracks))
                    print('h["em_sr_em"].SetBinContent({},{})'.format(ctrs,s.nems))
                    print('h["em_sr_calo"].SetBinContent({},{})'.format(ctrs,s.ncalos))
                    print('h["em_sr_mu"].SetBinContent({},{})'.format(ctrs,s.nmus))
        ctre=0
        ctrs=0
        print("Simulation")
        for e in sim_events:
            ctre += 1
            for l in e.large_regions:
                for s in l.subregions:
                    #print("SIM SR",ctrs)
                    ctrs += 1
                    print('h["sim_sr_track"].SetBinContent({},{})'.format(ctrs,s.ntracks))
                    print('h["sim_sr_em"].SetBinContent({},{})'.format(ctrs,s.nems))
                    print('h["sim_sr_calo"].SetBinContent({},{})'.format(ctrs,s.ncalos))
                    print('h["sim_sr_mu"].SetBinContent({},{})'.format(ctrs,s.nmus))
                    # print("  ntracks",s.ntracks)
                    # print("  nems",s.nems)
                    # print("  ncalos",s.ncalos)
                    # print("  nmus",s.nmus)
        
        if 0: #debugging
            for e in sim_events:
                print('new event with {} LRs'.format(len(e.large_regions)))
                for l in e.large_regions:
                    print(' new LR with {} SRs'.format(len(l.subregions)),l.reg_type)
                    for s in l.subregions:
                        print('  new SR {}'.format(s.ID),s.reg_type)
                        print('  ',len(s.tracks))
                        print('  ',len(s.ems   ))
                        print('  ',len(s.calos ))
                        print('  ',len(s.mus   ))
    #ComparePerEvent(em_events,sim_events,nevts=6)
    #CheckAllObjects(em_events,sim_events)

    #TrackEmInSim(em_events,sim_events)
    #exit(0)

    # for em_track in em_tracks:
    #     print("em_tracks", em_track,"\n")
    # print("em_tracks", em_tracks,"\n")
    # print("em_ems   ", em_ems   ,"\n")
    # print("em_calos ", em_calos ,"\n")
    # print("em_mus   ", em_mus   ,"\n")
    # print(len(em_tracks),len(em_tracks[0]))
    # sim_tracks, sim_ems, sim_calos, sim_mus = GetSimulationData(parser)
    # print("sim_tracks", sim_tracks,"\n")
    # exit(0)

    # print event-by-event comparisons
    #ComparePerEvent(em_tracks,sim_tracks)
    # cumulative over all events ()
    #CompareOverlap(em_events,sim_events)
    exit(0)
    # investigate common tracks
    CheckCommonTracks(em_tracks,sim_tracks)

    # dump a number of distributions to a root file
    #DumpHistograms(em_tracks,sim_tracks)

    # investigate the tracks that are sent through the input links
    #DumpInputs(in_tracks, in_ems, in_calos)


                
    # # truncate according to the emulator setup
    # NTRACK=15
    # NEM=15
    # NCALO=15
    # for ii in range(CLKMAX):
    #     print
    #     print("Event",ii)
    #     print("Tracks- EMU:", list(map(hex,  em_tracks[ii][:NTRACK])))
    #     print("        SIM:", list(map(hex, sim_tracks[ii][:NTRACK])))
    #     print("EM    - EMU:", list(map(hex,  em_ems   [ii][:NEM])))
    #     print("        SIM:", list(map(hex, sim_ems   [ii][:NEM])))
    #     print("CALO  - EMU:", list(map(hex,  em_calos [ii][:NCALO])))
    #     print("        SIM:", list(map(hex, sim_calos [ii][:NCALO])))

    


    # print(em_tracks[5])
    # print(em_ems[5])
    # print(em_calos[5])
    # print(sim_tracks[5])
    # print(sim_ems[5])
    # print(sim_calos[5])
