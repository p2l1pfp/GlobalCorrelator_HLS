import sys
import os
import argparse


TMUX_OUT=18
NEVENTS=1
#CLKMAX=150
#CLKMAX=108
CLKMAX=TMUX_OUT*NEVENTS
CLKMAX=140 # 134 lines in sim out

NTRACK=22
NEM=13
NCALO=15
NMU=2

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

class Event:
    def __init__(self, n_large = 1, n_small = 1):
        self.n_small_regions = n_small
        self.n_big_regions   = n_large
        self.large_regions   = []
        for i in range(self.n_big_regions):
            bigR = Region(n_sub = self.n_small_regions)
            bigR.reg_type = "large"
            self.large_regions.append(bigR)
    # def AddLargeRegion(self):
    #     pass #self.small_regions = []

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
    ordering = [(0,0), (0,2), (0,1), (0,3), (1,2), (0,4), (1,3), (0,5), (1,4), (0,6), (1,5), (0,7), (1,6), (0,8), (1,7), (1,0), (1,8), (1,1)]

    tracks=[]
    ems=[]
    calos=[]
    mus=[]

    # local vars
    EM_NEVENTS=6
    EM_NEVENTS=21
    EM_NLARGE=1
    EM_NSMALL=18

    evts = []

    # setup for the first region
    evt = Event(EM_NLARGE, EM_NSMALL)
    # large_reg = Region(n_sub=EM_NSMALL)
    # evt.large_regions.append( large_reg )

    # for l in evt.large_regions:
    #     print(' new LR with {} SRs'.format(len(l.subregions)),l.reg_type)
    #     for s in l.subregions:
    #         print('  new SR {}'.format(s.ID),s.reg_type)
            
    # exit(0)
    # for em_track in em_tracks:


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
            sr.ID     = step

            # increment counters
            # (add support for large region ctr later)
            ctr = ctr+1
            sreg_ctr = ctr % EM_NSMALL
            if ctr and (ctr % EM_NSMALL)==0: 
                evt_ctr = evt_ctr+1
                # create new event
                evts.append( evt )
                evt = Event(EM_NLARGE, EM_NSMALL)
                # large_reg = Region(n_sub=EM_NSMALL)
                # evt.large_regions.append( large_reg )


    # add empty events for now
    while len(evts) < EM_NEVENTS:
        evt = Event(EM_NLARGE, EM_NSMALL)
        evts.append( evt )
        
    #'serialized' list of 18 tmux regions * however many events...
    #return tracks, ems, calos, mus
    return evts

def GetSimulationData(parser):
    '''
    Simulation data is stored across multiple files,
    each corresponding to a time series of single-object outputs.
    #files is #em+#calo+#track, which may equal 15+20+25=60.
      --> 22+13+15+2=52?
    Each line corresonds to a 'step' of the output and...
    TODO add comment on length of each file.
    '''

    # local vars
    SIM_NEVENTS=6
    SIM_NEVENTS=21
    SIM_NLARGE=1
    SIM_NSMALL=18

    evts = []
    for i in range(SIM_NEVENTS):
        evts.append( Event(SIM_NLARGE, SIM_NSMALL) )

    sreg_ctr = 0
    lreg_ctr = 0
    evt_ctr  = 0

    for ii in range(NEM+NCALO+NTRACK+NMU):
        with open("{}/sim_HLS_input_object_{}.dat".format(parser.sim_output_dir,ii),'r') as f:
            ctr=0
            for l in f:
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
                        sr.ems.append(val)
                    elif ii<NEM+NCALO: 
                        sr.calos.append(val)
                    elif ii<NEM+NCALO+NTRACK:
                        sr.tracks.append(val)
                        sr.track_clks[val]=ctr
                    elif ii<NEM+NCALO+NTRACK+NMU:
                        sr.mus.append(val)
                ctr += 1

    return evts


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
    pt   =  SelectBits(x,16, 0)
    empt =  SelectBits(x,16,16)
    eta  =  SelectBits(x,10,32)#-2**10 # allow negative
    phi  =  SelectBits(x,10,42)
    isEM =  SelectBits(x, 1,52)
    return pt, empt, eta, phi, isEM
    #print("calo: pt={}, empt={}, eta={}, phi={}, isEM={} ({})".format(pt, empt, eta, phi, isEM, x))

def GetEMParams(x):
    pt   =  SelectBits(x,16, 0)
    pte  =  SelectBits(x,16,16)
    eta  =  SelectBits(x,10,32)#-2**10 # allow negative
    phi  =  SelectBits(x,10,42)
    return pt, pte, eta, phi
    #print("EM: pt={}, pte={}, eta={}, phi={} ({})".format(pt, pte, eta, phi, x))


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

def CompareOverlap(em_tracks,sim_tracks):
    '''
    Loop over each event / tmuxed region and dump a number of outputs into a histogram
    '''

    # add all non-zero entries to the total track lists
    all_em=[]
    all_sim=[]
    # for ii in range(CLKMAX):
    #     for tk in [x for x in em_tracks[ii] if x]: all_em.append(tk)
    #     for tk in [x for x in sim_tracks[ii] if x]: all_sim.append(tk)
    for em in em_tracks:
        #print(em)
        for tk in [x for x in em if x]: 
            all_em.append(tk)
            #print(hex(tk))
    for sim in sim_tracks:
        for tk in [x for x in sim if x]: all_sim.append(tk)

    exit(0)

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


if __name__ == "__main__":

    # GetTrackParams(0x40411333000B0012)
    # GetCaloParams(0x40411333000B0012)
    # GetEMParams(0x40411333000B0012)
    # exit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="",  dest = "input_file", help = "regionizer input from emulator")
    parser.add_argument("-e", "--emulation-output", default="outputs.txt",  dest = "emulator_output", help = "regionizer output from emulator")
    parser.add_argument("-s", "--simulation-data-dir", default="/home/therwig/sandbox/otsdaq-cms-firmware/regionizer_full_small/sim/sim_data/", dest = "sim_output_dir", 
                        help = "regionizer output directory from simulation")
    parser = parser.parse_args(sys.argv[1:])
    print()
    print("Reading from emulation output: "+parser.emulator_output)
    print("Reading from simulation output: "+parser.sim_output_dir)
    if len(parser.input_file): print("Checking inputs from: "+parser.input_file)

    if len(parser.input_file):
        in_tracks, in_ems, in_calos = GetInputs(parser)
        print (in_tracks)

    em_events  = GetEmulationData(parser)
    sim_events = GetSimulationData(parser)

    if 0: #debugging
        for e in em_events:
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

    TrackEmInSim(em_events,sim_events)
    exit(0)

    # for em_track in em_tracks:
    #     print("em_tracks", em_track,"\n")
    # print("em_tracks", em_tracks,"\n")
    # print("em_ems   ", em_ems   ,"\n")
    # print("em_calos ", em_calos ,"\n")
    # print("em_mus   ", em_mus   ,"\n")
    # print(len(em_tracks),len(em_tracks[0]))
    sim_tracks, sim_ems, sim_calos, sim_mus = GetSimulationData(parser)
    # print("sim_tracks", sim_tracks,"\n")
    # exit(0)

    # print event-by-event comparisons
    #ComparePerEvent(em_tracks,sim_tracks)
    # cumulative over all events ()
    CompareOverlap(em_tracks,sim_tracks)
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
