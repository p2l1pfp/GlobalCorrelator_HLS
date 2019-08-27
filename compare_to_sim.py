import sys
import os
import argparse

# optionally dump some kinematics
plotROOT=False

#CLKMAX=150
#CLKMAX=108
CLKMAX=18
#CLKMAX=1

def GetEmulationData(parser):
  
    NTRACK=15
    NEM=15
    NCALO=15

    tracks=[]
    ems=[]
    calos=[]
    with open(parser.emulator_output,'r') as f:
        ctr=0
        for l in f:
            arr = list(map(lambda x : int(x,16), l.split()))
            step = arr[0]
            tracks += [arr[1:NTRACK+1]]
            ems    += [arr[1+NTRACK:NTRACK+NEM+1]]
            calos  += [arr[1+NTRACK+NEM:NTRACK+NEM+NCALO+1]]
            ctr+=1
            if ctr>CLKMAX:
                #print("em: reset clock max")
                pass
            # print( (tracks), (ems), (calos) )

    return tracks, ems, calos

def GetSimulationData(parser):
    NEM=15
    NCALO=20
    NTRACK=25

    ems=[]
    calos=[]
    tracks=[]
    for i in range(CLKMAX):
        ems.append([0]*NEM)
        calos.append([0]*NCALO)
        tracks.append([0]*NTRACK)

    for ii in range(NEM+NCALO+NTRACK):
        with open("{}/sim_HLS_input_object_{}.dat".format(parser.sim_output_dir,ii),'r') as f:
            ctr=0
            for l in f:
                val = int(l,base=16)
                if ctr >= CLKMAX:
                    #print("sim: reset clock max")
                    continue
                    #break
                if ii<NEM:
                    #print(ctr,ii)
                    ems[ctr][ii]=val
                elif ii<NEM+NCALO:
                    calos[ctr][ii-NEM]=val
                elif ii<NEM+NCALO+NTRACK:
                    tracks[ctr][ii-NEM-NCALO]=val
                # if ctr==5: 
                #     print(l)
                ctr += 1

    # hack to shorten the max # of sim outputs to match the emulation
    for x in ems: x = x[:15]
    for x in calos: x = x[:15]
    for x in tracks: x = x[:15]

    return tracks, ems, calos


def SelectBits(x, nbits, shift):
    return ((2**nbits-1 << shift) & x) >> shift
def BitsToInt(x,nbits):
    # max_uint=2**nbits-1
    # min_uint=0
    # max_int=2**(nbits-1)-1
    # min_int=-2**(nbits-1)
    # if x <= max_int: return x
    # else: return x - max_uint
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
    #print("track: pt={}, pte={}, eta={}, phi={}, z0={}, qual={} ({})".format(pt, pte, eta, phi, z0, qual, x))
    return pt, pte, eta, phi, z0, qual

def GetCaloParams(x):
    pt   =  SelectBits(x,16, 0)
    empt =  SelectBits(x,16,16)
    eta  =  SelectBits(x,10,32)#-2**10 # allow negative
    phi  =  SelectBits(x,10,42)
    isEM =  SelectBits(x, 1,52)
    #print("calo: pt={}, empt={}, eta={}, phi={}, isEM={} ({})".format(pt, empt, eta, phi, isEM, x))
    return pt, empt, eta, phi, isEM

def GetEMParams(x):
    pt   =  SelectBits(x,16, 0)
    pte  =  SelectBits(x,16,16)
    eta  =  SelectBits(x,10,32)#-2**10 # allow negative
    phi  =  SelectBits(x,10,42)
    #print("EM: pt={}, pte={}, eta={}, phi={} ({})".format(pt, pte, eta, phi, x))
    return pt, pte, eta, phi


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

if __name__ == "__main__":

    # GetTrackParams(0x40411333000B0012)
    # GetCaloParams(0x40411333000B0012)
    # GetEMParams(0x40411333000B0012)
    # exit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="",  dest = "input_file", help = "regionizer input from emulator")
    parser.add_argument("-e", "--emulation-output", default="output.txt",  dest = "emulator_output", help = "regionizer output from emulator")
    parser.add_argument("-s", "--simulation-data-dir", default="/home/therwig/sandbox/otsdaq-cms-firmware/regionizer_full/sim/sim_data/", dest = "sim_output_dir", 
                        help = "regionizer output directory from simulation")
    parser = parser.parse_args(sys.argv[1:])
    print("Reading from emulation output: "+parser.emulator_output)
    print("Reading from simulation output: "+parser.sim_output_dir)
    if len(parser.input_file): print("Checking inputs from: "+parser.input_file)

    in_tracks, in_ems, in_calos = GetInputs(parser)
    em_tracks, em_ems, em_calos = GetEmulationData(parser)
    sim_tracks, sim_ems, sim_calos = GetSimulationData(parser)

    em_tups=[]
    sim_tups=[]

    all_em=[]
    all_sim=[]
    
    if plotROOT:
        import ROOT
        f = ROOT.TFile("out.root","RECREATE")
        hPhi_EM = ROOT.TH1D("hPhi_EM","",40,0,2**10)
        hEta_EM = ROOT.TH1D("hEta_EM","",40,0,2**10)
        hPhi_SIM = ROOT.TH1D("hPhi_SIM","",40,0,2**10)
        hEta_SIM = ROOT.TH1D("hEta_SIM","",40,0,2**10)
        hPhi_EMonly = ROOT.TH1D("hPhi_EMonly","",40,0,2**10)
        hEta_EMonly = ROOT.TH1D("hEta_EMonly","",40,0,2**10)
    
    for ii in range(CLKMAX):
        print
        print("Event {}".format(ii))
        # remove zeros
        em_tracks[ii] = [x for x in em_tracks[ii] if x]
        sim_tracks[ii] = [x for x in sim_tracks[ii] if x]
        # print(em_tracks[ii])
        # print(sim_tracks[ii])

        # Compare tracks
        common_tks = set(em_tracks[ii]).intersection(set(sim_tracks[ii]))
        em_only = set(em_tracks[ii]).difference(common_tks)
        sim_only = set(sim_tracks[ii]).difference(common_tks)
        print("{} common, {} only emulation, {} only simulation".format(len(common_tks),len(em_only),len(sim_only)))
        for tk in em_tracks[ii]:
            pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
            if tk:
                all_em.append(tk)
                print("EM tk",pt,eta,phi,"{:0>16x}".format(tk))
                em_tups.append( (pt,eta,phi) )
                if plotROOT:
                    hPhi_EM.Fill(phi)
                    hEta_EM.Fill(eta)
        for tk in sim_tracks[ii]:
            pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
            if tk:
                all_sim.append(tk)
                print("SIM tk",pt,eta,phi,"{:0>16x}".format(tk))
                sim_tups.append( (pt,eta,phi) )
                if plotROOT:
                    hPhi_SIM.Fill(phi)
                    hEta_SIM.Fill(eta)

    common_tks = set(all_em).intersection(set(all_sim))
    em_only = set(all_em).difference(common_tks)
    sim_only = set(all_sim).difference(common_tks)
    print("\nOVERALL: {} common, {} only emulation, {} only simulation".format(len(common_tks),len(em_only),len(sim_only)))
    for tk in em_only:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        print("Emulation-only track",pt,eta,phi, "{:0>16x}".format(tk))
        if plotROOT:
            hPhi_EMonly.Fill(phi)
            hEta_EMonly.Fill(eta)
    for tk in common_tks:
        pt, pte, eta, phi, z0, qual = GetTrackParams(tk)
        print("Common tracks", pt, eta, phi, "{:0>16x}".format(tk))


    if plotROOT:
        hPhi_EM.Write()
        hEta_EM.Write()
        hPhi_SIM.Write()
        hEta_SIM.Write()
        hPhi_EMonly.Write()
        hEta_EMonly.Write()
        f.Close()
                
        
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
        #print("Input tracks", pt, eta, phi, "{:0>16x}".format(tk))

    # do we find EM only in the inputs?
    all_in = set(in_tracks)
    em_no_input = em_only.difference(all_in)
    print("We find a total of {} input tracks ({} unique),".format(len(in_tracks),len(all_in)))
    print(" and we find that {} of the {} emulated-but-not-simulated tracks".format(len(em_no_input),len(em_only)))
    print(" do not match any input track!!")


    # print(len(em_tups))
    # print(len(set(em_tups)))
    # print(set(em_tups))

    # print(len(sim_tups))
    # print(len(set(sim_tups)))
    # print(set(sim_tups))
                
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
