import sys
import os
import argparse

# CLKMAX=150
# CLKMAX=108
CLKMAX=18
    


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
                print("em: reset clock max")
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

def PrintTrack(x):
    # MSB -> LSP
    # pt_mask = 2**16-1
    # pte_mask = 2**16-1 << 16
    # eta_mask = 2**10-1 << 32
    # phi_mask = 2**10-1 << 42
    # z0_mask = 2**10-1 << 52
    # qual_mask = 2**1-1 << 62

    pt   =  SelectBits(x,16, 0)
    pte  =  SelectBits(x,16,16)
    eta  =  SelectBits(x,10,32)-2**10 # allow negative
    phi  =  SelectBits(x,10,42)
    z0   =  SelectBits(x,10,52)
    qual =  SelectBits(x, 1,62)

    print("pt   = ",pt  )
    print("pte  = ",pte )
    print("eta  = ",eta ) # allow negative
    print("phi  = ",phi )
    print("z0   = ",z0  )
    print("qual = ",qual)
    return

if __name__ == "__main__":

    # PrintTrack(0x40411333000B0012)
    # exit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--emulation-output", default="output.txt",  dest = "emulator_output", help = "regionizer output from emulator")
    parser.add_argument("-s", "--simulation-data-dir", default="/home/therwig/sandbox/otsdaq-cms-firmware/regionizer_full/sim/sim_data/", dest = "sim_output_dir", 
                        help = "regionizer output directory from simulation")
    parser = parser.parse_args(sys.argv[1:])
    print("Reading from emulation output:",parser.emulator_output)
    print("Reading from simulation output:",parser.sim_output_dir)

    em_tracks, em_ems, em_calos = GetEmulationData(parser)
    sim_tracks, sim_ems, sim_calos = GetSimulationData(parser)

    # truncate according to the emulator setup
    NTRACK=15
    NEM=15
    NCALO=15

    for ii in range(CLKMAX):
        print("Event",ii)
        print("Tracks- EMU:", list(map(hex,  em_tracks[ii][:NTRACK])))
        print("        SIM:", list(map(hex, sim_tracks[ii][:NTRACK])))
        print("EM    - EMU:", list(map(hex,  em_ems   [ii][:NEM])))
        print("        SIM:", list(map(hex, sim_ems   [ii][:NEM])))
        print("CALO  - EMU:", list(map(hex,  em_calos [ii][:NCALO])))
        print("        SIM:", list(map(hex, sim_calos [ii][:NCALO])))
        print()

    


    # print(em_tracks[5])
    # print(em_ems[5])
    # print(em_calos[5])
    # print(sim_tracks[5])
    # print(sim_ems[5])
    # print(sim_calos[5])
