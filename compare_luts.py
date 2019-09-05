import sys
import os


emu_file='lut.txt'
sim_file='/home/therwig/sandbox/otsdaq-cms-firmware/regionizer_full_sw/gen_level1_lookup_table/tNegBarrel.txt'
#sim_file='/home/therwig/sandbox/otsdaq-cms-firmware/regionizer_full_sw/gen_level1_lookup_table/tPosBarrel.txt'

emu_i=[]
emu_eta=[]
emu_phi=[]
emu_eta_overlap=[]
emu_phi_overlap=[]

sim_i=[]
sim_eta=[]
sim_phi=[]
sim_eta_overlap=[]
sim_phi_overlap=[]

with open(emu_file,'r') as f:
    ctr=0
    for l in f:
        arr = l.split()
        emu_i.append( int(arr[0]) )
        emu_eta.append( int(arr[2]) )
        emu_eta_overlap.append( int(arr[3]) )
        emu_phi.append( int(arr[5]) )
        emu_phi_overlap.append( int(arr[6]) )

with open(sim_file,'r') as f:
    ctr=0
    for l in f:
        if len(l)==0: break
        arr = l.split()
        if len(arr)<6: break
        sim_i.append( int(arr[0]) )
        sim_eta.append( int(arr[2]) )
        sim_eta_overlap.append( int(arr[3]) )
        sim_phi.append( int(arr[5]) )
        sim_phi_overlap.append( int(arr[6]) )

assert(len(emu_i)           == len(sim_i))
assert(len(emu_eta)         == len(sim_eta))
assert(len(emu_phi)         == len(sim_phi))
assert(len(emu_eta_overlap) == len(sim_eta_overlap))
assert(len(emu_phi_overlap) == len(sim_phi_overlap))

print("Comparing entry differences for emulation / sim ")
for i in range(len(emu_i)):
    if(emu_i[i]           != sim_i[i]
       or ((emu_eta[i]    != sim_eta[i]) and not (emu_eta[i] ==-1 and sim_eta[i]==2))
       or emu_phi[i]      != sim_phi[i]
       or ((emu_eta_overlap[i] != sim_eta_overlap[i]) and not (emu_eta_overlap[i] ==-1 and sim_eta_overlap[i]==2))
       or emu_phi_overlap[i] != sim_phi_overlap[i]
    ):

        ignore=False
        # haven't implemented proper eta/phi regions for emulation in the overlap case
        # if (emu_eta_overlap[i] and sim_eta_overlap[i]) and (emu_eta[i] != sim_eta[i]): ignore=True
        # if (emu_phi_overlap[i] and sim_phi_overlap[i]) and (emu_phi[i] != sim_phi[i]): ignore=True
        # temporarily ignore eta
        # if(emu_eta[i] != sim_eta[i] or emu_eta_overlap[i] != sim_eta_overlap[i]): ignore=True
        # temporarily ignore phi
        #if(emu_phi[i] != sim_phi[i] or emu_phi_overlap[i] != sim_phi_overlap[i]): ignore=True
 
        # different junk values written in usued case
        if((emu_eta[i] ==-1 and sim_eta[i]==2) or emu_eta_overlap[i] != sim_eta_overlap[i]): ignore=True
        if ignore: continue
        if emu_i[i] != sim_i[i]: print("Index mismatch {} vs {} !!!".format(emu_i[i], sim_i[i]))

        print("Index {}: Eta {}/{} ({}/{}), Phi {}/{} ({}/{}) --> ".format( emu_i[i],
                                                                            emu_eta[i],         sim_eta[i],        
                                                                            emu_eta_overlap[i], sim_eta_overlap[i], 
                                                                            emu_phi[i],         sim_phi[i],        
                                                                            emu_phi_overlap[i], sim_phi_overlap[i]),
              end='')
        if(emu_i[i]           != sim_i[i]): print(" INDEX",end='')
        if(emu_eta[i]         != sim_eta[i]): print(" ETA",end='')
        if(emu_eta_overlap[i] != sim_eta_overlap[i]): print(" ETAOLP",end='')
        if(emu_phi[i]         != sim_phi[i]): print(" PHI",end='')
        if(emu_phi_overlap[i] != sim_phi_overlap[i]): print(" PHIOLP",end='')
        print("")
        

# for i in range(len(emu_i)):
#     print(emu_i[i], emu_eta[i], emu_eta_overlap[i], emu_phi[i], emu_phi_overlap[i])
