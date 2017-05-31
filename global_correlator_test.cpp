/*******************************************************************************
 Vendor: Xilinx 
Associated Filename: hier_func_test.c
Purpose:Vivado HLS Coding Style example 
Device: All 
Revision History: May 30, 2008 - initial release
                                                
*******************************************************************************
#-  (c) Copyright 2011-2016 Xilinx, Inc. All rights reserved.
#-
#-  This file contains confidential and proprietary information
#-  of Xilinx, Inc. and is protected under U.S. and
#-  international copyright and other intellectual property
#-  laws.
#-
#-  DISCLAIMER
#-  This disclaimer is not a license and does not grant any
#-  rights to the materials distributed herewith. Except as
#-  otherwise provided in a valid license issued to you by
#-  Xilinx, and to the maximum extent permitted by applicable
#-  law: (1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND
#-  WITH ALL FAULTS, AND XILINX HEREBY DISCLAIMS ALL WARRANTIES
#-  AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, INCLUDING
#-  BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-
#-  INFRINGEMENT, OR FITNESS FOR ANY PARTICULAR PURPOSE; and
#-  (2) Xilinx shall not be liable (whether in contract or tort,
#-  including negligence, or under any other theory of
#-  liability) for any loss or damage of any kind or nature
#-  related to, arising under or in connection with these
#-  materials, including for any direct, or any indirect,
#-  special, incidental, or consequential loss or damage
#-  (including loss of data, profits, goodwill, or any type of
#-  loss or damage suffered as a result of any action brought
#-  by a third party) even if such damage or loss was
#-  reasonably foreseeable or Xilinx had been advised of the
#-  possibility of the same.
#-
#-  CRITICAL APPLICATIONS
#-  Xilinx products are not designed or intended to be fail-
#-  safe, or for use in any application requiring fail-safe
#-  performance, such as life-support or safety devices or
#-  systems, Class III medical devices, nuclear facilities,
#-  applications related to the deployment of airbags, or any
#-  other applications that could lead to death, personal
#-  injury, or severe property or environmental damage
#-  (individually and collectively, "Critical
#-  Applications"). Customer assumes the sole risk and
#-  liability of any use of Xilinx products in Critical
#-  Applications, subject only to applicable laws and
#-  regulations governing limitations on product liability.
#-
#-  THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS
#-  PART OF THIS FILE AT ALL TIMES. 
#- ************************************************************************


This file contains confidential and proprietary information of Xilinx, Inc. and 
is protected under U.S. and international copyright and other intellectual 
property laws.

DISCLAIMER
This disclaimer is not a license and does not grant any rights to the materials 
distributed herewith. Except as otherwise provided in a valid license issued to 
you by Xilinx, and to the maximum extent permitted by applicable law: 
(1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX 
HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, 
INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR 
FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether 
in contract or tort, including negligence, or under any other theory of 
liability) for any loss or damage of any kind or nature related to, arising under 
or in connection with these materials, including for any direct, or any indirect, 
special, incidental, or consequential loss or damage (including loss of data, 
profits, goodwill, or any type of loss or damage suffered as a result of any 
action brought by a third party) even if such damage or loss was reasonably 
foreseeable or Xilinx had been advised of the possibility of the same.

CRITICAL APPLICATIONS
Xilinx products are not designed or intended to be fail-safe, or for use in any 
application requiring fail-safe performance, such as life-support or safety 
devices or systems, Class III medical devices, nuclear facilities, applications 
related to the deployment of airbags, or any other applications that could lead 
to death, personal injury, or severe property or environmental damage 
(individually and collectively, "Critical Applications"). Customer assumes the 
sole risk and liability of any use of Xilinx products in Critical Applications, 
subject only to applicable laws and regulations governing limitations on product 
liability. 

THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT 
ALL TIMES.

*******************************************************************************/
#include "global_correlator.h"
#include <fstream>
#include <iostream>


int main() {

	// 2D Array Version
	/////////////////////////

        int10 calos_pt[ieta_max][iphi_max] = {0};
        int10 calos_h_over_e[ieta_max][iphi_max] = {0};
	std::ifstream read_calo_file("tb_data/in_calo.dat", std::ifstream::in);
	std::cout << "calo in: " << read_calo_file.is_open() << std::endl;
	std::string read_calo_object;
	int read_calo_ieta, read_calo_iphi, read_calo_pt, read_calo_h_over_e;
	while(read_calo_file>>read_calo_object>>read_calo_ieta>>read_calo_iphi>>read_calo_pt>>read_calo_h_over_e){
		calos_pt[read_calo_ieta][read_calo_iphi] = read_calo_pt + 1;
		calos_h_over_e[read_calo_ieta][read_calo_iphi] = read_calo_h_over_e;
	}
	read_calo_file.close();
  
	int10 tracks_pt[ieta_max][iphi_max] = {0};
	std::ifstream read_track_file("tb_data/in_track.dat", std::ifstream::in);
	std::cout << "track in: " << read_track_file.is_open() << std::endl;
	std::string read_track_object;
	int read_track_ieta, read_track_iphi, read_track_pt;
	while(read_track_file>>read_track_object>>read_track_ieta>>read_track_iphi>>read_track_pt){
	    tracks_pt[read_track_ieta][read_track_iphi] = read_track_pt + 4;
	}
	read_track_file.close();


	int11 pf_neutral_pt[ieta_max][iphi_max] = {0};
	calo_track_linking_grid(calos_pt, calos_h_over_e, tracks_pt, pf_neutral_pt);

	std::ofstream out_pf("tb_data/out_pf.dat", std::ofstream::out);
	std::cout << "out pf " << out_pf.is_open() << std::endl;
        for(int e=0; e<ieta_max; e++){
	  for(int p=0; p<iphi_max; p++){
	    out_pf << pf_neutral_pt[e][p] << std::endl;
	  }
        }
        out_pf.close();


	//List version
	//////////////////////////////

	int10 track_list_pt[ntrack_max] = {0};
	int10 track_list_ieta[ntrack_max] = {0};
	int10 track_list_iphi[ntrack_max] = {0};
	int10 calo_list_pt[ncalo_max] = {0};
	int10 calo_list_h_over_e[ncalo_max] = {0};
	int10 calo_list_ieta[ncalo_max] = {0};
	int10 calo_list_iphi[ncalo_max] = {0};

	std::ifstream read_calo_list_file("tb_data/in_calo_list.dat");
	std::string read_calo_list_object;
	int read_calo_list_ieta, read_calo_list_iphi, read_calo_list_pt, read_calo_list_h_over_e;
	int i=-1;
	while(read_calo_list_file>>read_calo_list_object>>read_calo_list_ieta>>read_calo_list_iphi>>read_calo_list_pt>>read_calo_list_h_over_e){
		i++;
		if(i>=ncalo_max) break;
		calo_list_pt[i]=read_calo_list_pt;
		calo_list_h_over_e[i]=read_calo_list_h_over_e;
		calo_list_ieta[i]=read_calo_list_ieta;
		calo_list_ieta[i]=read_calo_list_iphi;
	}
	read_calo_list_file.close();

	std::ifstream read_track_list_file("tb_data/in_track_list.dat");
	std::string read_track_list_object;
	int read_track_list_ieta, read_track_list_iphi, read_track_list_pt;
	i=-1;
	while(read_track_list_file>>read_track_list_object>>read_track_list_ieta>>read_track_list_iphi>>read_track_list_pt){
		i++;
		if(i>=ntrack_max) break;
		track_list_pt[i]=read_track_list_pt;
		track_list_ieta[i]=read_track_list_ieta;
		track_list_ieta[i]=read_track_list_iphi;
	}
	read_track_list_file.close();

	int10 pf_neutral_list_pt[npf_neutral_max] = {0};
	int10 pf_neutral_list_ieta[npf_neutral_max] = {0};
	int10 pf_neutral_list_iphi[npf_neutral_max] = {0};

	calo_track_linking_list(track_list_pt, track_list_ieta, track_list_iphi,
							calo_list_pt, calo_list_h_over_e, calo_list_ieta, calo_list_iphi,
							pf_neutral_list_pt, pf_neutral_list_ieta, pf_neutral_list_iphi);

	return 0;
}

// XSIP watermark, do not delete 67d7842dbbe25473c3c32b93c0da8047785f30d78e8a024de1b57352245f9689
