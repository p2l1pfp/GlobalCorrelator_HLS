#include "global_correlator.h"

void calo_track_linking_list(int10 track_list_pt[ntrack_max],
							 int10 track_list_ieta[ntrack_max],
							 int10 track_list_iphi[ntrack_max],
							 int10 calo_list_pt[ncalo_max],
							 int10 calo_list_h_over_e[ncalo_max],
							 int10 calo_list_ieta[ncalo_max],
							 int10 calo_list_iphi[ncalo_max],
							 int10 pf_neutral_list_pt[npf_neutral_max],
							 int10 pf_neutral_list_ieta[npf_neutral_max],
							 int10 pf_neutral_list_iphi[npf_neutral_max]){

	//algo here

	//track loop is serial
	for(int t=0; t<ntrack_max; t++){

		//calo loop is parallel (pragma in first line)
		for(int c=0; c<ncalo_max; c++){
#pragma HLS UNROLL
			//int10 delta_pt = track_list_pt[t] - calo_list_pt[c];
		}//calo loop
	}//track loop

}


void calo_track_linking_grid(int10 calos_pt[ieta_max][iphi_max],
			     int10 calos_h_over_e[ieta_max][iphi_max],
			     int10 tracks_pt[ieta_max][iphi_max],
			     int10 pf_neutral_pt[ieta_max][iphi_max]){

#pragma HLS ARRAY_PARTITION variable=calos_pt complete dim=0
#pragma HLS ARRAY_PARTITION variable=calos_h_over_e complete dim=0
#pragma HLS ARRAY_PARTITION variable=tracks_pt complete dim=0
#pragma HLS ARRAY_PARTITION variable=pf_neutral_pt complete dim=0

	for(int i=0; i<ieta_max; i++){
#pragma HLS UNROLL
		for(int j=0; j<iphi_max; j++){
//#pragma HLS UNROLL
			int10 delta_pt = calos_pt[i][j] - tracks_pt[i][j];
                        int delta_pt_int = calos_pt[i][j] - tracks_pt[i][j];
			std::cout << "calo " << calos_pt[i][j] << " - track " << tracks_pt[i][j] << std::endl;
			std::cout << "delta int: " << delta_pt_int << ", int10: " << delta_pt << std::endl;
			if(calos_h_over_e[i][j]>100){
				pf_neutral_pt[i][j] = delta_pt;
			}
			std::cout << "out " << calos_pt[i][j] << " " << tracks_pt[i][j] << " " << calos_h_over_e[i][j] << " " << pf_neutral_pt[i][j] << std::endl; 
		}
	}

}

