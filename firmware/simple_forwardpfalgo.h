#ifndef SIMPLE_PFALGO3_H
#define SIMPLE_PFALGO3_H

#include "data.h"

bool match_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, etaphi_t boxSize) ;
etaphi_t dr_box(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) ;
template<int NB> ap_uint<NB>  dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) ;

void pfalgo3_forward_ref(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
void pfalgo3_forward_ref_set_debug(bool debug);
void pfalgo3_forward(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU], em2calo_dr_t drvals_em2calo[NPHOTON][NSELCALO]) ;
//void pfalgo3_forward(EmCaloObj calo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
void mp7wrapped_pack_in(EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) ;
void mp7wrapped_unpack_in(MP7DataWord data[MP7_NCHANN], EmCaloObj emcalo[NEMCALO], HadCaloObj hadcalo[NCALO], MuObj mu[NMU]) ;
void mp7wrapped_pack_out(PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU], MP7DataWord data[MP7_NCHANN]) ;
void mp7wrapped_pack_out_necomb(PFNeutralObj outne_all[NNEUTRALS], PFChargedObj outmu[NMU], MP7DataWord data[MP7_NCHANN]) ;
void mp7wrapped_unpack_out(MP7DataWord data[MP7_NCHANN], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
void mp7wrapped_unpack_out_necomb(MP7DataWord data[MP7_NCHANN], PFNeutralObj outpho[NPHOTON], PFNeutralObj outne[NSELCALO], PFChargedObj outmu[NMU]) ;
void mp7wrapped_pfalgo3_forward(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN]) ;

#endif

#ifndef DRVALSET
#define DRVALSET
#define PFALGO3_DR2MAX_EM_CALO 525
#define PFPUPPI_DR2MAX 8405

//#define PFALGO3_DR2MAX_EM_CALO 23

#endif
