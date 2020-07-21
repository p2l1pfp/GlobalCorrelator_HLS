#ifndef TEST_H
#define TEST_H
//
#include "ap_int.h"
#include "ap_fixed.h"
#include "../firmware/tk_input_converter.h"
#include <stdlib.h>
#include <fstream>
#include <string>

using std::cout;
using std::endl;

//test functions
void test_input_converter();
void test_tk_pack();
void test_pf_pack();
void test_tanlambda_to_eta();
void test_pt_inversion();
void test_prop_phi();
void test_prop_tanlamnda();
void test_resolution();

//helpers
l1tk_word_t get_random_l1track(unsigned int ntest);
inline float urand(float lo, float hi){return lo+(hi-lo)*float(rand())/RAND_MAX;}


#endif
