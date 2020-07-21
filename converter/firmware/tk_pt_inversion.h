
#include "tk_input_converter.h"

template<class pt_T> 
void init_pt_inv_table(pt_T table_out[(1<<PT_INV_TAB_SIZE)]);

template<class pt_inv_T, class pt_T> 
void convert_pt(pt_inv_T inv, pt_T &pt);
