#ifndef SIMPLE_PFLOW_SIMPLE_VTX_H
#define SIMPLE_PFLOW_SIMPLE_VTX_H

#include "../../firmware/data.h"

void simple_vtx_hwopt(TkObj track0[NALLTRACK], VtxObj *outvtx);
void simple_vtx_ref  (TkObj track [NALLTRACK], VtxObj *outch);

#endif
