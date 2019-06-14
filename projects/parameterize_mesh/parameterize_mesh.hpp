#ifndef _parameterize_mesh_hpp

#define _parameterize_mesh_hpp


#include "my_types.h"

#ifdef __cplusplus
extern "C" {
#endif
#define _DLLAPI extern "C" __declspec(dllexport) 

int cut_mesh(unsigned int mesh);
int parameterize_mesh(unsigned int mesh, int res, GIfloat stretch_weight);

#include "parameterize_mesh.h"
#include "poly2freesrf.h"
#include "c_poly2freesrf.h"

#ifdef __cplusplus
}
#endif


#endif