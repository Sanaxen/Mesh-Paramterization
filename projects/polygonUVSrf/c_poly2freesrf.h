#ifndef _C_POLY2FREESRF_H
#define _C_POLY2FREESRF_H

#ifdef __cplusplus
extern "C" {
#endif

void ReadMeshUV( char* filename, int* vtxnum, int* trinum, double** vertex, double** uv, unsigned int** tri );

int Polygon2SbezFitC(int* chkgridA, double *ctrlpnt, int upatch, int vpatch, double rectA[][2], 
					 int vtxnum, int trinum, double* vertex, double* uv, unsigned int* tri, double* normal, 
					double* start, double* step, double* (*fit)(int, double*) );

#ifdef __cplusplus
};
#endif

#endif
