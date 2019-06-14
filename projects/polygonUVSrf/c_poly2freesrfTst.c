#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "c_poly2freesrf.h"

int GetGmProgress(){return 0;}
void gm_msg_loop(){return;}
gmProgressCreateForm(){return 1;}
gmProgressDeleteForm(void* p){return 1;}
gmProgressSetPosition(int i, void* p){return 1;}
gmProgressSetFormTitle(char* a, void* p){return 1;}
gmProgressSetMinMax(int i, int j, void* p){return 1;}
gmProgressShowForm(void* p){return 1;}

static FILE* DumpSrfStart(char* filename)
{
	FILE* fp = fopen(filename, "w");
	if ( fp == NULL ){
		return NULL;
	}
	fprintf(fp, "pconst_num(10,10);\n");

	return fp;
}
static void DumpSrfEnd(FILE* fp)
{
	if ( fp == NULL ){
		return;
	}
	fprintf(fp, "end;\n");
	fclose(fp);
}

static void DumpSrf(FILE*fp, int id, int un, int vn, double* cp)
{	


	int i, j, k;
/*
	k = 0;
	for ( i = 0; i < 3*vn+1; i++ ){
		for ( j = 0; j < 3*un+1; j++ ){
			fprintf(fp, "Pcp=[%f,%f,%f];\n", cp[3*k], cp[3*k+1], cp[3*k+2]);
			k++;
		}
	}
*/
	k = 0;
	fprintf(fp, "Sbez.%d = %d %d(\n", id, un, vn);
	for ( i = 0; i < 3*vn+1; i++ ){
		for ( j = 0; j < 3*un+1; j++ ){
			fprintf(fp, "[%f,%f,%f]", cp[3*k], cp[3*k+1], cp[3*k+2]);
			if ( i < 3*vn || j < 3*un){
				fprintf(fp, ",\n");
			}else{
				fprintf(fp, "\n");
			}
			k++;
		}
	}

	fprintf(fp, ");\n");
	fprintf(fp, "#uvmesh(120,120);\n");
	fprintf(fp, "#etrgmesh(Sbez.%d);\n", id);
}
int main(int argc, char** argv)
{
	int vtxnum;
	int trinum;
	double* vertex;
	double* uv;
	double* normal;
	int* tri;
	double* step;
	double* start;
	int upatch;
	int vpatch;
	int chkgrid[3];
	int grid[3];
	double* ctrlpnt;
    int stat;
	FILE* fp;

	//ReadMeshUV( "G:\\devwork\\Recon\\Recon.940503b\\src\\uvmap_mesh_maxplanck_optmz.dat", &vtxnum, &trinum, &vertex, &uv, &tri );
 	ReadMeshUV( argv[1], &vtxnum, &trinum, &vertex, &uv, &tri );



	normal = NULL;
	start = NULL;
	step = NULL;

	chkgrid[0] = 20;
	chkgrid[1] = 20;
	chkgrid[2] = 20;

	upatch = 45;
	vpatch = 45;

	grid[0] = 3*upatch+1;
	grid[1] = 3*vpatch+1;
	grid[2] = 1;

	ctrlpnt = (double*)malloc( 3*grid[0]*grid[1]*sizeof(double));
	if ( ctrlpnt == NULL ){
		return -1;
	}
	fp = DumpSrfStart("sbezer.txt");

	stat = Polygon2SbezFitC(chkgrid, ctrlpnt, upatch, vpatch, NULL,
					 vtxnum, trinum, vertex, uv, tri, normal, start, step, NULL);
	DumpSrf(fp, 0, upatch, vpatch, ctrlpnt);

	free(ctrlpnt);
	DumpSrfEnd(fp);
	return 0;
}
