#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "poly2freesrf.h"
#include "c_poly2freesrf.h"

#ifdef USE_GUI_DLL
extern "C" {
#include "ProgressDLL.h"
int GetGmProgress();
void gm_msg_loop();
};
#endif


void ReadMeshUV( char* filename, int* vtxnum, int* trinum, double** vertex, double** uv, unsigned int** tri )
{
	FILE* fp = fopen( filename, "r");
	if ( fp == NULL ){
		return;
	}
	char buf[256];

	fgets(buf, 256, fp);
	sscanf(buf, "%d", vtxnum);

	fgets(buf, 256, fp);
	sscanf(buf, "%d", trinum);

	*vertex = (double*)malloc( 3**vtxnum*sizeof(double));
	*uv = (double*)malloc( 3**vtxnum*sizeof(double));
	*tri = (unsigned int*)malloc( 3**trinum*sizeof(unsigned int));
	memset(*tri, '\0', 3**trinum*sizeof(int));

	for ( int i = 0; i < *vtxnum; i++ ){
		fgets(buf, 256, fp);
		sscanf(buf, "%lf %lf %lf %lf %lf",
			&((*vertex)[3*i]), &((*vertex)[3*i+1]), &((*vertex)[3*i+2]),
			&((*uv)[3*i]),	&((*uv)[3*i+1]) );
		(*uv)[3*i+2] = 0.0;
	}
	for ( int j = 0; j < *trinum; j++ ){
		fgets(buf, 256, fp);
		sscanf(buf, "%d %d %d",&((*tri)[3*j]), &((*tri)[3*j+1]), &((*tri)[3*j+2]));
	}
	fclose(fp);
}



static double vector_area(double v1[3], double v2[3] )
{
	double area;
	
	area = fabs(v1[0]*v2[1] - v1[1]*v2[0])*0.5;

	return area;
}

static void getTri(int id, double* vertex, unsigned int* triid, double vtx[3][3] ){
	int vid[3];

	vid[0] = triid[3*id  ];
	vid[1] = triid[3*id+1];
	vid[2] = triid[3*id+2];

	for (int i = 0; i < 3; i++ ){
		vtx[i][0] = vertex[3*vid[i]];
		vtx[i][1] = vertex[3*vid[i]+1];
		vtx[i][2] = vertex[3*vid[i]+2];
	}
}

static double lineParameter(double* stp, double* edp, double* pnt, double* onp, int* stat)
{
	double v0[3], v1[3];
	double  tt;

	v0[0] = edp[0] - stp[0];
	v0[1] = edp[1] - stp[1];
	v0[2] = edp[2] - stp[2];

	v1[0] = pnt[0] - stp[0];
	v1[1] = pnt[1] - stp[1];
	v1[2] = pnt[2] - stp[2];
	double wt = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2] );

	*stat = -2;
	if  ( wt >= 1.0e-16 )  {
		tt = (v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2] ) / wt;
		if ( onp != NULL ){
			onp[0] = stp[0] + tt * v0[0];
			onp[1] = stp[1] + tt * v0[1];
			onp[2] = stp[2] + tt * v0[2];
		}
		*stat = 0;
	}  else  {
		tt = 0.0;
		if ( onp != NULL ){
			onp[0] = stp[0];
			onp[1] = stp[1];
			onp[2] = stp[2];
		}
		*stat = -1;
	}
	return tt;
}

static void TrianglePointUV(double vtx0[3], double vtx1[3], double vtx2[3], double* uv, double* pnt)
{
	double xx[3], yy[3];

	xx[0] = vtx1[0] - vtx0[0];
	xx[1] = vtx1[1] - vtx0[1];
	xx[2] = vtx1[2] - vtx0[2];

	yy[0] = vtx2[0] - vtx0[0];
	yy[1] = vtx2[1] - vtx0[1];
	yy[2] = vtx2[2] - vtx0[2];

/*
	double t = uv[0]*uv[1];

	pnt[0] = (uv[1]-t)*xx[0] + t*yy[0] + vtx0[0];
	pnt[1] = (uv[1]-t)*xx[1] + t*yy[1] + vtx0[1];
	pnt[2] = (uv[1]-t)*xx[2] + t*yy[2] + vtx0[2];
*/
	pnt[0] = uv[0]*xx[0] + uv[1]*yy[0] + vtx0[0];
	pnt[1] = uv[0]*xx[1] + uv[1]*yy[1] + vtx0[1];
	pnt[2] = uv[0]*xx[2] + uv[1]*yy[2] + vtx0[2];
}

static int TriangleUV(double vtx0[3], double vtx1[3], double vtx2[3], double *pnt, double *uv)
{

	double w;
	double a,b,c;

	double xx[3], yy[3];

	xx[0] = vtx1[0] - vtx0[0];
	xx[1] = vtx1[1] - vtx0[1];
	xx[2] = vtx1[2] - vtx0[2];

	yy[0] = vtx2[0] - vtx0[0];
	yy[1] = vtx2[1] - vtx0[1];
	yy[2] = vtx2[2] - vtx0[2];

	w = xx[0]*yy[1] - yy[0]*xx[1];

	if ( w == 0.0){
		return 0;
	}


	a = yy[0]*(vtx0[1] - pnt[1]) + yy[1]*(-vtx0[0] + pnt[0]);

	b = xx[0]*(-vtx0[1] + pnt[1]) + xx[1]*(vtx0[0] - pnt[0]);

	uv[0] = a/w;
	uv[1] = b/w;
	return 1;
}


void PolygonSurface::Create(int chkgrid[3], int upatch, int vpatch, double rect[2][2],
			int vtxnum, int trinum, double* vertex, double* uv, unsigned int* tri, double* normal, 
			double* start, double* step)
{
	this->chkgrid_[0] = chkgrid[0];
	this->chkgrid_[1] = chkgrid[1];
	this->chkgrid_[2] = chkgrid[2];
	this->vtxnum_ = vtxnum;
	this->trinum_ = trinum;
	this->vertex_ = vertex;
	this->uv_ = uv;
	this->tri_ = tri;
	this->normal_ = normal;
	if ( start != NULL ){
		this->start_[0] = start[0];
		this->start_[1] = start[1];
		this->start_[2] = start[2];
	}
	if ( step != NULL ){
		this->step_[0] = step[0];
		this->step_[1] = step[1];
		this->step_[2] = step[2];
	}

	// UV側のメッシュでグリッド管理する
	// mgridは頂点の代わりにUV座標が対象となる
	this->mgrid_.Init( vtxnum_, trinum_, uv_, tri_);
	this->mgrid_.mesh_to_grid(chkgrid_);
}


int PolygonSurface::Evalue(double* uvparam, double point[3] )
{
	double uvpoint[3];
	double onpoint[3];

	uvpoint[0] = uvparam[0];
	uvpoint[1] = uvparam[1];
	uvpoint[2] = 0.0;
	onpoint[2] = 0.0;

	double dist;
	int		id;

	this->mgrid_.overlap_coef_ = 0.01;
	do{
		// UV格子上の点が乗る三角形idとUV三角形上のUV点を求める
		id = this->mgrid_.OnTri(uvpoint, dist, onpoint);
		if ( id < 0 ){
			this->mgrid_.overlap_coef_ *= 1.005;
			continue;
		}
		break;
	}while(this->mgrid_.overlap_coef_ < 2.0);

	if (id < 0 ){
		return -1;
	}
	//printf("id %d onpoint %f %f %f\n", id, onpoint[0], onpoint[1], onpoint[2]);

	//id から実際の３D 頂点を得る
	double vtx3[3][3];
	getTri(id, vertex_, tri_, vtx3 );
/*
	printf("vtx3[0] %f %f %f\n", vtx3[0][0], vtx3[0][1], vtx3[0][2]);
	printf("vtx3[1] %f %f %f\n", vtx3[1][0], vtx3[1][1], vtx3[1][2]);
	printf("vtx3[2] %f %f %f\n", vtx3[2][0], vtx3[2][1], vtx3[2][2]);
*/

	//id から実際の2D(UV) 頂点を得る
	double vtx[3][3];
	this->mgrid_.GetTri(id, vtx);
/*
	printf("vtx[0] %f %f %f\n", vtx[0][0], vtx[0][1], vtx[0][2]);
	printf("vtx[1] %f %f %f\n", vtx[1][0], vtx[1][1], vtx[1][2]);
	printf("vtx[2] %f %f %f\n", vtx[2][0], vtx[2][1], vtx[2][2]);
*/

	/*         vtx[2]
	                   b
	                 /    ＼
	                /       ＼
	               /   p      ＼
	              /             ＼ a
	            c ----------------  vtx[1] 
	            vtx[0]
	*/
	double tt[2];
	int	stat;

	// UV三角形上のUV点の2D(UV) 頂点の２辺におけるパラメータを求める
	/*
	tt[0] = lineParameter(vtx[0], vtx[1], onpoint, NULL, &stat);
	tt[1] = lineParameter(vtx[0], vtx[2], onpoint, NULL, &stat);
	*/
	stat = TriangleUV(vtx[0], vtx[1], vtx[2], onpoint, tt);
	if ( stat == 0 ){
		tt[0] = lineParameter(vtx[0], vtx[1], onpoint, NULL, &stat);
		tt[1] = lineParameter(vtx[0], vtx[2], onpoint, NULL, &stat);
	}

//	printf("tt %f %f\n", tt[0], tt[1] );
	// パラメータから実際の3D三角形上の点を求める
	TrianglePointUV(vtx3[0], vtx3[1], vtx3[2], tt, point);

#if 0
	printf("Line=[%f,%f],[%f,%f];\n", 
		vtx[0][0], vtx[0][1], vtx[1][0], vtx[1][1]);
	printf("Line=[%f,%f],[%f,%f];\n", 
		vtx[0][0], vtx[0][1], vtx[2][0], vtx[2][1]);
	printf("Line=[%f,%f],[%f,%f];\n", 
		vtx[1][0], vtx[1][1], vtx[2][0], vtx[2][1]);
	printf("Pon=[%f,%f];\n", onpoint[0], onpoint[1]);


	printf("Line=[%f,%f,%f],[%f,%f,%f];\n", 
		vtx3[0][0], vtx3[0][1], vtx3[0][2], vtx3[1][0], vtx3[1][1], vtx3[1][2] );
	printf("Line=[%f,%f,%f],[%f,%f,%f];\n", 
		vtx3[0][0], vtx3[0][1], vtx3[0][2], vtx3[2][0], vtx3[2][1], vtx3[2][2] );
	printf("Line=[%f,%f,%f],[%f,%f,%f];\n", 
		vtx3[1][0], vtx3[1][1], vtx3[1][2], vtx3[2][0], vtx3[2][1], vtx3[2][2] );

	printf("Pon=[%f,%f,%f];\n", point[0], point[1], point[2]);
	printf("\n");
#endif

	return 0;
}

void PolygonSurface::DumpUVLines(FILE* fp, int un, int vn)
{
	double uvpoint[3];
	double point[3];
	double pre_point[3];
	double step[2];
	int		ist;

	step[0] = 1.0/(double)(un-1);
	step[1] = 1.0/(double)(vn-1);

	int i, j;
#if 0
	for ( j = 0; j < vn; j++ ){
		for ( i = 0; i < un; i++ ){
			uvpoint[0] = (double)i*step[0];
			uvpoint[1] = (double)j*step[1];

			// パラメータから実際の3Dの点を求める
			ist = Evalue(uvpoint, point);
			if ( i == 0 ){
				pre_point[0] = point[0];
				pre_point[1] = point[1];
				pre_point[2] = point[2];
				continue;
			}
			fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n", 
						pre_point[0], pre_point[1], pre_point[2],
						point[0], point[1], point[2] );
			pre_point[0] = point[0];
			pre_point[1] = point[1];
			pre_point[2] = point[2];
		}
	}

	for ( j = 0; j < un; j++ ){
		for ( i = 0; i < vn; i++ ){
			uvpoint[0] = (double)j*step[0];
			uvpoint[1] = (double)i*step[1];

			// パラメータから実際の3Dの点を求める
			ist = Evalue(uvpoint, point);
			if ( i == 0 ){
				pre_point[0] = point[0];
				pre_point[1] = point[1];
				pre_point[2] = point[2];
				continue;
			}
			fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n", 
						pre_point[0], pre_point[1], pre_point[2],
						point[0], point[1], point[2] );
			pre_point[0] = point[0];
			pre_point[1] = point[1];
			pre_point[2] = point[2];
		}
	}
#else

	un = 100;
	step[0] = 1.0/(double)(un-1);

	for ( j = 0; j < un; j++ ){
		uvpoint[1] = (double)j*step[0];		
		uvpoint[0] = 0.5;

		// パラメータから実際の3Dの点を求める
		ist = Evalue(uvpoint, point);
		if ( j == 0 ){
			pre_point[0] = point[0];
			pre_point[1] = point[1];
			pre_point[2] = point[2];
			continue;
		}
		fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n", 
					pre_point[0], pre_point[1], pre_point[2],
					point[0], point[1], point[2] );
		pre_point[0] = point[0];
		pre_point[1] = point[1];
		pre_point[2] = point[2];
	}
#endif

	fprintf(fp, "end;\n");
}


double* Polygon2SbezFit( PolygonSurface* ps, int upatch, int vpatch, double rect[2][2], double* (*fit)(int, double*) )
{

	double	uvpoint[3];
	double	point[3];
	double* ctrlpnt;
	int		grid[3];

	ctrlpnt = NULL;

	uvpoint[0] = 0.0;
	uvpoint[1] = 0.0;
	uvpoint[2] = 0.0;
	point[0] = 0.0;
	point[1] = 0.0;
	point[2] = 0.0;

	grid[0] = 3*upatch+1;
	grid[1] = 3*vpatch+1;
	grid[2] = 1;

	ps->alcflg_ = 0;
	if ( ps->GetCPArea() == NULL ){
		ctrlpnt = (double*)malloc( 3*grid[0]*grid[1]*sizeof(double));
		if ( ctrlpnt == NULL ){
			return NULL;
		}
		ps->alcflg_ = 1;
	}else{
		ctrlpnt = ps->GetCPArea();
	}

	double* wkcp, *wkpnt;
	int	k_s, k_e;

	wkcp = NULL;
	wkpnt = NULL;
	if ( fit != NULL ){
		if ( grid[0] > grid[1] ){
			wkpnt = (double*)malloc(3*grid[0]*sizeof(double));
		}else{
			wkpnt = (double*)malloc(3*grid[1]*sizeof(double));
		}
	}

	double ru, rv;
	int i, j, k;
	int ii = 0;
	int jj = 0;
	int nn;
	k = 0;

       	printf("C0-リーゼンフェルトスプライン作成 START ");
#ifdef USE_GUI_DLL
	void* progress = NULL;
	if( GetGmProgress() == 1 ){
		progress = gmProgressCreateForm();
        	gmProgressSetFormTitle("C0-リーゼンフェルトスプライン作成", progress);
		gmProgressSetMinMax(0, grid[0]*grid[1], progress);
		gmProgressShowForm(progress);
	}		
#endif
	for ( j = 0; j < grid[1]; j++ ){
		for ( i = 0; i < grid[0]; i++ ){
			ru = (double)i/(grid[0]-1);
			rv = (double)j/(grid[1]-1);
			uvpoint[0] = (rect[0][1] - rect[0][0])*ru + rect[0][0]; 
			uvpoint[1] = (rect[1][1] - rect[1][0])*rv + rect[1][0];

			// パラメータから実際の3Dの点を求める
			int ist = ps->Evalue(uvpoint, point);
			if ( ist != 0 ){
			//	printf("i %d j %d error u %f v %f\n", i, j, uvpoint[0], uvpoint[1] );
			}
			//printf("i %d j %d   u %f v %f %f %f %f\n", i, j, uvpoint[0], uvpoint[1],point[0], point[1], point[2] );
			ctrlpnt[3*k  ] = point[0];
			ctrlpnt[3*k+1] = point[1];
			ctrlpnt[3*k+2] = point[2];
			k++;
#ifdef USE_GUI_DLL
			if( GetGmProgress() == 1 && (k < 100 || k % 100 == 0 )){
				gmProgressSetPosition(k, progress);
				gm_msg_loop();
			}
#endif
		}
	}
    printf("END\n\n ");


#ifdef USE_GUI_DLL
	if( GetGmProgress() == 1 ){
		gmProgressDeleteForm(progress);
	}
#endif
	if ( fit == NULL ) goto rtn;

	
       	printf("C1-スプライン化 STEP1 START ");
#ifdef USE_GUI_DLL
	if( GetGmProgress() == 1 ){
		progress = gmProgressCreateForm();
        	gmProgressSetFormTitle("C1-スプライン化 STEP1", progress);
		gmProgressSetMinMax(0, grid[1], progress);
		gmProgressShowForm(progress);
	}		
#endif
	// U-方向フィッティング
	for ( j = 0; j <= grid[1]; j += 3){
		ii = 0;
		for ( i = 0; i <= grid[0]; i += 3 ){
			wkpnt[3*ii  ] = ctrlpnt[3*(grid[0]*j+i)];
			wkpnt[3*ii+1] = ctrlpnt[3*(grid[0]*j+i)+1];
			wkpnt[3*ii+2] = ctrlpnt[3*(grid[0]*j+i)+2];
			ii++;
		}
		nn = ii;
		assert(nn==(grid[0]-1)/3+1);

		wkcp = (fit)(nn, wkpnt);
		if ( wkcp != NULL ){
			ii = 0;
			for ( i = 0; i < grid[0]; i += 1 ){
				ctrlpnt[3*(grid[0]*j+i)  ] = wkcp[3*ii];
				ctrlpnt[3*(grid[0]*j+i)+1] = wkcp[3*ii+1];
				ctrlpnt[3*(grid[0]*j+i)+2] = wkcp[3*ii+2];
				ii++;
			}
			(fit)(-1, wkcp);
		}
#ifdef USE_GUI_DLL
		if( GetGmProgress() == 1 && (j < 100 || j % 100  == 0)){
			gmProgressSetPosition(j, progress);
			gm_msg_loop();
		}
#endif
	}	
#ifdef USE_GUI_DLL
	if( GetGmProgress() == 1 ){
		gmProgressDeleteForm(progress);
	}		
#endif
	printf("END\n\n ");


       	printf("C1-スプライン化 STEP2 START ");
#ifdef USE_GUI_DLL
	if( GetGmProgress() == 1 ){
		progress = gmProgressCreateForm();
        	gmProgressSetFormTitle("C1-スプライン化 STEP2", progress);
		gmProgressSetMinMax(0, grid[0], progress);
		gmProgressShowForm(progress);
	}		
#endif
	// V-方向フィッティング
	for ( i = 0; i < grid[0]; i += 1){
		ii = 0;
		for ( j = 0; j <= grid[1]; j += 3 ){
			wkpnt[3*ii  ] = ctrlpnt[3*(grid[0]*j+i)];
			wkpnt[3*ii+1] = ctrlpnt[3*(grid[0]*j+i)+1];
			wkpnt[3*ii+2] = ctrlpnt[3*(grid[0]*j+i)+2];
			ii++;
		}
		nn = ii;
		assert(nn==(grid[1]-1)/3+1);

		wkcp = (fit)(nn, wkpnt);
		if ( wkcp != NULL ){
			ii = 0;
			for ( j = 0; j < grid[1]; j += 1 ){
				ctrlpnt[3*(grid[0]*j+i)  ] = wkcp[3*ii];
				ctrlpnt[3*(grid[0]*j+i)+1] = wkcp[3*ii+1];
				ctrlpnt[3*(grid[0]*j+i)+2] = wkcp[3*ii+2];
				ii++;
			}
			(fit)(-1, wkcp);
		}
#ifdef USE_GUI_DLL
		if( GetGmProgress() == 1 && (i < 100 || i % 100  == 0)){
			gmProgressSetPosition(i, progress);
			gm_msg_loop();
		}
#endif
	}
#ifdef USE_GUI_DLL
	if( GetGmProgress() == 1 ){
		gmProgressDeleteForm(progress);
	}		
#endif
       	printf("END\n\n ");

	if ( wkpnt != NULL ) free( wkpnt);

rtn: ;

#if 0
	printf("end;\n");
#endif
	return ctrlpnt;
}

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

int Polygon2SbezFitC(int* chkgridA, double *ctrlpnt, int upatch, int vpatch, double rectA[][2], 
					 int vtxnum, int trinum, double* vertex, double* uv, unsigned int* tri, double* normal, 
					double* start, double* step, double* (*fit)(int, double*) )
{

	double* cp;
	int chkgrid[3];
	double rect[2][2];

	cp = ctrlpnt;
	if ( chkgridA == NULL ){
		chkgrid[0] = 20;
		chkgrid[1] = 20;
		chkgrid[2] = 20;
	}else{
		chkgrid[0] = chkgridA[0];
		chkgrid[1] = chkgridA[1];
		chkgrid[2] = chkgridA[2];
	}
	if ( rectA == NULL ){
		rect[0][0] = 0.0;
		rect[0][1] = 1.0;
		rect[1][0] = 0.0;
		rect[1][1] = 1.0;
	}else{
		rect[0][0] = rectA[0][0];
		rect[0][1] = rectA[0][1];
		rect[1][0] = rectA[1][0];
		rect[1][1] = rectA[1][1];
	}

//	FILE* fp = DumpSrfStart("sbezer.txt");

	PolygonSurface* ps = new PolygonSurface;
	ps->SetCPArea(ctrlpnt);

	ps->Create(chkgrid, upatch, vpatch, rect,
					vtxnum, trinum, vertex, uv, tri, NULL, NULL, NULL );

/*
	FILE* ffp = fopen("uvline.txt", "w");
	ps->DumpUVLines(ffp, 40,40);
	fclose(ffp);
//	exit(0);
*/


	cp = Polygon2SbezFit( ps, upatch, vpatch, rect, fit );
	if (ps->alcflg_ ) free(cp);

//	DumpSrf(fp, 0, upatch, vpatch, ctrlpnt);
//	DumpSrfEnd(fp);

	delete ps;
	return 0;
}

#ifdef TEST
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

	ReadMeshUV( "G:\\devwork\\Recon\\Recon.940503b\\src\\uvmap_mesh_maxplanck_optmz.dat", &vtxnum, &trinum, &vertex, &uv, &tri );

	int upatch;
	int vpatch;
	int chkgrid[3];
	int grid[3];
	double* ctrlpnt;

	normal = NULL;
	start = NULL;
	step = NULL;

	chkgrid[0] = 20;
	chkgrid[1] = 20;
	chkgrid[2] = 20;

	upatch = 5;
	vpatch = 5;

	grid[0] = 3*upatch+1;
	grid[1] = 3*vpatch+1;
	grid[2] = 1;

	ctrlpnt = (double*)malloc( 3*grid[0]*grid[1]*sizeof(double));
	if ( ctrlpnt == NULL ){
		return NULL;
	}
	FILE* fp = DumpSrfStart("sbezer.txt");

	int stat = Polygon2SbezFitC(chkgrid, ctrlpnt, upatch, vpatch, NULL, 
					 vtxnum, trinum, vertex, uv, tri, normal, start, step, NULL);

	DumpSrf(fp, 0, upatch, vpatch, ctrlpnt);
	free(ctrlpnt);
	DumpSrfEnd(fp);
	return 0;
}
#endif

#ifdef TEST0
int main(int argc, char** argv)
{
	int vtxnum;
	int trinum;
	double* vertex;
	double* uv;
	int* tri;

	ReadMeshUV( "G:\\devwork\\Recon\\Recon.940503b\\src\\uvmap_mesh_maxplanck_optmz.dat", &vtxnum, &trinum, &vertex, &uv, &tri );

	double* cp;

	int upatch;
	int vpatch;
	int chkgrid[3];
	double rect[2][2];

	chkgrid[0] = 20;
	chkgrid[1] = 20;
	chkgrid[2] = 20;


	FILE* fp = DumpSrfStart("sbezer.txt");

	int srfnum[2];
	srfnum[0] = 49;
	srfnum[1] = 49;

	double step[2];
	int kk = 0;

	PolygonSurface* ps = new PolygonSurface;
	
	ps->Create(chkgrid, upatch, vpatch, rect,
					vtxnum, trinum, vertex, uv, tri, NULL, NULL, NULL );

/*
	FILE* ffp = fopen("uvline.txt", "w");
	ps->DumpUVLines(ffp, 40,40);
	fclose(fp);
	exit(0);
*/
#if 0
	upatch = 9;
	vpatch = 9;

	step[0] = 1.0/((double)srfnum[0]);
	step[1] = 1.0/((double)srfnum[1]);
	for ( int ii = 0; ii < srfnum[1]; ii++ ){
		for ( int jj = 0; jj < srfnum[0]; jj++ ){
		
			rect[0][0] = (double)ii     * step[0];
			rect[0][1] = (double)(ii+1) * step[0];
			rect[1][0] = (double)jj     * step[1];
			rect[1][1] = (double)(jj+1) * step[1];
			printf("[%f,%f] x [%f,%f]\n", rect[0][0], rect[1][0], rect[0][1], rect[1][1]);
			cp = Polygon2SbezFit( ps, upatch, vpatch, rect  );
			printf("%dx%d end\n", jj, ii);
			
			DumpSrf(fp, kk, upatch, vpatch, cp);
			kk++;
			free(cp);
		}
	}
#else
	upatch = 35;
	vpatch = 35;

	rect[0][0] = 0.0;
	rect[0][1] = 1.0;
	rect[1][0] = 0.0;
	rect[1][1] = 1.0;
	cp = Polygon2SbezFit( ps, upatch, vpatch, rect );
	DumpSrf(fp, 0, upatch, vpatch, cp);
	free(cp);
#endif
	DumpSrfEnd(fp);

	delete ps;
	return 0;
}
#endif

	
	
