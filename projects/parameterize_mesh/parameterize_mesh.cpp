#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <GI/gi.h>
#include "parameterize_mesh.hpp"

#ifndef M_PI
	#define M_PI			3.1415926535897932
#endif
#define M_TWO_PI			6.2831853071795865

#define VEC3_ADD(d,v,w)		(d)[0]=(v)[0]+(w)[0]; (d)[1]=(v)[1]+(w)[1]; (d)[2]=(v)[2]+(w)[2]
#define VEC3_SUB(d,v,w)		(d)[0]=(v)[0]-(w)[0]; (d)[1]=(v)[1]-(w)[1]; (d)[2]=(v)[2]-(w)[2]
#define VEC3_CROSS(d,v,w)	(d)[0]=(v)[1]*(w)[2] - (v)[2]*(w)[1]; \
							(d)[1]=(v)[2]*(w)[0] - (v)[0]*(w)[2]; \
							(d)[2]=(v)[0]*(w)[1] - (v)[1]*(w)[0]
#define VEC3_NORMALIZE(v)	{ float n=sqrt((v)[0]*(v)[0]+(v)[1]*(v)[1]+(v)[2]*(v)[2]); \
							if(fabs(n)>1e-6) { float m=1.0f/n; (v)[0]*=m; (v)[1]*=m; (v)[2]*=m; } }

unsigned int uiMesh, uiGIM[3], uiList, uiPattern;

void* initialize_parameterize_mesh()
{
	GIcontext pContext;
	pContext = giCreateContext();
	giMakeCurrent(pContext);

	return pContext;
}

void terminate_parameterize_mesh(void* env)
{
	giDestroyContext((GIcontext*)env);
}


void delete_mesh(MeshVertexList_t* mesh)
{
	if ( mesh == NULL )
	{
		return;
	}

	if ( mesh->pIndices ) free( mesh->pIndices);
	if ( mesh->pNormals ) free( mesh->pNormals);
	if ( mesh->pVertices ) free( mesh->pVertices);
	if ( mesh->pUV ) free( mesh->pUV);
	if ( mesh->boundaryLine )
	{
		N_Cylinder* p = (N_Cylinder*)(mesh->boundaryLine);
		delete p;
	}

	mesh->pIndices = NULL;
	mesh->pNormals = NULL;
	mesh->pVertices = NULL;
	mesh->pUV = NULL;
	mesh->boundaryLine = NULL;
}

void save_mesh(char* filename, MeshVertexList_t* mesh)
{
	if ( mesh == NULL )
	{
		return;
	}

	if ( mesh->iNumIndices == 0 )
	{
		return;
	}

	FILE* fp;

	if (filename != NULL )
	{
		fp = fopen( filename, "w");
	}else
	{
		return;
	}

	if ( fp ) fprintf(fp, "mtllib %s.mtl\n", filename);
	for ( int i = 0; i < mesh->iNumVertices; i++ )
	{
		fprintf(fp, "v %f %f %f\n", mesh->pVertices[3*i], mesh->pVertices[3*i+1], mesh->pVertices[3*i+2]);
	}

	if ( mesh->pUV )
	{
		for ( int i = 0; i < mesh->iNumVertices; i++ )
		{
			fprintf(fp, "vt %.10f %.10f\n", mesh->pUV[2*i], mesh->pUV[2*i+1]);
		}
	}
	if ( mesh->pNormals )
	{
		for ( int i = 0; i < mesh->iNumVertices; i++ )
		{
			fprintf(fp, "vn %f %f %f\n", mesh->pNormals[3*i], mesh->pNormals[3*i+1], mesh->pNormals[3*i+2]);
		}
	}

	if ( fp ) fprintf(fp, "usemtl mat\n");
	for ( int i = 0; i < mesh->iNumIndices; i++ )
	{
		if ( mesh->pUV && mesh->pNormals )
		{
			fprintf( fp, "f %d/%d/%d",   mesh->pIndices[3*i]+1, mesh->pIndices[3*i]+1, mesh->pIndices[3*i]+1);
			fprintf( fp, "  %d/%d/%d",   mesh->pIndices[3*i+1]+1, mesh->pIndices[3*i+1]+1, mesh->pIndices[3*i+1]+1);
			fprintf( fp, "  %d/%d/%d\n", mesh->pIndices[3*i+2]+1, mesh->pIndices[3*i+2]+1, mesh->pIndices[3*i+2]+1);
		}else
		if (mesh->pUV == NULL && mesh->pNormals )
		{
			fprintf( fp, "f %d//%d",   mesh->pIndices[3*i]+1, mesh->pIndices[3*i]+1);
			fprintf( fp, "  %d//%d",   mesh->pIndices[3*i+1]+1, mesh->pIndices[3*i+1]+1);
			fprintf( fp, "  %d//%d\n", mesh->pIndices[3*i+2]+1, mesh->pIndices[3*i+2]+1);
		}
		else if ( mesh->pUV && mesh->pNormals == NULL )
		{
			fprintf( fp, "f %d/%d/ ",   mesh->pIndices[3*i]+1, mesh->pIndices[3*i]+1, mesh->pIndices[3*i]+1);
			fprintf( fp, "  %d/%d/ ",   mesh->pIndices[3*i+1]+1, mesh->pIndices[3*i+1]+1, mesh->pIndices[3*i+1]+1);
			fprintf( fp, "  %d/%d/ \n", mesh->pIndices[3*i+2]+1, mesh->pIndices[3*i+2]+1, mesh->pIndices[3*i+2]+1);
		}else
		{
			fprintf( fp, "f %d %d %d\n",   mesh->pIndices[3*i]+1, mesh->pIndices[3*i+1]+1, mesh->pIndices[3*i+2]+1);
		}
	}
	if ( fp ) fclose(fp);
	
	{
		char f[256];
		strcpy(f, filename);
		strcat(f, ".mtl");

		fp = fopen(f, "w");
		if ( fp ) fprintf(fp, "newmtl mat\n"
					"map_Kd checker_1k.bmp\n"
					"Ka 0.25000 0.25000 0.25000\n"
					"Kd 1.00000 1.00000 1.00000\n"
					"Ks 1.00000 1.00000 1.00000\n"
					"Ns 5.00000\n");
		if ( fp ) fclose(fp);
	}
}

void save_mesh2(char* filename, MeshVertexList_t* mesh)
{
	if ( mesh == NULL )
	{
		return;
	}

	if ( mesh->iNumIndices == 0 )
	{
		return;
	}

	FILE* fp;

	if (filename != NULL )
	{
		fp = fopen( filename, "w");
	}else
	{
		return;
	}

	fprintf(fp, "%d\n", mesh->iNumVertices);
	fprintf(fp, "%d\n", mesh->iNumIndices);
	for ( int i = 0; i < mesh->iNumVertices; i++ )
	{
		fprintf(fp, "%f %f %f %f %f\n", mesh->pVertices[3*i], mesh->pVertices[3*i+1], mesh->pVertices[3*i+2],
			mesh->pUV[2*i], mesh->pUV[2*i+1]);
	}

	for ( int i = 0; i < mesh->iNumIndices; i++ )
	{
		fprintf( fp, "%d %d %d\n",   mesh->pIndices[3*i]+1, mesh->pIndices[3*i+1]+1, mesh->pIndices[3*i+2]+1);
	}
	if ( fp ) fclose(fp);
}

void save_mesh_off(char* filename, MeshVertexList_t* mesh)
{
	if (mesh == NULL)
	{
		return;
	}

	if (mesh->iNumIndices == 0)
	{
		return;
	}

	FILE* fp;

	if (filename != NULL)
	{
		fp = fopen(filename, "w");
	}
	else
	{
		return;
	}

	if (fp) fprintf(fp, "OFF\n");
	fprintf(fp, "%d %d 0\n", mesh->iNumVertices, mesh->iNumIndices);
	for (int i = 0; i < mesh->iNumVertices; i++)
	{
		fprintf(fp, "%.10f %.10f %.10f\n", mesh->pVertices[3 * i], mesh->pVertices[3 * i + 1], mesh->pVertices[3 * i + 2]);
	}
	for (int i = 0; i < mesh->iNumIndices; i++)
	{
			fprintf(fp, "3 %d %d %d\n", mesh->pIndices[3 * i], mesh->pIndices[3 * i + 1], mesh->pIndices[3 * i + 2]);
	}
	if (fp) fclose(fp);
}

#define SQR(a) ((a)*(a))
void save_boundary(char* filename, MeshVertexList_t* mesh, double radius, double scale)
{
	N_Cylinder *p = (N_Cylinder*)(mesh->boundaryLine);

	double r = radius;
	if ( radius < 0.0 )
	{
		r = sqrt(SQR(mesh->b_max[0]-mesh->b_min[0])+SQR(mesh->b_max[1]-mesh->b_min[1])+SQR(mesh->b_max[2]-mesh->b_min[2]))/400.0;
	}
	p->put(filename, r / scale);
}

unsigned int create_mesh_io(MeshVertexList_t* mesh)
{
	unsigned int uiMesh;

	/* set attribute arrays */
	giBindAttrib(GI_POSITION_ATTRIB, 0);
	giBindAttrib(GI_PARAM_ATTRIB, 2);
	giAttribPointer(0, 3, GI_FALSE, 0, mesh->pVertices);
	giAttribPointer(1, 3, GI_TRUE, 0, mesh->pNormals);
	giEnableAttribArray(0);
	giEnableAttribArray(1);

	/* create mesh */
	uiMesh = giGenMesh();
	giBindMesh(uiMesh);
	giGetError();
	giIndexedMesh(0, mesh->iNumVertices-1, mesh->iNumIndices*3, mesh->pIndices);
	mesh->pUV = NULL;
	calcBoundaryBox(mesh->b_min, mesh->b_max);

	if(giGetError() != GI_NO_ERROR)
	{
		if ( mesh ) delete_mesh( mesh );
		return 0;
	}
	return uiMesh;

}

unsigned int create_mesh(const char *filename, MeshVertexList_t* mesh, double scale)
{
	unsigned int uiMesh;
	FILE *pFile = fopen(filename, "r");
	GLfloat *pVertices, *pNormals;
	unsigned int *pIndices;
	int i, j, f, iNumVertices, iNumIndices;
	unsigned int i0, i1, i2;
	GLfloat v1[3], v2[3], n[3];
	double val;

	if(!pFile)
	{
		fprintf(stderr, "cannot open file!\n");
		return 0;
	}

	/* read vertex data */
	fscanf(pFile, "%d\n%d\n", &iNumVertices, &iNumIndices);
	iNumVertices *= 3;
	iNumIndices *= 3;
	pVertices = (GLfloat*)malloc(iNumVertices*sizeof(GLfloat));
	pNormals = (GLfloat*)calloc(iNumVertices, sizeof(GLfloat));
	pIndices = (unsigned int*)malloc(iNumIndices*sizeof(unsigned int));
	for(i=0; i<iNumVertices; ++i)
	{
		fscanf(pFile, "%lf\n", &val);
		pVertices[i] = val*scale;
	}
	/* read index data and average face normals */
	for(i=0; i<iNumIndices; i+=3)
	{
		fscanf(pFile, "%d\n", &f);
		for(j=0; j<f; ++j)
			fscanf(pFile, "%d\n", pIndices+i+j);
		i0 = 3 * pIndices[i];
		i1 = 3 * pIndices[i+1];
		i2 = 3 * pIndices[i+2];
		VEC3_SUB(v1, pVertices+i1, pVertices+i0);
		VEC3_SUB(v2, pVertices+i2, pVertices+i0);
		VEC3_CROSS(n, v1, v2);
		VEC3_ADD(pNormals+i0, pNormals+i0, n);
		VEC3_ADD(pNormals+i1, pNormals+i1, n);
		VEC3_ADD(pNormals+i2, pNormals+i2, n);
	}

	/* normalize vertex normals */
	for(i=0; i<iNumVertices; i+=3)
		VEC3_NORMALIZE(pNormals+i);

	/* set attribute arrays */
	giBindAttrib(GI_POSITION_ATTRIB, 0);
	giBindAttrib(GI_PARAM_ATTRIB, 2);
	giAttribPointer(0, 3, GI_FALSE, 0, pVertices);
	giAttribPointer(1, 3, GI_TRUE, 0, pNormals);
	giEnableAttribArray(0);
	giEnableAttribArray(1);

	/* create mesh */
	uiMesh = giGenMesh();
	giBindMesh(uiMesh);
	giGetError();
	giIndexedMesh(0, iNumVertices-1, iNumIndices, pIndices);

	if ( mesh == NULL )
	{
		/* clean up */
		free(pVertices);
		free(pNormals);
		free(pIndices);
	}else
	{
		calcBoundaryBox(mesh->b_min, mesh->b_max);
		mesh->id = uiMesh;
		mesh->boundaryLine = NULL;
		mesh->iNumIndices = iNumIndices/3;
		mesh->iNumVertices = iNumVertices/3;
		mesh->pIndices = pIndices;
		mesh->pNormals = pNormals;
		mesh->pVertices = pVertices;
		mesh->pUV = NULL;
	}

	if(giGetError() != GI_NO_ERROR)
	{
		if ( mesh ) delete_mesh( mesh );
		return 0;
	}
	return uiMesh;

}

unsigned int create_mesh_obj(const char *filename, MeshVertexList_t* mesh, double scale)
{
	unsigned int uiMesh;
	FILE *pFile = fopen(filename, "r");
	GLfloat *pVertices, *pNormals;
	unsigned int *pIndices;
	int i, j, f, iNumVertices, iNumIndices;
	unsigned int i0, i1, i2;
	GLfloat v1[3], v2[3], n[3];
	char buf[256];

	char tmpfile[256];
	FILE* fp;
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	if(!pFile)
	{
		fprintf(stderr, "cannot open file!\n");
		return 0;
	}
	
	_splitpath(filename, drive_, dir_, fname_, ext_ );
	sprintf(tmpfile, "%s%s%s_wrk.obj", drive_, dir_, fname_);

	fp = fopen(tmpfile, "w");

	iNumVertices = 0;
	iNumIndices = 0;
	while( fgets(buf, 256, pFile) != NULL )
	{
		if ( buf[0] == 'v' && !(buf[1] == 'n' || buf[1] == 't')) iNumVertices++;
		if (buf[0] == 'f') iNumIndices++;
	}
	fclose(pFile);

	pFile = fopen(filename, "r");

	printf("vertex %d face %d\n", iNumVertices, iNumIndices);
	/* read vertex data */
	iNumVertices *= 3;
	iNumIndices *= 3;
	pVertices = (GLfloat*)malloc(iNumVertices*sizeof(GLfloat));
	pNormals = (GLfloat*)calloc(iNumVertices, sizeof(GLfloat));
	pIndices = (unsigned int*)malloc(iNumIndices*sizeof(unsigned int));

	for(i=0; i<iNumVertices/3; ++i)
	{
		double x, y, z;
		fgets(buf, 256, pFile);
		while (buf[0] != 'v' || (buf[0] == 'v' && (buf[1] == 'n' || buf[1] == 't')))
		{
			fgets(buf, 256, pFile);
		}
		sscanf(buf, "v %lf %lf %lf", &x, &y, &z);
		pVertices[3 * i] = x*scale;
		pVertices[3 * i + 1] = y*scale;
		pVertices[3 * i + 2] = z*scale;
		//pVertices[3*i] = ((int)(1000000.0*pVertices[3*i]))/1000000.0;
		//pVertices[3*i+1] = ((int)(1000000.0*pVertices[3*i+1]))/1000000.0;
		//pVertices[3*i+2] = ((int)(1000000.0*pVertices[3*i+2]))/1000000.0;
		fprintf(fp, "v %.6f %.6f %.6f\n", pVertices[3*i], pVertices[3*i+1], pVertices[3*i+2]);
	}

	/* read index data and average face normals */
	for(i=0; i<iNumIndices; i+=3)
	{
		fgets(buf, 256, pFile);
		while( buf[0] != 'f' )
		{
			fgets(buf, 256, pFile);
		}
		if (sscanf(buf, "f %d %d %d", pIndices + i, pIndices + i + 1, pIndices + i + 2) != 3)
		{
			int dmy[3];
			if (sscanf(buf, "f %d//%d %d//%d %d//%d", pIndices + i, dmy, pIndices + i + 1, dmy + 1, pIndices + i + 2, dmy + 2) != 6)
			{
				int dmy2[3];
				if (sscanf(buf, "f %d/%d/%d %d/%d/%d %d/%d/%d", pIndices + i, dmy2, dmy, pIndices + i + 1, dmy2 + 1, dmy + 1, pIndices + i + 2, dmy2 + 2, dmy + 2) != 9)
				{
					fprintf(stderr, "no support format!!\n");
					return 0;
				}
			}
		}
		fprintf(fp, "f %d %d %d\n", pIndices[i], pIndices[i+1], pIndices[i+2]);
		pIndices[i] += -1;
		pIndices[i+1] += -1;
		pIndices[i+2] += -1;
		i0 = 3 * pIndices[i];
		i1 = 3 * pIndices[i+1];
		i2 = 3 * pIndices[i+2];
		VEC3_SUB(v1, pVertices+i1, pVertices+i0);
		VEC3_SUB(v2, pVertices+i2, pVertices+i0);
		VEC3_CROSS(n, v1, v2);
		VEC3_ADD(pNormals+i0, pNormals+i0, n);
		VEC3_ADD(pNormals+i1, pNormals+i1, n);
		VEC3_ADD(pNormals+i2, pNormals+i2, n);
	}
	fclose(pFile);
	fclose(fp);

	/* normalize vertex normals */
	for(i=0; i<iNumVertices; i+=3)
		VEC3_NORMALIZE(pNormals+i);

	/* set attribute arrays */
	giBindAttrib(GI_POSITION_ATTRIB, 0);
	giBindAttrib(GI_PARAM_ATTRIB, 2);
	giAttribPointer(0, 3, GI_FALSE, 0, pVertices);
	giAttribPointer(1, 3, GI_TRUE, 0, pNormals);
	giEnableAttribArray(0);
	giEnableAttribArray(1);

	/* create mesh */
	uiMesh = giGenMesh();
	giBindMesh(uiMesh);
	giGetError();
	giIndexedMesh(0, iNumVertices-1, iNumIndices, pIndices);

	if ( mesh == NULL )
	{
		/* clean up */
		free(pVertices);
		free(pNormals);
		free(pIndices);
	}else
	{
		calcBoundaryBox(mesh->b_min, mesh->b_max);
		mesh->id = uiMesh;
		mesh->boundaryLine = NULL;
		mesh->iNumIndices = iNumIndices/3;
		mesh->iNumVertices = iNumVertices/3;
		mesh->pIndices = pIndices;
		mesh->pNormals = pNormals;
		mesh->pVertices = pVertices;
		mesh->pUV = NULL;
	}

	if(giGetError() != GI_NO_ERROR)
	{
		if ( mesh ) delete_mesh( mesh );
		return 0;
	}
	return uiMesh;
}

unsigned int create_mesh_m(const char *filename, MeshVertexList_t* mesh, double scale)
{
	unsigned int uiMesh;
	FILE *pFile = fopen(filename, "r");
	GLfloat *pVertices, *pNormals;
	unsigned int *pIndices;
	int i, j, f, iNumVertices, iNumIndices;
	unsigned int i0, i1, i2;
	GLfloat v1[3], v2[3], n[3];
	char buf[256];
	int id;

	char tmpfile[256];
	FILE* fp;
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	if(!pFile)
	{
		fprintf(stderr, "cannot open file!\n");
		return 0;
	}

	_splitpath(filename, drive_, dir_, fname_, ext_ );
	sprintf(tmpfile, "%s%s%s_wrk.obj", drive_, dir_, fname_);

	fp = fopen(tmpfile, "w");

	iNumVertices = 0;
	iNumIndices = 0;
	while( fgets(buf, 256, pFile) != NULL )
	{
		if ( buf[0] == 'V' ) iNumVertices++;
		if ( buf[0] == 'F' ) iNumIndices++;
	}
	fclose(pFile);
	

	pFile = fopen(filename, "r");

	/* read vertex data */
	iNumVertices *= 3;
	iNumIndices *= 3;
	pVertices = (GLfloat*)malloc(iNumVertices*sizeof(GLfloat));
	pNormals = (GLfloat*)calloc(iNumVertices, sizeof(GLfloat));
	pIndices = (unsigned int*)malloc(iNumIndices*sizeof(unsigned int));

	for(i=0; i<iNumVertices/3; ++i)
	{
		double x, y, z;
		fgets(buf, 256, pFile);
		while( buf[0] != 'V' )
		{
			fgets(buf, 256, pFile);
		}
		sscanf(buf, "Vertex %d %lf %lf %lf", &id, &x, &y, &z);
		pVertices[3*(id-1)] = x*scale;
		pVertices[3*(id-1)+1] = y*scale;
		pVertices[3*(id-1)+2] = z*scale;
	}
	for(i=0; i<iNumVertices/3; ++i)
	{
		fprintf(fp, "v %.6f %.6f %.6f\n", pVertices[3*i], pVertices[3*i+1], pVertices[3*i+2]);
	}

	/* read index data and average face normals */
	for(i=0; i<iNumIndices; i+=3)
	{
		fgets(buf, 256, pFile);
		while( buf[0] != 'F' )
		{
			fgets(buf, 256, pFile);
		}
		sscanf(buf, "Face %d %d %d %d", &id, pIndices+i, pIndices+i+1, pIndices+i+2);
		fprintf(fp, "f %d %d %d\n", pIndices[i], pIndices[i+1], pIndices[i+2]);
		pIndices[i] += -1;
		pIndices[i+1] += -1;
		pIndices[i+2] += -1;

		i0 = 3 * pIndices[i];
		i1 = 3 * pIndices[i+1];
		i2 = 3 * pIndices[i+2];
		VEC3_SUB(v1, pVertices+i1, pVertices+i0);
		VEC3_SUB(v2, pVertices+i2, pVertices+i0);
		VEC3_CROSS(n, v1, v2);
		VEC3_ADD(pNormals+i0, pNormals+i0, n);
		VEC3_ADD(pNormals+i1, pNormals+i1, n);
		VEC3_ADD(pNormals+i2, pNormals+i2, n);
	}
	fclose(pFile);
	fclose(fp);

	/* normalize vertex normals */
	for(i=0; i<iNumVertices; i+=3)
		VEC3_NORMALIZE(pNormals+i);

	/* set attribute arrays */
	giBindAttrib(GI_POSITION_ATTRIB, 0);
	giBindAttrib(GI_PARAM_ATTRIB, 2);
	giAttribPointer(0, 3, GI_FALSE, 0, pVertices);
	giAttribPointer(1, 3, GI_TRUE, 0, pNormals);
	giEnableAttribArray(0);
	giEnableAttribArray(1);

	/* create mesh */
	uiMesh = giGenMesh();
	giBindMesh(uiMesh);
	giGetError();
	giIndexedMesh(0, iNumVertices-1, iNumIndices, pIndices);

	if ( mesh == NULL )
	{
		/* clean up */
		free(pVertices);
		free(pNormals);
		free(pIndices);
	}else
	{
		calcBoundaryBox(mesh->b_min, mesh->b_max);

		mesh->id = uiMesh;
		mesh->boundaryLine = NULL;
		mesh->iNumIndices = iNumIndices/3;
		mesh->iNumVertices = iNumVertices/3;
		mesh->pIndices = pIndices;
		mesh->pNormals = pNormals;
		mesh->pVertices = pVertices;
		mesh->pUV = NULL;
	}

	if(giGetError() != GI_NO_ERROR)
	{
		if ( mesh ) delete_mesh( mesh );
		return 0;
	}
	return uiMesh;
}

unsigned int create_mesh_off(const char *filename, MeshVertexList_t* mesh, double scale)
{
	unsigned int uiMesh;
	FILE *pFile = fopen(filename, "r");
	GLfloat *pVertices, *pNormals;
	unsigned int *pIndices;
	int i, j, f, iNumVertices, iNumIndices;
	unsigned int i0, i1, i2;
	GLfloat v1[3], v2[3], n[3];
	char buf[256];
	int id;

	char tmpfile[256];
	FILE* fp;
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	if (!pFile)
	{
		fprintf(stderr, "cannot open file!\n");
		return 0;
	}

	_splitpath(filename, drive_, dir_, fname_, ext_);
	sprintf(tmpfile, "%s%s%s_wrk.obj", drive_, dir_, fname_);

	fp = fopen(tmpfile, "w");

	fgets(buf, 256, pFile);	// skipp "OFF"
	fgets(buf, 256, pFile);

	int dmy;
	iNumVertices = 0;
	iNumIndices = 0;
	sscanf(buf, "%d %d %d", &iNumVertices, &iNumIndices, &dmy);
	printf("%d %d\n", iNumVertices, iNumIndices);

	/* read vertex data */
	iNumVertices *= 3;
	iNumIndices *= 3;
	pVertices = (GLfloat*)malloc(iNumVertices * sizeof(GLfloat));
	pNormals = (GLfloat*)calloc(iNumVertices, sizeof(GLfloat));
	pIndices = (unsigned int*)malloc(iNumIndices * sizeof(unsigned int));

	for (i = 0; i<iNumVertices / 3; ++i)
	{
		double x, y, z;
		fgets(buf, 256, pFile);
		sscanf(buf, "%lf %lf %lf", &x, &y, &z);
		pVertices[3 * i] = x*scale;
		pVertices[3 * i + 1] = y*scale;
		pVertices[3 * i + 2] = z*scale;
	}
	for (i = 0; i<iNumVertices / 3; ++i)
	{
		fprintf(fp, "v %.6f %.6f %.6f\n", pVertices[3 * i], pVertices[3 * i + 1], pVertices[3 * i + 2]);
	}

	/* read index data and average face normals */
	for (i = 0; i<iNumIndices; i += 3)
	{
		int nf = 0;
		fgets(buf, 256, pFile);
		sscanf(buf, "%d %d %d %d", &nf, pIndices + i, pIndices + i + 1, pIndices + i + 2);
		fprintf(fp, "f %d %d %d\n", pIndices[i], pIndices[i + 1], pIndices[i + 2]);
		//pIndices[i] += -1;
		//pIndices[i + 1] += -1;
		//pIndices[i + 2] += -1;

		i0 = 3 * pIndices[i];
		i1 = 3 * pIndices[i + 1];
		i2 = 3 * pIndices[i + 2];
		VEC3_SUB(v1, pVertices + i1, pVertices + i0);
		VEC3_SUB(v2, pVertices + i2, pVertices + i0);
		VEC3_CROSS(n, v1, v2);
		VEC3_ADD(pNormals + i0, pNormals + i0, n);
		VEC3_ADD(pNormals + i1, pNormals + i1, n);
		VEC3_ADD(pNormals + i2, pNormals + i2, n);
	}
	fclose(pFile);
	fclose(fp);
	printf("OK!!!\n");

	/* normalize vertex normals */
	for (i = 0; i<iNumVertices; i += 3)
		VEC3_NORMALIZE(pNormals + i);

	/* set attribute arrays */
	giBindAttrib(GI_POSITION_ATTRIB, 0);
	giBindAttrib(GI_PARAM_ATTRIB, 2);
	giAttribPointer(0, 3, GI_FALSE, 0, pVertices);
	giAttribPointer(1, 3, GI_TRUE, 0, pNormals);
	giEnableAttribArray(0);
	giEnableAttribArray(1);

	/* create mesh */
	uiMesh = giGenMesh();
	giBindMesh(uiMesh);
	giGetError();
	giIndexedMesh(0, iNumVertices - 1, iNumIndices, pIndices);

	if (mesh == NULL)
	{
		/* clean up */
		free(pVertices);
		free(pNormals);
		free(pIndices);
	}
	else
	{
		calcBoundaryBox(mesh->b_min, mesh->b_max);

		mesh->id = uiMesh;
		mesh->boundaryLine = NULL;
		mesh->iNumIndices = iNumIndices / 3;
		mesh->iNumVertices = iNumVertices / 3;
		mesh->pIndices = pIndices;
		mesh->pNormals = pNormals;
		mesh->pVertices = pVertices;
		mesh->pUV = NULL;
	}

	if (giGetError() != GI_NO_ERROR)
	{
		if (mesh) delete_mesh(mesh);
		return 0;
	}
	return uiMesh;
}

/* cut mesh */
int cut_mesh(unsigned int mesh)
{
	GIboolean bCut;

	/* cut mesh */ 
	giBindMesh(mesh);
	giCutterParameteri(GI_CUTTER, GI_INITIAL_GIM);
	giCut();
	//putCut("cutline.obj", 0.025);

	/* check success */
	giGetMeshbv(GI_HAS_CUT, &bCut);
	if(!bCut)
		return -1;
	return 0;
}

/* parameterize mesh */
int parameterize_mesh(unsigned int mesh, int res, GIfloat stretch_weight, int type)
{
	GIboolean bParams;

	/* configure parameterization and parameterize mesh */ 
	giBindMesh(mesh);
	//giParameterizerParameteri(GI_PARAMETERIZER, GI_STRETCH_MINIMIZING);		/* default */
	giParameterizerParameteri(GI_PARAMETERIZER, GI_GIM);

	switch (type)
	{
	case 1:
		giParameterizerParameteri(GI_INITIAL_PARAMETERIZATION, GI_MEAN_VALUE);		/* default */
		break;
	case 0:
		giParameterizerParameteri(GI_INITIAL_PARAMETERIZATION, GI_SHAPE_PRESERVING);
		break;
	case 2:
		giParameterizerParameteri(GI_INITIAL_PARAMETERIZATION, GI_TUTTE_BARYCENTRIC);
		break;
	case 3:
		giParameterizerParameteri(GI_INITIAL_PARAMETERIZATION, GI_DISCRETE_HARMONIC);
		break;
	case 4:
		giParameterizerParameteri(GI_INITIAL_PARAMETERIZATION, GI_DISCRETE_AUTHALIC);
		break;
	}
	
	//giParameterizerParameteri(GI_INITIAL_PARAMETERIZATION, GI_DISCRETE_AUTHALIC);
	//giParameterizerParameteri(GI_INITIAL_PARAMETERIZATION, GI_DISCRETE_HARMONIC);

	//giParameterizerParameteri(GI_STRETCH_METRIC, GI_MAX_GEOMETRIC_STRETCH);
	//giParameterizerParameteri(GI_STRETCH_METRIC, GI_RMS_GEOMETRIC_STRETCH);	/* default */
	//giParameterizerParameteri(GI_STRETCH_METRIC, GI_COMBINED_STRETCH);

	
	//giParameterizerParameterf(GI_AREA_WEIGHT, 0.5);
	//giParameterizerParameterf(GI_CONFORMAL_WEIGHT, 0.5);
	giParameterizerParameterf(GI_STRETCH_WEIGHT, stretch_weight);
	giParameterizerParameteri(GI_PARAM_RESOLUTION, res);
	giParameterize();

	/* check success */
	giGetMeshbv(GI_HAS_PARAMS, &bParams);
	if(!bParams)
		return -1;
	return 0;
}

int create_mesh_i(char* filename, MeshVertexList_t* mesh, double scale)
{
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	printf("topology check ...\n");

	_splitpath(filename, drive_, dir_, fname_, ext_ );

	if (strcmp(ext_, ".off") == 0)
	{
		printf("load .off\n");
		if (!(uiMesh = create_mesh_off(filename, mesh, scale)))
		{
			fprintf(stderr, "topology error!\n");
			return -1;
		}
	}
	else
		if ( strcmp(ext_, ".ply2") == 0 )
	{
		if(!(uiMesh = create_mesh(filename, mesh, scale)))
		{
			fprintf(stderr, "topology error!\n");
			return -1;
		}
	}else
	if ( strcmp(ext_, ".obj") == 0 )
	{
		if(!(uiMesh = create_mesh_obj(filename, mesh, scale)))
		{
			fprintf(stderr, "Error of the mesh topology!\n");
			return -1;
		}
	}else
	if ( strcmp(ext_, ".m") == 0 )
	{
		if(!(uiMesh = create_mesh_m(filename, mesh, scale)))
		{
			fprintf(stderr, "Error of the mesh topology!\n");
			return -1;
		}
	}else
	{
		fprintf(stderr, "Error of the mesh topology!\n");
		return -1;
	}
	giDeleteMesh(uiMesh);

	return 0;
}

int parameterize_mesh_file(char* infile, char* outfile, int Iter, double stretch_coef, int type )
{
	MeshVertexList_t mesh;

	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	printf("topology check...\n");

	_splitpath(infile, drive_, dir_, fname_, ext_ );

	if ( strcmp(ext_, ".ply2") == 0 )
	{
		if(!(uiMesh = create_mesh(infile, &mesh, 1.0)))
		{
			fprintf(stderr, "topology error!\n");
			return -1;
		}
	}else
	if ( strcmp(ext_, ".obj") == 0 )
	{
		if(!(uiMesh = create_mesh_obj(infile, &mesh, 1.0)))
		{
			fprintf(stderr, "topology error!\n");
			return -1;
		}
	}else
	if ( strcmp(ext_, ".m") == 0 )
	{
		if(!(uiMesh = create_mesh_m(infile, &mesh, 1.0)))
		{
			fprintf(stderr, "topology error!\n");
			return -1;
		}
	}else
	{
		fprintf(stderr, "topology error!\n");
		return -1;
	}
	GIdouble bmin[3], bmax[3];
	calcBoundaryBox2(bmin, bmax);
	delete_mesh(&mesh);

	printf("cutting...\n");
	if(cut_mesh(uiMesh) == -1)
	{
		fprintf(stderr, "cutting error!\n");
		return -1;
	}
	printf("parameterizing...\n");
	if(parameterize_mesh(uiMesh, Iter, stretch_coef, type) == -1)
	{
		fprintf(stderr, "parameterization error!\n");
		{
			char outfile[256];

			sprintf(outfile, "%s_err.obj", fname_);
			putMesh(outfile, NULL, NULL, NULL, NULL, NULL, NULL, 1.0);
		}
		return -1;
	}

	giGLAttribRenderParameteri(2, GI_TEXTURE_COORD_DOMAIN, GI_UNIT_SQUARE);

	/* draw mesh into display list */
	giGLAttribRenderParameteri(0, GI_GL_RENDER_SEMANTIC, GI_GL_VERTEX);
	giGLAttribRenderParameteri(1, GI_GL_RENDER_SEMANTIC, GI_GL_NORMAL);
	giGLAttribRenderParameteri(2, GI_GL_RENDER_SEMANTIC, GI_GL_TEXTURE_COORD);
	giGLAttribRenderParameteri(2, GI_GL_RENDER_CHANNEL, 0);

	putMesh(outfile, NULL, NULL, NULL, NULL, NULL, NULL, 1.0);

	N_Cylinder ncyl;
	putCut(0.025, ncyl, 1.0);
	{
		char outfile[256];

		sprintf(outfile, "%s_boundary.obj", fname_);
		ncyl.put(outfile, sqrt(SQR(bmax[0]-bmin[0])+SQR(bmax[1]-bmin[1])+SQR(bmax[2]-bmin[2]))/400.0);
	}

	/* clean up */
	giDeleteMesh(uiMesh);
	return 0;
}

int parameterize_mesh_io(MeshVertexList_t* mesh, int Iter, double stretch_coef, int type, double scale)
{
	printf("topology check...\n");
	if(!(uiMesh = create_mesh_io(mesh)))
	{
		fprintf(stderr, "topology error!\n");
		return -1;
	}

	printf("cutting...\n");
	if(cut_mesh(uiMesh) == -1)
	{
		fprintf(stderr, "cutting error!\n");
		return -1;
	}

	printf("parameterizing...\n");
	if(parameterize_mesh(uiMesh, Iter, stretch_coef, type) == -1)
	{
		fprintf(stderr, "parameterization error!\n");
		{
			char outfile[256];

			sprintf(outfile, "parameterize_err.obj");
			putMesh(outfile, NULL, NULL, NULL, NULL, NULL, NULL, scale);
		}
		return -1;
	}

	giGLAttribRenderParameteri(2, GI_TEXTURE_COORD_DOMAIN, GI_UNIT_SQUARE);

	/* draw mesh into display list */
	giGLAttribRenderParameteri(0, GI_GL_RENDER_SEMANTIC, GI_GL_VERTEX);
	giGLAttribRenderParameteri(1, GI_GL_RENDER_SEMANTIC, GI_GL_NORMAL);
	giGLAttribRenderParameteri(2, GI_GL_RENDER_SEMANTIC, GI_GL_TEXTURE_COORD);
	giGLAttribRenderParameteri(2, GI_GL_RENDER_CHANNEL, 0);

	delete_mesh(mesh);
	putMesh(NULL, &(mesh->iNumIndices), &(mesh->iNumVertices), &(mesh->pIndices), &(mesh->pVertices), &(mesh->pUV), &(mesh->pNormals), scale);

	N_Cylinder* ncyl = new N_Cylinder;
	const double r =
		(mesh->b_max[0] - mesh->b_min[0])*(mesh->b_max[0] - mesh->b_min[0]) +
		(mesh->b_max[1] - mesh->b_min[1])*(mesh->b_max[1] - mesh->b_min[1]) +
		(mesh->b_max[2] - mesh->b_min[2])*(mesh->b_max[2] - mesh->b_min[2]);

	putCut(sqrt(r)*0.001, *ncyl, scale);
	mesh->boundaryLine = (void*)ncyl;

	/* clean up */
	giDeleteMesh(uiMesh);
	return 0;
}






void ParameterCurvs(PolygonSurface_t* ps, int n,  int polydivn, std::vector<N_Cylinder>& lines );
PolygonSurface_t* CreateSurfaceMesh(MeshVertexList_t* mesh)
{
	int vtxnum;
	int trinum;
	double* vertex;
	double* uv;
	double* normal;
	int* tri;
	int chkgrid[3];
	double rect[2][2];
	chkgrid[0] = 20;
	chkgrid[1] = 20;
	chkgrid[2] = 20;
	rect[0][0] = 0.0;
	rect[0][1] = 1.0;
	rect[1][0] = 0.0;
	rect[1][1] = 1.0;

	PolygonSurface_t* p = new PolygonSurface_t;
	p->mesh_ = mesh;
	p->ps = new PolygonSurface;
	p->parameterCurves = NULL;

	uv = (double*)malloc(3*mesh->iNumVertices*sizeof(double));
	for ( int i = 0; i < mesh->iNumVertices; i++ )
	{
		uv[3*i] = mesh->pUV[2*i];
		uv[3*i+1] = mesh->pUV[2*i+1];
		uv[3*i+2] = 0.0;
	}
	PolygonSurface* ps = (PolygonSurface*)(p->ps);
	ps->Create(chkgrid, 0, 0, rect,
		mesh->iNumVertices, mesh->iNumIndices, mesh->pVertices, uv, mesh->pIndices, NULL, NULL, NULL );

	return p;
}

void DeleteSurfaceMesh(PolygonSurface_t* p)
{
	PolygonSurface* ps = (PolygonSurface*)(p->ps);
	free( ps->uv_);
	std::vector<N_Cylinder>* cv = (std::vector<N_Cylinder>*)p->parameterCurves;
	delete cv;

	delete ps;
	delete p;
}

void SurfaceMesh(PolygonSurface_t* ps, double* uv, double point[3])
{
	PolygonSurface* p = (PolygonSurface*)(ps->ps);
	p->Evalue(uv, point);
}


void ParameterCurv(PolygonSurface_t* ps, double param, int u_or_v, int polydivn, N_Cylinder& lines )
{
	double s = 1.0/(double)polydivn;

	double p1[3], p2[3];
	double uv[2];

	uv[0] = uv[1] = 0.0;
	if ( u_or_v == 0 ) uv[0] = param;
	else uv[1] = param;

	SurfaceMesh(ps, uv, p1);
	for ( int i = 1; i <= polydivn; i++ )
	{
		if ( u_or_v == 0 ) uv[1] = s*i;
		else uv[0] = s*i;
		SurfaceMesh(ps, uv, p2);

		Cylinder cyl(p1, p2, -1.0);
		p1[0] = p2[0];
		p1[1] = p2[1];
		p1[2] = p2[2];
		lines.Add(cyl);
	}
}

void ParameterCurvs(PolygonSurface_t* ps, int n,  int polydivn, std::vector<N_Cylinder>* lines )
{

	int uvflag = 0;

	double ds = 1.0/(double)n;

	for ( int i = 0; i <= n; i++ )
	{
		N_Cylinder polyline;
		double param = i*ds;
		ParameterCurv( ps, param, uvflag, polydivn, polyline );
		lines->push_back( polyline );
	}
	uvflag = 1;
	for ( int i = 0; i <= n; i++ )
	{
		N_Cylinder polyline;
		double param = i*ds;
		ParameterCurv( ps, param, uvflag, polydivn, polyline );
		lines->push_back( polyline );
	}
}


void CreateParameterCurvs(PolygonSurface_t* ps, int n,  int polydivn)
{
	ps->parameterCurves = (void*)(new std::vector<N_Cylinder>);
	std::vector<N_Cylinder>* cv = (std::vector<N_Cylinder>*)ps->parameterCurves;
	ParameterCurvs(ps, n, polydivn, cv);

}

void save_ParameterCurvs(char* filename, PolygonSurface_t* ps, double radius)
{
	double r = radius;
	if ( radius < 0.0 )
	{
		r = sqrt(SQR(ps->mesh_->b_max[0]-ps->mesh_->b_min[0])+SQR(ps->mesh_->b_max[1]-ps->mesh_->b_min[1])+SQR(ps->mesh_->b_max[2]-ps->mesh_->b_min[2]))/650.0;
	}

	std::vector<N_Cylinder>* cv = (std::vector<N_Cylinder>*)ps->parameterCurves;

	FILE* fp = fopen(filename, "w");
	int vcont = 1;
	for ( int i = 0; i < cv->size(); i++ )
	{
		(*cv)[i].put(fp, r, vcont);
	}
	fclose(fp);
}


void ParameterMesh(PolygonSurface_t* ps, int m,  int n, MeshVertexList_t* mesh)
{

	int uvflag = 0;

	double ds = 1.0/(double)m;
	double dt = 1.0/(double)n;

	mesh->iNumVertices = (n+1)*(m+1);
	mesh->iNumIndices = 2*(n)*(m);
	mesh->pVertices = (double*)malloc(3*mesh->iNumVertices*sizeof(double));
	mesh->pUV = (double*)malloc(2*mesh->iNumVertices*sizeof(double));
	mesh->pIndices = (unsigned int*)malloc(3*mesh->iNumIndices*sizeof(unsigned int));

	int k = 0;
#if 0
	double param[2];
	for ( int j = 0; j < n+1; j++ )
	{
		for ( int i = 0; i < m+1; i++ )
		{
			param[0] = i*ds;
			param[1] = j*dt;
			SurfaceMesh(ps, param, &(mesh->pVertices[3*k]));
			mesh->pUV[2*k] = param[0];
			mesh->pUV[2*k+1] = param[1];
			k++;
		}
	}

	k = 0;
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
		{
			mesh->pIndices[3 * k] = (m + 1)*j + i;
			mesh->pIndices[3 * k + 1] = (m + 1)*j + (i + 1);
			mesh->pIndices[3 * k + 2] = (m + 1)*(j + 1) + (i + 1);
			k++;
			mesh->pIndices[3 * k] = (m + 1)*j + i;
			mesh->pIndices[3 * k + 1] = (m + 1)*(j + 1) + (i + 1);
			mesh->pIndices[3 * k + 2] = (m + 1)*(j + 1) + i;
			k++;
		}
	}

#else
#pragma omp parallel
	{
#pragma omp for
		for (int k = 0; k < (n + 1)*(m + 1); k++)
		{
			int j = k / (m + 1);
			int i = k % (m + 1);

			double param[2];
			param[0] = i*ds;
			param[1] = j*dt;
			SurfaceMesh(ps, param, &(mesh->pVertices[3 * k]));
			mesh->pUV[2 * k] = param[0];
			mesh->pUV[2 * k + 1] = param[1];
		}

#pragma omp for
		for (int k = 0; k < n*m; k += 1)
		{
			int j = k / m;
			int i = k % m;

			int kk = 2*k;
			mesh->pIndices[3 * kk] = (m + 1)*j + i;
			mesh->pIndices[3 * kk + 1] = (m + 1)*j + (i + 1);
			mesh->pIndices[3 * kk + 2] = (m + 1)*(j + 1) + (i + 1);

			kk++;
			mesh->pIndices[3 * kk] = (m + 1)*j + i;
			mesh->pIndices[3 * kk + 1] = (m + 1)*(j + 1) + (i + 1);
			mesh->pIndices[3 * kk + 2] = (m + 1)*(j + 1) + i;
		}
	}
#endif
}


void* new_mesh_p(void)
{
	MeshVertexList_t*p = (MeshVertexList_t*)malloc(sizeof(MeshVertexList_t));
	memset(p, '\0', sizeof(MeshVertexList_t));
	return (void*)p;
}
void delete_mesh_p(void* mesh)
{
	MeshVertexList_t* m = (MeshVertexList_t*)mesh;
	delete_mesh(m);
}
int create_mesh_i_p(char* filename, void* mesh)
{
	MeshVertexList_t* m = (MeshVertexList_t*)mesh;
	return create_mesh_i(filename, m, 1.0);
}
void save_mesh_p(char* filename, void* mesh)
{
	MeshVertexList_t* m = (MeshVertexList_t*)mesh;
	save_mesh(filename, m);
}
void save_boundary_p(char* filename, void* mesh, double radius)
{
	MeshVertexList_t* m = (MeshVertexList_t*)mesh;
	save_boundary(filename, m, radius, 1.0);
}
int parameterize_mesh_io_p(void* mesh, int Iter, double stretch_coef )
{
	MeshVertexList_t* m = (MeshVertexList_t*)mesh;
	return parameterize_mesh_io(m, Iter, stretch_coef, 0, 1.0);
}


void* CreateSurfaceMesh_p(void* mesh)
{
	MeshVertexList_t* m = (MeshVertexList_t*)mesh;

	return (void*)CreateSurfaceMesh(m);
}

void DeleteSurfaceMesh_p(void* p)
{
	PolygonSurface_t* ps = (PolygonSurface_t*)p;
	DeleteSurfaceMesh(ps);
}

 void SurfaceMesh_p(void* p, double* uv, double point[3])
 {
	PolygonSurface_t* ps = (PolygonSurface_t*)p;
	SurfaceMesh(ps, uv, point);
 }

void CreateParameterCurvs_p(void* p, int n,  int polydivn)
{
	PolygonSurface_t* ps = (PolygonSurface_t*)p;
	CreateParameterCurvs(ps, n, polydivn);
}

void save_ParameterCurvs_p(char* filename, void* p, double radius)
{
	PolygonSurface_t* ps = (PolygonSurface_t*)p;
	save_ParameterCurvs(filename, ps, radius);
}

void ParameterMesh_p(void* ps, int im,  int in, void* mesh)
{
	PolygonSurface_t* p = (PolygonSurface_t*)ps;
	MeshVertexList_t* m = (MeshVertexList_t*)mesh;
	ParameterMesh(p, im, in, m);
}
