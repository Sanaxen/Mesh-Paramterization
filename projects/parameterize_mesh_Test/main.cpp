#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "parameterize_mesh.h"


int main(int argc, char *argv[])
{
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	char* inputfile = NULL;
	int iter = 1024;
	double stretch_coef = 1.0;
	
	int debug = 0;
	int type = 0;
	double scale = 1.0;
	int mesh_output = 0;
	int uvline_output = 0;
	int output_crv_num = 5;
	int output_crv_div_num = 80;
	int udim = 300;
	int vdim = 300;
	fname_[0] = '\0';

	if ( argc < 2 )
	{
		printf("[-input] filename\n");
		printf("-iter [1 <= N ]\n");
		printf("-stretch  [0 < # <= 1.0]\n");
		printf("-mesh(M,N) 2 <= M, 2 <= N\n");
		printf("-uvline(M,N) 2 <= M, 2 <= N\n");
		printf("-scale [0 < N]\n");
		printf("-type [0:Floater's Mean Value 1:Floater's Shape Preserving 2:Tutte's Barycentric]\n");
		return 0;
	}

	for ( int i = 1; i < argc; i++ )
	{
		if (strcmp(argv[i], "-debug") == 0)
		{
			debug = atoi(argv[i + 1]);
			i++;
			continue;
		}
		if (strcmp(argv[i], "-type") == 0)
		{
			type = atoi(argv[i + 1]);
			i++;
			continue;
		}
		if (strcmp(argv[i], "-scale") == 0)
		{
			scale = atof(argv[i + 1]);
			i++;
			continue;
		}
		if ( strcmp(argv[i], "-input") == 0 )
		{
			inputfile = argv[i+1];
			i++;
			continue;
		}
		if ( strcmp(argv[i], "-iter") == 0 )
		{
			iter = atoi(argv[i+1]);
			i++;
			continue;
		}
		if ( strcmp(argv[i], "-stretch") == 0 )
		{
			stretch_coef = atof(argv[i+1]);
			i++;
			continue;
		}
		if ( strstr(argv[i], "-mesh(")  )
		{
			mesh_output = 1;
			int nn, mm;
			int n = sscanf(argv[i], "-mesh(%d,%d)", &nn, &mm);
			if ( n == 2 )
			{
				udim = nn;
				vdim = mm;
			}
			continue;
		}
		if ( strstr(argv[i], "-mesh")  )
		{
			mesh_output = 1;
			continue;
		}

		if ( strstr(argv[i], "-uvline(") )
		{
			int nn, mm;
			int n = sscanf(argv[i], "-uvline(%d,%d)", &nn, &mm);
			if ( n == 2 )
			{
				output_crv_num = nn;
				output_crv_div_num = mm;
			}
			uvline_output = 1;
			continue;
		}
		if ( strstr(argv[i], "-uvline") )
		{
			uvline_output = 1;
			continue;
		}
		if ( inputfile == NULL )
		{
			inputfile = argv[i];
			mesh_output = 1;
			uvline_output = 1;
			continue;
		}
	}
	printf("-type %d-> [0:Floater's Shape Preserving 1:Floater's Mean Value 2:Tutte's Barycentric]\n", type);

	if ( inputfile == NULL )
	{
		return 0;
	}
	_splitpath(inputfile, drive_, dir_, fname_, ext_ );

	char outfile[256];
	sprintf(outfile, "%s%s%s_out.obj", drive_, dir_, fname_);

	void* env = initialize_parameterize_mesh();

#if 0
	parameterize_mesh_file(argv[1], outfile, iter, stretch_coef);

#else
	MeshVertexList_t mesh;
	int stat = create_mesh_i(inputfile, &mesh, scale);
	printf("stat %d\n", stat);
	if (debug)
	{
		sprintf(outfile, "%s%s%s_in.obj", drive_, dir_, fname_);
		save_mesh(outfile, &mesh);

		sprintf(outfile, "%s%s%s_in.off", drive_, dir_, fname_);
		save_mesh_off(outfile, &mesh);
	}
	try
	{
		stat = parameterize_mesh_io(&mesh, iter, stretch_coef, type, scale);
		printf("stat %d\n", stat);
		if (stat == -1)
		{
			return -1;
		}
	}
	catch (...)
	{
		return -1;
	}

	sprintf(outfile, "%s%s%s_out.obj", drive_, dir_, fname_);
	save_mesh(outfile, &mesh);

	sprintf(outfile, "%s%s%s_boundary.obj", drive_, dir_, fname_);
	save_boundary(outfile, &mesh, -1.0, scale);
#endif

	//sprintf(outfile, "%s_out.ply2", fname_);
	//save_mesh2(outfile, &mesh);

	double param[2]={0,0};
	double p[3];

	printf("CreateSurfaceMesh start\n"); fflush(stdout);
	PolygonSurface_t* ps = CreateSurfaceMesh(&mesh);
	//SurfaceMesh(ps, param, p);
	//printf("%f %f %f\n", p[0], p[1], p[2]);
	printf("CreateSurfaceMesh end\n"); fflush(stdout);

	if ( uvline_output )
	{
		printf("CreateParameterCurvs start\n"); fflush(stdout);
		CreateParameterCurvs(ps, output_crv_num, output_crv_div_num);
		sprintf(outfile, "%s%s%s_parameter.obj", drive_, dir_, fname_);
		save_ParameterCurvs(outfile, ps, -1.0);
		printf("CreateParameterCurvs end\n"); fflush(stdout);
	}
	if ( mesh_output )
	{
		printf("ParameterMesh start\n"); fflush(stdout);
		MeshVertexList_t* mesh2 = (MeshVertexList_t*)new_mesh_p();
		ParameterMesh( ps, udim,  vdim, mesh2);
		sprintf(outfile, "%s%s%s_parameter_mesh.obj", drive_, dir_, fname_);
		save_mesh(outfile, mesh2);

		delete_mesh_p((void*)mesh2);
		printf("ParameterMesh end\n"); fflush(stdout);
	}
	DeleteSurfaceMesh(ps);

	terminate_parameterize_mesh(env);
	return 0;
}