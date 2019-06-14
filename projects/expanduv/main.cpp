#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Windows.h>


#ifndef M_PI
	#define M_PI			3.1415926535897932
#endif
#define M_TWO_PI			6.2831853071795865

#define VEC3_ADD(d,v,w)		(d)[0]=(v)[0]+(w)[0]; (d)[1]=(v)[1]+(w)[1]; (d)[2]=(v)[2]+(w)[2]
#define VEC3_SUB(d,v,w)		(d)[0]=(v)[0]-(w)[0]; (d)[1]=(v)[1]-(w)[1]; (d)[2]=(v)[2]-(w)[2]
#define VEC3_CROSS(d,v,w)	(d)[0]=(v)[1]*(w)[2] - (v)[2]*(w)[1]; \
							(d)[1]=(v)[2]*(w)[0] - (v)[0]*(w)[2]; \
							(d)[2]=(v)[0]*(w)[1] - (v)[1]*(w)[0]
#define VEC3_NORMALIZE(v)	{ double n=sqrt((v)[0]*(v)[0]+(v)[1]*(v)[1]+(v)[2]*(v)[2]); \
							if(fabs(n)>1e-16) { double m=1.0/n; (v)[0]*=m; (v)[1]*=m; (v)[2]*=m; } }


#define SQR(a)	((a)*(a))

#define DOTPRODUCT(a)		((a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2])
#define DOTPRODUCT2(a,b)	((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

struct Image {    
	unsigned long sizeX;  //横    
	unsigned long sizeY; //縦    
	char *data;  //画像データ格納
};
typedef struct Image Image;


int ImageLoad(char *filename, Image *image) 
{    
	FILE *file;    
	unsigned long size;               // イメージのバイトサイズ    
	unsigned long i;                    // カウンタ    
	unsigned short int planes;        //デバイス面数    
	unsigned short int bpp;            // ピクセル数    
	char temp;                         // カラー関係作業用    
	//ファイルオープン    
	if ((file = fopen(filename, "rb"))==NULL)    
	{            
		printf("File Not Found : %s\n",filename);            
		return 0;    
	}    
	//ビットマップの幅データ部分へ移動    
	fseek(file, 18, SEEK_CUR);    
	//横幅を読み込む    
	if ((i = fread(&image->sizeX, 4, 1, file)) != 1) 
	{       
		printf("Error reading width from %s.\n", filename);        
		return 0;   
	}   
	//縦幅を読み込む    
	if ((i = fread(&image->sizeY, 4, 1, file)) != 1) 
	{        
		printf("Error reading height from %s.\n", filename);        
		return 0;    
	}     
	//画像サイズの計算     
	size = image->sizeX * image->sizeY * 3;     
	//プレーン数を読み込む    
	if ((fread(&planes, 2, 1, file)) != 1) 
	{  
		//bmpは「1」になる        
		printf("Error reading planes from %s.\n", filename);        
		return 0;    
	}    
	if (planes != 1) 
	{         
		printf("Planes from %s is not 1: %u\n", filename, planes);         
		return 0;   
	}    
	//ピクセル値を読み込む    
	if ((i = fread(&bpp, 2, 1, file)) != 1) 
	{        
		printf("Error reading bpp from %s.\n", filename);        
		return 0;
	}    
	if (bpp != 24) 
	{   
		//24bppでなければエラー         
		printf("Bpp from %s is not 24: %u\n", filename, bpp);        
		return 0;    
	}    
	//24ビット飛ばして、カラーデータ(RGB)部分へ    
	fseek(file, 24, SEEK_CUR);    
	//データ読み込み    
	image->data = (char *) malloc(size);    
	if (image->data == NULL) 
	{        
		printf("Error allocating memory for color-corrected image data");        
		return 0;          
	}    
	if ((i = fread(image->data, size, 1, file)) != 1) 
	{        
		printf("Error reading image data from %s.\n", filename);        
		return 0;   
	}   
	for (i=0;i<size;i+=3) 
	{          
		// reverse all of the colors. (bgr -> rgb)        
		temp = image->data[i];        
		image->data[i] = image->data[i+2];        
		image->data[i+2] = temp;    
	}    
	return 1;
}

	
	
class Color
{
public:
	int rgb[3];
};

class MeshV
{
public:
	unsigned int *pIndices;
	int iNumVertices;
	int iNumIndices;
	double *pVertices;
	double *pNormals;
	double* pUV;
	Color* pColors;

	MeshV()
	{
		pIndices = NULL;
		pVertices = NULL;
		pNormals = NULL;
		pColors = NULL;
		pUV = NULL;
		iNumVertices = 0;
		iNumIndices = 0;
	}

	~MeshV()
	{
		if ( pIndices ) delete [] pIndices;
		if ( pVertices ) delete [] pVertices;
		if ( pNormals ) delete [] pNormals;
		if ( pColors ) delete [] pColors;
		if ( pUV ) delete [] pUV;
		pIndices = NULL;
		pVertices = NULL;
		pNormals = NULL;
		pColors = NULL;
		iNumVertices = 0;
		iNumIndices = 0;
	}

	double* VertexCoord(int id)
	{
		return &(pVertices[3*id]);
	}
	double* UVCoord(int id)
	{
		return &(pUV[2*id]);
	}
	double min[3];
	double max[3];
};

int create_mesh_obj(const char *filename, MeshV& mesh)
{
	int i, j, f;
	unsigned int i0, i1, i2;
	double v1[3], v2[3], n[3];
	char buf[256];

	FILE *pFile = fopen(filename, "r");
	if(!pFile)
	{
		fprintf(stderr, "cannot open file!\n");
		return -1;
	}
	
	while( fgets(buf, 256, pFile) != NULL )
	{
		if ( buf[0] == 'v' && (buf[1] == 't' || buf[1] == 'n') ) continue;
		if ( buf[0] == 'v' ) mesh.iNumVertices++;
		if ( buf[0] == 'f' ) mesh.iNumIndices++;
	}
	fclose(pFile);

	pFile = fopen(filename, "r");

	/* read vertex data */
	mesh.iNumVertices *= 3;
	mesh.iNumIndices *= 3;
	mesh.pVertices = new double[mesh.iNumVertices];
	mesh.pUV = new double[2*mesh.iNumVertices/3];
	mesh.pNormals = new double [mesh.iNumVertices];
	mesh.pIndices = new unsigned int[mesh.iNumIndices];
	mesh.pColors = new Color[mesh.iNumVertices/3];

	bool use_color = true;
	for(i=0; i<mesh.iNumVertices/3; ++i)
	{
		double x, y, z;
		int rgb[3];
		float rgbf[3];
		fgets(buf, 256, pFile);
		while( buf[0] != 'v' || (buf[0] == 'v' && (buf[1] == 't' || buf[1] == 'n')))
		{
			fgets(buf, 256, pFile);
		}
		if ( sscanf(buf, "v %lf %lf %lf %d %d %d", &x, &y, &z, rgb, rgb+1, rgb+2) != 6 )
		{
			if ( sscanf(buf, "v %lf %lf %lf %f %f %f", &x, &y, &z, rgbf, rgbf+1, rgbf+2) != 6 )
			{
				sscanf(buf, "v %lf %lf %lf", &x, &y, &z);
				use_color = false;
			}else
			{
				rgb[0] = (int)(255.0*rgbf[0]);
				rgb[1] = (int)(255.0*rgbf[1]);
				rgb[2] = (int)(255.0*rgbf[2]);
			}
		}
		mesh.pVertices[3*i] = x;
		mesh.pVertices[3*i+1] = y;
		mesh.pVertices[3*i+2] = z;
		mesh.pNormals[3*i] = 0.0;
		mesh.pNormals[3*i+1] = 0.0;
		mesh.pNormals[3*i+2] = 0.0;
		if ( use_color )
		{
			mesh.pColors[i].rgb[0] = rgb[0];
			mesh.pColors[i].rgb[1] = rgb[1];
			mesh.pColors[i].rgb[2] = rgb[2];
		}
	}

	if ( !use_color )
	{
		delete mesh.pColors;
		mesh.pColors = NULL;
	}

	for(i=0; i<mesh.iNumVertices/3; ++i)
	{
		double u, v;
		fgets(buf, 256, pFile);
		while( buf[0] != 'v' || (buf[0] == 'v' && ( buf[1] == 'n')))
		{
			fgets(buf, 256, pFile);
		}
		sscanf(buf, "vt %lf %lf", &u, &v);
		mesh.pUV[2*i] = u;
		mesh.pUV[2*i+1] = v;
	}

	for(i=0; i<mesh.iNumVertices/3; ++i)
	{
		double x, y, z;
		fgets(buf, 256, pFile);
		while( buf[0] != 'v')
		{
			fgets(buf, 256, pFile);
		}
		sscanf(buf, "vn %lf %lf %lf", &x, &y, &z);
		mesh.pNormals[3*i] = x;
		mesh.pNormals[3*i+1] = y;
		mesh.pNormals[3*i+2] = z;
	}

	/* read index data and average face normals */
	for(i=0; i<mesh.iNumIndices; i+=3)
	{
		fgets(buf, 256, pFile);
		while( buf[0] != 'f' )
		{
			fgets(buf, 256, pFile);
		}
		bool s = (sscanf(buf, "f %d %d %d", mesh.pIndices+i, mesh.pIndices+i+1, mesh.pIndices+i+2) == 3 );
		if ( !s)
		{
			int nrm[3];
			s = (sscanf(buf, "f %d//%d %d//%d %d//%d", mesh.pIndices+i, nrm, mesh.pIndices+i+1, nrm+1, mesh.pIndices+i+2, nrm+2) == 6 );
		}
		if ( !s)
		{
			int nrm[3], uv[3];
			s = (sscanf(buf, "f %d/%d/%d %d/%d/%d %d/%d/%d", mesh.pIndices+i, uv, nrm, mesh.pIndices+i+1, uv+1, nrm+1, mesh.pIndices+i+2, uv+2, nrm+2) == 9 );
		}
		if ( !s )
		{
			printf("obj file read error.\n");
			return -2;
		}
		mesh.pIndices[i] += -1;
		mesh.pIndices[i+1] += -1;
		mesh.pIndices[i+2] += -1;
		//i0 = 3 * mesh.pIndices[i];
		//i1 = 3 * mesh.pIndices[i+1];
		//i2 = 3 * mesh.pIndices[i+2];
		//VEC3_SUB(v1, mesh.pVertices+i1, mesh.pVertices+i0);
		//VEC3_SUB(v2, mesh.pVertices+i2, mesh.pVertices+i0);
		//VEC3_CROSS(n, v1, v2);
		//VEC3_ADD(mesh.pNormals+i0, mesh.pNormals+i0, n);
		//VEC3_ADD(mesh.pNormals+i1, mesh.pNormals+i1, n);
		//VEC3_ADD(mesh.pNormals+i2, mesh.pNormals+i2, n);
	}
	fclose(pFile);

	/* normalize vertex normals */
	for(i=0; i<mesh.iNumVertices; i+=3)
		VEC3_NORMALIZE(mesh.pNormals+i);
	mesh.iNumIndices /= 3;
	mesh.iNumVertices /= 3;


	mesh.min[0] = mesh.max[0] = mesh.pVertices[0];
	mesh.min[1] = mesh.max[1] = mesh.pVertices[1];
	mesh.min[2] = mesh.max[2] = mesh.pVertices[2];

	for ( int i = 1; i < mesh.iNumVertices; i++ )
	{
		if ( mesh.pVertices[3*i] < mesh.min[0] ) mesh.min[0] = mesh.pVertices[3*i];
		if ( mesh.pVertices[3*i+1] < mesh.min[1] ) mesh.min[1] = mesh.pVertices[3*i+1];
		if ( mesh.pVertices[3*i+2] < mesh.min[2] ) mesh.min[2] = mesh.pVertices[3*i+2];
		if ( mesh.pVertices[3*i] > mesh.max[0] ) mesh.max[0] = mesh.pVertices[3*i];
		if ( mesh.pVertices[3*i+1] > mesh.max[1] ) mesh.max[1] = mesh.pVertices[3*i+1];
		if ( mesh.pVertices[3*i+2] > mesh.max[2] ) mesh.max[2] = mesh.pVertices[3*i+2];
	}
	return 0;
}

void save_mesh(char* filename, MeshV* mesh, Image* image)
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
		double rgb[3];
		if ( image != NULL )
		{
			double u = mesh->pUV[2*i];
			double v = mesh->pUV[2*i+1];

			int iu = (int)(u*(double)(image->sizeX-1));
			int iv = (int)(v*(double)(image->sizeY-1));
			rgb[0] = (double)image->data[(3*image->sizeX*iv + 3*iu)];
			rgb[1] = (double)image->data[(3*image->sizeX*iv + 3*iu)+1];
			rgb[2] = (double)image->data[(3*image->sizeX*iv + 3*iu)+2];

			if ( rgb[0] < 0 ) rgb[0] = 128.0;
			if ( rgb[1] < 0 ) rgb[1] = 128.0;
			if ( rgb[2] < 0 ) rgb[2] = 128.0;
			fprintf(fp, "v %f %f %f %d %d %d\n", 
				mesh->pVertices[3*i], mesh->pVertices[3*i+1], mesh->pVertices[3*i+2],
				(int)(rgb[0]/128.0)*255, (int)(rgb[1]/128.0)*255, (int)(rgb[2]/128.0)*255);
		}else
		{
#if 10
			rgb[0] = 0.5*mesh->pNormals[3*i]+0.5;
			rgb[1] = 0.5*mesh->pNormals[3*i+1]+0.5;
			rgb[2] = 0.5*mesh->pNormals[3*i+2]+0.5;
			fprintf(fp, "v %f %f %f %d %d %d\n", 
				mesh->pVertices[3*i], mesh->pVertices[3*i+1], mesh->pVertices[3*i+2],
				(int)(rgb[0]*255.0), (int)(rgb[1]*255.0), (int)(rgb[2]*255.0));
#else
			fprintf(fp, "v %f %f %f\n", 
				mesh->pVertices[3*i], mesh->pVertices[3*i+1], mesh->pVertices[3*i+2]);
#endif
		}
	}

	if ( mesh->pUV )
	{
		for ( int i = 0; i < mesh->iNumVertices; i++ )
		{
			fprintf(fp, "vt %f %f\n", mesh->pUV[2*i], mesh->pUV[2*i+1]);
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
int main(int argc, char** argv)
{
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	_splitpath(argv[1], drive_, dir_, fname_, ext_ );
	char* fn = argv[1];
	//fn = "C:\\Documents and Settings\\neut_yamato\\デスクトップ\\yama_tmp\\aaa\\opengi\\projects\\vs2005\\genus2_out.obj";
	//_splitpath(fn, drive_, dir_, fname_, ext_ );

	char outfile[256];

	MeshV mesh;
	int stat = create_mesh_obj(fn, mesh);
	if ( stat != 0 )
	{
		return -1;
	}

	double org[3];
	double u, v;
	double len[3];

	double min[3], max[3];

	len[0] = (mesh.max[0] - mesh.min[0]);
	len[1] = (mesh.max[1] - mesh.min[1]);
	len[2] = (mesh.max[2] - mesh.min[2]);

	min[0] = mesh.min[0]-len[0]*1.2;
	min[1] = mesh.min[1]-len[1]*1.2;
	min[2] = mesh.min[2]-len[2]*1.2;

	max[0] = mesh.max[0]+len[0]*1.2;
	max[1] = mesh.max[1]+len[1]*1.2;
	max[2] = mesh.max[2]+len[2]*1.2;

	len[0] = (max[0] - min[0]);
	len[1] = (max[1] - min[1]);
	len[2] = (max[2] - min[2]);

	org[0] = min[0];
	org[1] = min[1];
	org[2] = (max[2] + min[2])*0.5;

	MeshV mesh2;

	mesh2.iNumIndices = mesh.iNumIndices;
	mesh2.iNumVertices = mesh.iNumVertices;
	mesh2.pIndices = mesh.pIndices;
	mesh2.pColors = mesh.pColors;
	mesh2.pUV = mesh.pUV;
	mesh2.pVertices = new double[mesh.iNumVertices*3];
	mesh2.pNormals = mesh.pNormals;

	char dirname[256];
	sprintf(dirname, "%s%s\\%s", drive_, dir_, fname_); 
	CreateDirectoryA(dirname , NULL);
	printf("CreateDirectory => [%s]\n", dirname);

	int n = 30;
	double t = 1.0/(double)n;
	for ( int j = 0; j <= n; j++ )
	{
		Image* image = NULL;
		if ( j == 0 )
		{
			//image = new Image;
			//ImageLoad("c:\\tmp\\checker_1k.bmp", image);
		}
		for ( int i = 0; i < mesh.iNumVertices; i++ )
		{
			double *uv = mesh.UVCoord(i);
			double* p = mesh.VertexCoord(i);;

			double end[3];
			end[0] = org[0] + uv[0]*len[0];
			end[1] = org[1] + uv[1]*len[1];
			end[2] = org[2];

			double q[3];
			double tt = (double)j*t;
			q[0] = (1.0 - tt)*p[0] +  tt*end[0];
			q[1] = (1.0 - tt)*p[1] +  tt*end[1];
			q[2] = (1.0 - tt)*p[2] +  tt*end[2];

			mesh2.pVertices[3*i] = q[0];
			mesh2.pVertices[3*i+1] = q[1];
			mesh2.pVertices[3*i+2] = q[2];
		}
		sprintf(outfile, "%s\\%s_(%03d).obj", dirname, fname_, j);
		save_mesh(outfile, &mesh2, image);
	}
	mesh2.pIndices = 0;
	mesh2.pColors = 0;
	mesh2.pUV = 0;
	mesh2.pNormals = 0;

	return 0;
}
