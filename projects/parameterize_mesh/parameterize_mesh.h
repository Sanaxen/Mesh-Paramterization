#ifndef _parameterize_mesh_h

#define _parameterize_mesh_h

typedef struct MeshVertexList
{
	unsigned int id;
	double b_min[3];
	double b_max[3];
	unsigned int *pIndices;
	int iNumVertices;
	int iNumIndices;
	double *pVertices;
	double *pNormals;
	double *pUV;
	void* boundaryLine;
} MeshVertexList_t;

typedef struct PolygonSurface_type
{
	MeshVertexList_t* mesh_;
	void* ps;
	void* parameterCurves;
} PolygonSurface_t;

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _DLLAPI
#define _DLLAPI extern "C" __declspec(dllimport) 
#endif

_DLLAPI void* initialize_parameterize_mesh();
_DLLAPI void terminate_parameterize_mesh(void* env);

_DLLAPI void delete_mesh(MeshVertexList_t* mesh);

_DLLAPI int create_mesh_i(char* filename, MeshVertexList_t* mesh, double scale);
_DLLAPI unsigned int create_mesh_io(MeshVertexList_t* mesh);
_DLLAPI unsigned int create_mesh(const char *filename, MeshVertexList_t* mesh, double scale);
_DLLAPI unsigned int create_mesh_obj(const char *filename, MeshVertexList_t* mesh, double scale);
_DLLAPI unsigned int create_mesh_m(const char *filename, MeshVertexList_t* mesh, double scale);
_DLLAPI unsigned int create_mesh_off(const char *filename, MeshVertexList_t* mesh, double scale);

_DLLAPI void save_mesh(char* filename, MeshVertexList_t* mesh);
_DLLAPI void save_mesh_off(char* filename, MeshVertexList_t* mesh);
_DLLAPI void save_mesh2(char* filename, MeshVertexList_t* mesh);
_DLLAPI void save_boundary(char* filename, MeshVertexList_t* mesh, double radius, double scale);
_DLLAPI int parameterize_mesh_file(char* infile, char* outfile, int Iter, double stretch_coef, int type);
_DLLAPI int parameterize_mesh_io(MeshVertexList_t* mesh, int Iter, double stretch_coef, int type, double scale );

_DLLAPI PolygonSurface_t* CreateSurfaceMesh(MeshVertexList_t* mesh);
_DLLAPI void DeleteSurfaceMesh(PolygonSurface_t* ps);
_DLLAPI void SurfaceMesh(PolygonSurface_t* ps, double* uv, double point[3]);
_DLLAPI void CreateParameterCurvs(PolygonSurface_t* ps, int n,  int polydivn);
_DLLAPI void save_ParameterCurvs(char* filename, PolygonSurface_t* ps, double radius);
_DLLAPI void ParameterMesh(PolygonSurface_t* ps, int m,  int n, MeshVertexList_t* mesh);


_DLLAPI void* new_mesh_p(void);
_DLLAPI void delete_mesh_p(void* mesh);
_DLLAPI int create_mesh_i_p(char* filename, void* mesh);
_DLLAPI void save_mesh_p(char* filename, void* mesh);
_DLLAPI void save_boundary_p(char* filename, void* mesh, double radius);
_DLLAPI int parameterize_mesh_io_p(void* mesh, int Iter, double stretch_coef );

_DLLAPI void* CreateSurfaceMesh_p(void* mesh);
_DLLAPI void DeleteSurfaceMesh_p(void* p);
_DLLAPI void SurfaceMesh_p(void* p, double* uv, double point[3]);
_DLLAPI void CreateParameterCurvs_p(void* p, int n,  int polydivn);
_DLLAPI void save_ParameterCurvs_p(char* filename, void* p, double radius);
_DLLAPI void ParameterMesh_p(void* ps, int m,  int n, void* mesh);

#ifdef __cplusplus
}
#endif


#endif