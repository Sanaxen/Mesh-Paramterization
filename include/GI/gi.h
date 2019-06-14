/*
 *  OpenGI: Library for Parameterization and Geometry Image creation
 *  Copyright (C) 2008-2011  Christian Rau
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Christian Rau
 *
 *     rauy@users.sourceforge.net
 */

/** \file
 *  The main header file of the OpenGI library.
 */

#ifndef __GI_H__
#define __GI_H__

//#ifdef __cplusplus
//extern "C" {
//#endif

#define GI_VERSION_2_0
#define GI_VERSION_2_1

#ifdef _WIN32
	#if defined(__GNU_C__) || defined(OPENGI_STATIC)
		#define GIAPI extern
	#else
		#if defined(OPENGI_EXPORTS)
			#define GIAPI __declspec(dllexport)
		#else
			#define GIAPI __declspec(dllimport)
		#endif
	#endif
	#define GIAPIENTRY
	#define GICALLBACK
#else
	#define GIAPI extern
	#define GIAPIENTRY
	#define GICALLBACK
#endif

#include "../src/my_types.h"
/*************************************************************************/
/* Typedefs */

/* Standard data types */
typedef unsigned int   GIenum;						/**< Enumaration type. */
typedef unsigned char  GIboolean;					/**< Boolean type. */
typedef unsigned int   GIbitfield;					/**< Bitfield type. */
typedef char           GIchar;						/**< Character type. */
typedef signed char    GIbyte;						/**< 1 byte signed integer type. */
typedef short          GIshort;						/**< 2 byte signed integer type. */
typedef int            GIint;						/**< 4 byte signed integer type. */
typedef int            GIsizei;						/**< Size type. */
typedef unsigned char  GIubyte;						/**< 1 byte unsigned integer type. */
typedef unsigned short GIushort;					/**< 2 byte unsigned integer type. */
typedef unsigned int   GIuint;						/**< 4 byte unsigned integer type. */
typedef GLfloat         GIfloat;						/**< Single precission floating point type. */
typedef GLdouble         GIdouble;					/**< Double precission floating point type. */
typedef unsigned short GIhalf;						/**< Half precission floating point type. */
typedef void           GIvoid;						/**< Void type. */
typedef void*          GIcontext;					/**< Context type. */

/** Parameterization callback. */
typedef GIboolean (GICALLBACK *GIparamcb)(GIvoid*);

/** Sampler callback. */
typedef void (GICALLBACK *GIsamplercb)(GIdouble*, const GIdouble*, GIvoid*);

/** Error callback. */
typedef void (GICALLBACK *GIerrorcb)(GIenum, GIvoid*);


/*************************************************************************/
/* Constants */

/** \name Standard constants
 *  \{
 */
#define GI_TRUE                          1			/**< Boolean value true. */
#define GI_FALSE                         0			/**< Boolean value false. */
#define GI_NULL                          0			/**< Zero value. */
#define GI_NONE                          0			/**< Zero identifier. */
#define GI_ALL_PATCHES                   -1			/**< Whole mesh. */
/** \} */

/** \name Data types (OpenGL compatible)
 *  \{
 */
#define GI_BYTE                          0x1400		/**< Signed byte. */
#define GI_UNSIGNED_BYTE                 0x1401		/**< Unsigned byte. */
#define GI_SHORT                         0x1402		/**< Signed short integer. */
#define GI_UNSIGNED_SHORT                0x1403		/**< Unsigned short integer. */
#define GI_INT                           0x1404		/**< Signed integer. */
#define GI_UNSIGNED_INT                  0x1405		/**< Unsigned integer. */
#define GI_FLOAT                         0x1406		/**< Single precision float. */
#define GI_DOUBLE                        0x140A		/**< Double precision float. */
#define GI_HALF_FLOAT                    0x140B		/**< Half precision float. */
/** \} */

/** \name Enable flags
 *  \{
 */
#define GI_MULTITHREADING                0x0100		/**< Use multithreading. */
#define GI_EXACT_MAPPING_SUBSET          0x0110		/**< Vertices mapped to exact texels. */
#define GI_PARAM_CORNER_SUBSET           0x0111		/**< Vertices mapped to parameter domain corners. */
/** \} */

/** \name State parameters
 *  \{
 */
#define GI_VERSION                       0x0200		/**< OpenGI version. */
#define GI_MAX_ATTRIBS                   0x0201		/**< Number of supported custom attributes. */
#define GI_MESH_BINDING                  0x0202		/**< Current bound mesh. */
#define GI_IMAGE_BINDING                 0x0203		/**< Current bound image for active texture unit. */
#define GI_SAMPLED_ATTRIB_COUNT          0x0204		/**< Number of sampled attributes. */
#define GI_SAMPLED_ATTRIBS               0x0205		/**< Flags indicating sampled attributes. */
#define GI_EXACT_MAPPING_SUBSET_COUNT    0x0210		/**< Number of elements in exact mapping subset. */
#define GI_PARAM_CORNER_SUBSET_COUNT     0x0211		/**< Number of elements in corner subset. */
#define GI_EXACT_MAPPING_SUBSET_SORTED   0x0220		/**< Sorted flag for exact mapping subset. */
#define GI_PARAM_CORNER_SUBSET_SORTED    0x0221		/**< Sorted flag for corner subset. */
/** \} */

/** \name Special attribute semantics
 *  \{
 */
#define GI_POSITION_ATTRIB               0x0300		/**< Vertex position. */
#define GI_PARAM_ATTRIB                  0x0301		/**< Parameter coordinates. */
#define GI_PARAM_STRETCH_ATTRIB          0x0302		/**< Per-param stretch values. */
/** \} */

/** \name Attribute state
 *  \{
 */
#define GI_ATTRIB_ARRAY_ENABLED          0x0400		/**< Enabled flag. */
#define GI_ATTRIB_ARRAY_SIZE             0x0401		/**< Number of components. */
#define GI_ATTRIB_ARRAY_STRIDE           0x0402		/**< Stride value. */
#define GI_ATTRIB_ARRAY_NORMALIZED       0x0403		/**< Normalization flag. */
#define GI_ATTRIB_ARRAY_SEMANTIC         0x0404		/**< Semantic if bound to any. */
#define GI_ATTRIB_IMAGE                  0x0405		/**< Currently bound image. */
/** \} */

/** \name Mesh and patch properties
 *  \{
 */
#define GI_HAS_CUT                       0x0500		/**< Mesh has valid cut. */
#define GI_HAS_PARAMS                    0x0501		/**< Mesh has UVs. */
#define GI_PATCH_COUNT                   0x0502		/**< Number of patches. */
#define GI_FACE_COUNT                    0x0503		/**< Number of faces. */
#define GI_EDGE_COUNT                    0x0504		/**< Number of edges. */
#define GI_VERTEX_COUNT                  0x0505		/**< Number of vertices. */
#define GI_PARAM_COUNT			         0x0506		/**< Number of UVs. */
#define GI_AABB_MIN                      0x0507		/**< Minimal vertex of AABB. */
#define GI_AABB_MAX                      0x0508		/**< Maximal vertex of AABB. */
#define GI_RADIUS                        0x0509		/**< Radius of bounding sphere around origin. */
#define GI_ACTIVE_PATCH                  0x050A		/**< Active patch. */
#define GI_MIN_PARAM_STRETCH             0x050B		/**< Minimum param stretch value. */
#define GI_MAX_PARAM_STRETCH             0x050C		/**< Maximum param stretch value. */
#define GI_PARAM_STRETCH_METRIC          0x050D		/**< Metric of param stretch values. */
#define GI_TOPOLOGICAL_SIDEBAND_LENGTH   0x050E		/**< Length of topological sideband. */
#define GI_TOPOLOGICAL_SIDEBAND          0x050F		/**< Topological sideband (cut path information). */
/** \} */

/** \name Per-mesh attribute state
 *  \{
 */
#define GI_HAS_ATTRIB                    0x0600		/**< Mesh has attrib. */
#define GI_ATTRIB_SIZE                   0x0601		/**< Number of components. */
#define GI_ATTRIB_NORMALIZED             0x0602		/**< Normalization flag. */
#define GI_ATTRIB_SEMANTIC               0x0603		/**< Attribute semantic. */
/** \} */

/** \name Cutting parameters
 *  \{
 */
#define GI_CUTTER                        0x0700		/**< Cutting algorithm. */
#define GI_SUBDIVISION_ITERATIONS        0x0701		/**< Iterations for subdivision algorithms. */
//#define GI_STRAIGHTEN_CUT                0x0702		/**< Straighten cut paths. */
//#define GI_ORIENTATION_WEIGHT            0x0703		/**< Orientation bias weight for face clustering. */
//#define GI_SHAPE_WEIGHT                  0x0704		/**< Shape bias weight for face clustering. */
#define GI_INITIAL_GIM                   0x0710		/**< Initial cut of original GIM algorithm. */
#define GI_CATMULL_CLARK_SUBDIVISION     0x0711		/**< Catmull-Clark subdivision patchification. */
//#define GI_FACE_CLUSTERING               0x0712		/**< Quadric-based face clustering. */
/** \} */

/** \name Parameterization parameters
 *  \{
 */
#define GI_PARAMETERIZER                 0x0800		/**< Parameterization algorithm. */
#define GI_INITIAL_PARAMETERIZATION      0x0801		/**< Initial parameterization for stretch minimizer. */
#define GI_STRETCH_METRIC                0x0802		/**< Stretch metric to minimize. */
#define GI_CONFORMAL_WEIGHT              0x0803		/**< Conformal weight for intrinsic parameterization. */
#define GI_AUTHALIC_WEIGHT               0x0804		/**< Authalic weight for intrinsic parameterization. */
#define GI_STRETCH_WEIGHT                0x0805		/**< Eta parameter for stretch minimizer. */
#define GI_PARAM_RESOLUTION              0x0806		/**< Desired minimal resolution for sampling. */
#define GI_UNSYMMETRIC_SOLVER            0x0807		/**< Solver for unsymmetric systems. */
#define GI_AREA_WEIGHT                   0x0808		/**< Theta parameter for combined energy. */
#define GI_PARAM_SOURCE_ATTRIB           0x0809		/**< Attribute to use as parameter coords. */
#define GI_FROM_ATTRIB                   0x0810		/**< Set attrib as parameter coordinates. */
#define GI_TUTTE_BARYCENTRIC             0x0811		/**< Tutte's Barycentric parameterization. */
#define GI_SHAPE_PRESERVING              0x0812		/**< Floater's Shape Preserving parameterization. */
#define GI_DISCRETE_HARMONIC             0x0813		/**< Eck's Discrete Harmonic parameterization. */
#define GI_MEAN_VALUE                    0x0814		/**< Floater's Mean Value parameterization. */
#define GI_DISCRETE_AUTHALIC             0x0815		/**< Desbrun's Discrete Authalic parameterization. */
#define GI_INTRINSIC                     0x0816		/**< Desbrun's Intrinsic parameterization. */
#define GI_STRETCH_MINIMIZING            0x0817		/**< Yoshizawa's stretch minimizing parameterization. */
#define GI_GIM                           0x0818		/**< Gu's original Geometry Image parameterization. */
#define GI_SOLVER_BICGSTAB               0x0820		/**< BiCGStab solver. */
#define GI_SOLVER_GMRES                  0x0821		/**< GMRES(m) solver. */
#define GI_PARAM_STARTED                 0x0830		/**< Callback for parameterization start. */
#define GI_PARAM_CHANGED                 0x0831		/**< Callback for parameterization change. */
#define GI_PARAM_FINISHED                0x0832		/**< Callback for parameterization end. */
/** \} */

/** \name Stretch metrics
 *  \{
 */
#define GI_MAX_GEOMETRIC_STRETCH         0x0900		/**< Linf stretch as defined by Sander et al. */
#define GI_RMS_GEOMETRIC_STRETCH         0x0901		/**< L2 stretch as defined by Sander et al. */
#define GI_COMBINED_STRETCH              0x0902		/**< Combined energy as defined by Degener et al. */
/** \} */

/** \name Image properties
 *  \{
 */
#define GI_IMAGE_WIDTH                   0x0A00		/**< Width of image. */
#define GI_IMAGE_HEIGHT                  0x0A01		/**< Height of image. */
#define GI_IMAGE_COMPONENTS              0x0A02		/**< Number of components. */
#define GI_IMAGE_TYPE                    0x0A03		/**< Data type of image elements. */
#define GI_GL_IMAGE_TEXTURE              0x0A04		/**< OpenGL texture object of image (if texture storage). */
#define GI_GL_IMAGE_BUFFER               0x0A05		/**< OpenGL buffer object of image (if buffer storage). */
#define GI_IMAGE_DATA                    0x0A06		/**< RAM address of image data (if RAM storage). */
#define GI_IMAGE_STORAGE                 0x0A07		/**< Storage mode of image. */
#define GI_SUBIMAGE_X                    0x0A08		/**< X-offset of sub-image. */
#define GI_SUBIMAGE_Y                    0x0A09		/**< X-offset of sub-image. */
#define GI_SUBIMAGE_WIDTH                0x0A0A		/**< Width of sub-image. */
#define GI_SUBIMAGE_HEIGHT               0x0A0B		/**< Height of sub-image. */
#define GI_SUBIMAGE                      0x0A0C		/**< Sub-image rectangle as 4 ints. */
#define GI_EXTERNAL_DATA                 0x0A10		/**< Image stored in user-owned RAM. */
#define GI_GL_TEXTURE_DATA               0x0A11		/**< Image stored in user-owned GL texture. */
#define GI_GL_BUFFER_DATA                0x0A12		/**< Image stored in user-owned GL buffer. */
#define GI_NO_IMAGE_DATA                 0x0A13		/**< Image data not yet specified. */
/** \} */

/** \name Sampling parameters
 *  \{
 */
#define GI_SAMPLER                       0x0B00		/**< Sampler to use. */
#define GI_SAMPLER_USE_SHADER            0x0B01		/**< Use GLSL if supported. */
#define GI_SAMPLER_USE_RENDER_TO_TEXTURE 0x0B02		/**< Use framebuffer objects if supported. */
#define GI_SAMPLER_SOFTWARE              0x0B10		/**< Software sampling. */
#define GI_SAMPLER_OPENGL                0x0B11		/**< Hardware accelerated sampling using OpenGL. */
/** \} */

/** \name Per-attribute sampling parameters
 *  \{
 */
#define GI_SAMPLING_MODE                 0x0C00		/**< Sampling style. */
#define GI_SAMPLING_TRANSFORM            0x0C01		/**< Transformation matrix. */
#define GI_TEXTURE_DIMENSION             0x0C02		/**< Dimension of OpenGL texture. */
#define GI_GL_SAMPLE_TEXTURE             0x0C03		/**< OpenGL texture to resample. */
#define GI_SAMPLE_DEFAULT                0x0000		/**< Standard sampling. */
#define GI_SAMPLE_NORMALIZED             0x0C10		/**< Normalize after interpolation. */
#define GI_SAMPLE_TEXTURED               0x0C11		/**< Resample texture. */
/** \} */

/** \name Render parameters
 *  \{
 */
#define GI_GL_USE_VERTEX_TEXTURE         0x0001		/**< Use vertex texturing for GIM rendering. */
#define GI_GL_USE_GEOMETRY_SHADER        0x0002		/**< Use geometry shader for GIM rendering. */
#define GI_RENDER_RESOLUTION_U           0x0D03		/**< U-resolution of base mesh for GIM rendering. */
#define GI_RENDER_RESOLUTION_V           0x0D04		/**< V-resolution of base mesh for GIM rendering. */
#define GI_RENDER_CACHE_SIZE             0x0D05		/**< Size of cache for geometry image rendering. */
/** \} */

/** \name Per-attribute render parameters
 *  \{
 */
#define GI_GL_RENDER_SEMANTIC            0x0E00		/**< Render semantic. */
#define GI_GL_RENDER_CHANNEL             0x0E01		/**< Render channel (for attribs and texCoords). */
#define GI_TEXTURE_COORD_DOMAIN          0x0E02		/**< TexCoord domain for images. */
#define GI_GL_VERTEX                     0x0E10		/**< Use attribute as glVertex. */
#define GI_GL_NORMAL                     0x0E11		/**< Use attribute as glNormal. */
#define GI_GL_COLOR                      0x0E12		/**< Use attribute as glColor. */
#define GI_GL_SECONDARY_COLOR            0x0E13		/**< Use attribute as glSecondaryColor. */
#define GI_GL_FOG_COORD                  0x0E14		/**< Use attribute as glFogCoord. */
#define GI_GL_EVAL_COORD                 0x0E15		/**< Use attribute as glEvalCoord. */
#define GI_GL_TEXTURE_COORD              0x0E16		/**< Use attribute as gl[Multi]TexCoord. */
#define GI_GL_VERTEX_ATTRIB              0x0E17		/**< Use attribute as glVertexAttrib. */
#define GI_UNIT_SQUARE                   0x0E20		/**< 0-1-square. */
#define GI_HALF_TEXEL_INDENT             0x0E21		/**< 0-1-square indented by half a texel. */
/** \} */

/** \name Error codes
 *  \{
 */
#define GI_NO_ERROR                      0x0000		/**< No error encountered. */
#define GI_INVALID_ENUM                  0x1000		/**< Wrong enum value. */
#define GI_INVALID_OPERATION             0x1001		/**< Wrong operation in current state. */
#define GI_INVALID_VALUE                 0x1002		/**< Numeric value out of range. */
#define GI_INVALID_ID                    0x1003		/**< Object ID not existent. */
#define GI_INVALID_MESH                  0x1004		/**< Non-manifold mesh. */
#define GI_INVALID_CUT                   0x1005		/**< Cut incomplete. */
#define GI_NUMERICAL_ERROR               0x1006		/**< Error in numerical routine. */
#define GI_UNSUPPORTED_OPERATION         0x1007		/**< Feature not supported. */
#define GI_INVALID_PARAMETERIZATION      0x1008		/**< No unit square. */
/** \} */


/*************************************************************************/
/* Functions */

/** \name Context management
 *  \{
 */
GIAPI GIcontext     GIAPIENTRY giCreateContext();
GIAPI void          GIAPIENTRY giMakeCurrent(GIcontext context);
GIAPI GIcontext     GIAPIENTRY giGetCurrent();
GIAPI void          GIAPIENTRY giDestroyContext(GIcontext context);
/** \} */

/** \name State management
 *  \{
 */
GIAPI GIenum        GIAPIENTRY giGetEnumValue(const GIchar *name);
GIAPI void          GIAPIENTRY giGetBooleanv(GIenum pname, GIboolean *params);
GIAPI void          GIAPIENTRY giGetIntegerv(GIenum pname, GIint *params);
GIAPI void          GIAPIENTRY giGetFloatv(GIenum pname, GIfloat *params);
GIAPI void          GIAPIENTRY giEnable(GIenum pname);
GIAPI void          GIAPIENTRY giDisable(GIenum pname);
GIAPI GIboolean     GIAPIENTRY giIsEnabled(GIenum pname);
GIAPI void          GIAPIENTRY giVertexSubset(GIenum subset, GIsizei count, GIboolean sorted, const GIuint *indices);
GIAPI void          GIAPIENTRY giGetPointerv(GIenum pname, GIvoid **params);
GIAPI void          GIAPIENTRY giAttribImage(GIuint attrib, GIuint image);
GIAPI void          GIAPIENTRY giBindAttrib(GIenum semantic, GIuint attrib);
GIAPI void          GIAPIENTRY giEnableAttribArray(GIuint attrib);
GIAPI void          GIAPIENTRY giDisableAttribArray(GIuint attrib);
GIAPI void          GIAPIENTRY giAttribPointer(GIuint attrib, GIint size, GIboolean normalized, GIsizei stride, GIfloat *pointer);
GIAPI void          GIAPIENTRY giGetAttribbv(GIuint attrib, GIenum pname, GIboolean *params);
GIAPI void          GIAPIENTRY giGetAttribiv(GIuint attrib, GIenum pname, GIint *params);
GIAPI void          GIAPIENTRY giGetAttribfv(GIuint attrib, GIenum pname, GIfloat *params);
GIAPI void          GIAPIENTRY giGetAttribPointerv(GIuint attrib, GIvoid **params);
/** \} */

/** \name Mesh handling
 *  \{
 */
GIAPI GIuint        GIAPIENTRY giGenMesh();
GIAPI void          GIAPIENTRY giGenMeshes(GIsizei n, GIuint *meshes);
GIAPI void          GIAPIENTRY giBindMesh(GIuint mesh);
GIAPI void          GIAPIENTRY giDeleteMesh(GIuint mesh);
GIAPI void          GIAPIENTRY giDeleteMeshes(GIsizei n, const GIuint *meshes);
GIAPI void          GIAPIENTRY giIndexedMesh(GIuint start, GIuint end, GIsizei count, const GIuint *indices);
GIAPI void          GIAPIENTRY giNonIndexedMesh(GIint first, GIsizei count);
GIAPI void          GIAPIENTRY giCopyMesh(GIuint mesh);
GIAPI void          GIAPIENTRY giComputeParamStretch(GIenum metric);
GIAPI void          GIAPIENTRY giMeshActivePatch(GIint patch);
GIAPI void          GIAPIENTRY giGetMeshbv(GIenum pname, GIboolean *params);
GIAPI void          GIAPIENTRY giGetMeshiv(GIenum pname, GIint *params);
GIAPI void          GIAPIENTRY giGetMeshfv(GIenum pname, GIfloat *params);
GIAPI void          GIAPIENTRY giGetMeshAttribbv(GIuint attrib, GIenum pname, GIboolean *params);
GIAPI void          GIAPIENTRY giGetMeshAttribiv(GIuint attrib, GIenum pname, GIint *params);
GIAPI void          GIAPIENTRY giGetIndexedMesh(GIuint *vcount, GIuint *icount, GIuint *indices);
GIAPI void          GIAPIENTRY giGetNonIndexedMesh(GIuint *vcount);
/** \} */

/** \name Cutting
 *  \{
 */
//GIAPI void          GIAPIENTRY giCutterParameterb(GIenum pname, GIboolean param);
GIAPI void          GIAPIENTRY giCutterParameteri(GIenum pname, GIint param);
//GIAPI void          GIAPIENTRY giCutterParameterf(GIenum pname, GIfloat param);
GIAPI void          GIAPIENTRY giCut();
/** \} */

/** \name Parameterization
 *  \{
 */
GIAPI void          GIAPIENTRY giParameterizerParameterb(GIenum pname, GIboolean param);
GIAPI void          GIAPIENTRY giParameterizerParameteri(GIenum pname, GIint param);
GIAPI void          GIAPIENTRY giParameterizerParameterf(GIenum pname, GIfloat param);
GIAPI void          GIAPIENTRY giParameterizerCallback(GIenum which, GIparamcb fn, GIvoid *data);
GIAPI void          GIAPIENTRY giParameterize();
/** \} */

/** \name Image handling
 *  \{
 */
GIAPI GIuint        GIAPIENTRY giGenImage();
GIAPI void          GIAPIENTRY giGenImages(GIsizei n, GIuint *images);
GIAPI void          GIAPIENTRY giBindImage(GIuint image);
GIAPI void          GIAPIENTRY giDeleteImage(GIuint image);
GIAPI void          GIAPIENTRY giDeleteImages(GIsizei n, const GIuint *images);
GIAPI void          GIAPIENTRY giImageExternalData(GIsizei width, GIsizei height, GIsizei components, GIenum type, GIvoid *data);
GIAPI void          GIAPIENTRY giImageGLTextureData(GIsizei width, GIsizei height, GIsizei components, GIenum type, GIuint texture);
GIAPI void          GIAPIENTRY giImageGLBufferData(GIsizei width, GIsizei height, GIsizei components, GIenum type, GIuint buffer);
GIAPI void          GIAPIENTRY giSubImage(GIuint x, GIuint y, GIsizei width, GIsizei height);
GIAPI void          GIAPIENTRY giGetImageiv(GIenum pname, GIint *param);
/** \} */

/** \name Sampling
 *  \{
 */
GIAPI void          GIAPIENTRY giSamplerParameterb(GIenum pname, GIboolean param);
GIAPI void          GIAPIENTRY giSamplerParameteri(GIenum pname, GIint param);
GIAPI void          GIAPIENTRY giAttribSamplerParameteri(GIuint attrib, GIenum pname, GIint param);
GIAPI void          GIAPIENTRY giAttribSamplerParameterfv(GIuint attrib, GIenum pname, const GIfloat *params);
GIAPI void          GIAPIENTRY giSample();
/** \} */

/** \name Error handling
 *  \{
 */
GIAPI GIenum        GIAPIENTRY giGetError();
GIAPI const GIchar* GIAPIENTRY giErrorString(GIenum errorCode);
GIAPI void          GIAPIENTRY giErrorCallback(GIerrorcb fn, GIvoid *data);
/** \} */

/** \name OpenGL communication
 *  \{
 */
GIAPI void          GIAPIENTRY giGLRenderParameterb(GIenum pname, GIboolean param);
GIAPI void          GIAPIENTRY giGLRenderParameteri(GIenum pname, GIint param);
GIAPI void          GIAPIENTRY giGLAttribRenderParameteri(GIuint attrib, GIenum pname, GIint param);
GIAPI void          GIAPIENTRY giGLDrawMesh();
GIAPI void          GIAPIENTRY giGLDrawCut();
GIAPI void          GIAPIENTRY giGLDrawGIM();
GIAPI void          GIAPIENTRY giGLCleanUp();

GIAPI void			GIAPIENTRY putMesh(char* filename, int* iNumIndices, int* iNumVertices, unsigned int **pIndices, GIfloat **pVertices, GIfloat **pUV, GIfloat **pNormals, double scale);

GIAPI void			GIAPIENTRY calcBoundaryBox(double min[3], double max[3]);
GIAPI void			GIAPIENTRY calcBoundaryBox2(GIdouble min[3], GIdouble max[3]);

#include "../src/cylinder.h"
GIAPI void			GIAPIENTRY putCut(double r, N_Cylinder& boundaryLine, double scale);
/** \} */


//#ifdef __cplusplus
//}
//#endif

#endif
