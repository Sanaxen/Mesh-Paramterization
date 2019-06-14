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

/** \internal
 *  \file
 *  \brief Implementation of structures and functions for communication with OpenGL.
 */

#include "gi_gl.h"
#include "gi_context.h"
#include "gi_mesh.h"
#include "gi_cutter.h"
#include "gi_container.h"
#include "gi_memory.h"

#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <stdarg.h>

#include "my_gl.h"
#include <time.h>

#define GI_QUADS				0x80000000

//#define GI_SINGLE_TRISTRIP		1

#if 0
/** \internal
 *  \brief Vertex shader for sampling.
 *  \ingroup opengl
 */
static const GIchar *g_sample_vs = 
	"void main() { \n"
	"gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0; \n"
	"gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex; }";

/** \internal
 *  \brief Fragment shader for sampling.
 *  \ingroup opengl
 */
static const GIchar *g_sample_fs = 
	"#if defined(TEXTURE1D) \n"
	"	uniform sampler1D tex; \n"
	"#elif defined(TEXTURE2D) \n"
	"	uniform sampler2D tex; \n"
	"#elif defined(TEXTURE3D) \n"
	"	uniform sampler3D tex; \n"
	"#elif defined(TEXTURECUBE) \n"
	"	uniform samplerCube tex; \n"
	"#endif \n"
	"void main() { \n"
	"#if defined(TEXTURE1D) \n"
	"	gl_FragColor = texture1D(tex, gl_TexCoord[0].s); \n"
	"#elif defined(TEXTURE2D) \n"
	"	gl_FragColor = texture2D(tex, gl_TexCoord[0].st); \n"
	"#elif defined(TEXTURE3D) \n"
	"	gl_FragColor = texture3D(tex, gl_TexCoord[0].stp); \n"
	"#elif defined(TEXTURECUBE) \n"
	"	gl_FragColor = textureCube(tex, gl_TexCoord[0].stp); \n"
	"#elif defined(NORMALIZED) \n"
	"	gl_FragColor = vec4(normalize(gl_TexCoord[0].xyz), 1.0); \n"
	"#elif defined(PACKED) \n"
	"	gl_FragColor = vec4(0.5*normalize(gl_TexCoord[0].xyz)+vec3(0.5), 1.0); \n"
	"#else \n"
	"	gl_FragColor = gl_TexCoord[0]; \n"
	"#endif \n }";

/** \internal
 *  \brief Vertex shader for geometry image rendering.
 *  \ingroup opengl
 */
static const GIchar *g_render_vs = 
	"#ifdef TEXTURE \n"
	"	#define MultiTexCoord0 gl_Vertex \n"
	"	#ifndef GSHADER \n"
	"		uniform sampler2D geometryImage; \n"
	"	#endif \n"
	"#else \n"
	"	#define MultiTexCoord0 gl_MultiTexCoord0 \n"
	"	#define Vertex gl_Vertex \n"
	"#endif \n"
	"uniform vec4 texScale[4]; \n"
	"uniform vec4 texBias[4]; \n"
	"#ifndef GSHADER \n"
	"	varying vec3 ecPosition; \n"
	"#endif \n"
	"void main() { \n"
	"gl_TexCoord[0] = MultiTexCoord0.stst*texScale[0] + texBias[0]; \n"
	"gl_TexCoord[1] = MultiTexCoord0.stst*texScale[1] + texBias[1]; \n"
	"gl_TexCoord[2] = MultiTexCoord0.stst*texScale[2] + texBias[2]; \n"
	"#ifdef GSHADER \n"
	"	#ifdef TEXTURE \n"
	"		gl_Position = gl_Vertex*texScale[3] + texBias[3]; \n"
	"	#else \n"
	"		gl_Position = gl_Vertex; \n"
	"	#endif \n"
	"#else \n"
	"	#ifdef TEXTURE \n"
	"		vec4 Vertex = texture2D(geometryImage, gl_Vertex.st*texScale[3].st+texBias[3].st); \n"
	"	#endif \n"
	"	gl_Position = gl_ModelViewProjectionMatrix * Vertex; \n"
	"	vec4 eyePos = gl_ModelViewMatrix * Vertex; \n"
	"	ecPosition = eyePos.xyz / eyePos.w; \n"
	"#endif \n"
	"}";

/** \internal
 *  \brief Geometry shader for geometry image rendering.
 *  \ingroup opengl
 */
static const GIchar *g_render_gs = 
	"#extension GL_EXT_geometry_shader4 : require \n"
	"#ifdef TEXTURE \n"
	"	#define PositionIn vertex \n"
	"	uniform sampler2D geometryImage; \n"
	"	uniform vec4 offset[4]; \n"
	"	uniform vec2 sign[4]; \n"
	"#else \n"
	"	#define PositionIn gl_PositionIn \n"
	"#endif \n"
	"varying out vec3 ecPosition; \n"
	"void main() { \n"
	"int i, idx; \n"
	"vec4 eyePos; \n"
	"#ifdef TEXTURE \n"
	"	vec4 vertex[4]; \n"
	"	vertex[0] = texture2D(geometryImage, gl_PositionIn[0].st+sign[0]*offset[3].st); \n"
	"	vertex[2] = texture2D(geometryImage, gl_PositionIn[0].st+sign[2]*offset[3].st); \n"
	"	vertex[1] = texture2D(geometryImage, gl_PositionIn[0].st+sign[1]*offset[3].st); \n"
	"	vertex[3] = texture2D(geometryImage, gl_PositionIn[0].st+sign[3]*offset[3].st); \n"
	"#endif \n"
	"vec3 e = PositionIn[2].xyz - PositionIn[0].xyz; \n"
	"vec3 f = PositionIn[3].xyz - PositionIn[1].xyz; \n"
	"ivec4 strip = (dot(e, e)>dot(f, f)) ? ivec4(0, 1, 3, 2) : ivec4(3, 0, 2, 1); \n"
	"for(i=0; i<4; ++i) { \n"
	"	idx = strip[i]; \n"
	"#ifdef TEXTURE \n"
	"	gl_TexCoord[0] = gl_TexCoordIn[0][0] + sign[idx].xyxy*offset[0]; \n"
	"	gl_TexCoord[1] = gl_TexCoordIn[0][1] + sign[idx].xyxy*offset[1]; \n"
	"	gl_TexCoord[2] = gl_TexCoordIn[0][2] + sign[idx].xyxy*offset[2]; \n"
	"#else \n"
	"	gl_TexCoord[0] = gl_TexCoordIn[idx][0]; \n"
	"	gl_TexCoord[1] = gl_TexCoordIn[idx][1]; \n"
	"	gl_TexCoord[2] = gl_TexCoordIn[idx][2]; \n"
	"#endif \n"
	"	gl_Position = gl_ModelViewProjectionMatrix * PositionIn[idx]; \n"
	"	eyePos = gl_ModelViewMatrix * PositionIn[idx]; \n"
	"	ecPosition = eyePos.xyz / eyePos.w; \n"
	"	EmitVertex(); } }";

/** \internal
 *  \brief Fragment shader for geometry image rendering.
 *  \ingroup opengl
 */
static const GIchar *g_render_fs = 
	"uniform sampler2D normalImage; \n"
	"uniform sampler2D colorImage; \n"
	"uniform sampler2D textureImage[4]; \n"
	"uniform ivec4 enabledAttribsSBTS; \n"
	"uniform ivec4 enabledTextures; \n"
	"uniform vec4 baseColor; \n"
	"varying vec3 ecPosition; \n"
	"void main() { \n"
	"vec4 color = baseColor; \n"
	"if(enabledAttribsSBTS.y != 0) \n"
	"	color = texture2D(colorImage, gl_TexCoord[0].pq); \n"
	"if(enabledAttribsSBTS.x != 0) { \n"
	"	vec3 N = texture2D(normalImage, gl_TexCoord[0].st).rgb; \n"
	"	vec3 L = gl_LightSource[0].position.xyz; \n"
	"	if(gl_LightSource[0].position.w != 0.0) \n"
	"		L -= ecPosition; \n"
	"	L = normalize(L); \n"
	"	if(enabledAttribsSBTS.z != 0) \n"
	"		N = 2.0*N - vec3(1.0); \n"
	"	N = normalize(gl_NormalMatrix*N); \n"
	"	vec4 diffuse, specular; \n"
	"	float shininess; \n"
	"	if(gl_FrontFacing) { \n"
	"		color = gl_FrontLightModelProduct.sceneColor + gl_FrontLightProduct[0].ambient + gl_FrontMaterial.emission; \n"
	"		diffuse = gl_FrontLightProduct[0].diffuse; \n"
	"		specular = gl_FrontLightProduct[0].specular; \n"
	"		shininess = gl_FrontMaterial.shininess; \n"
	"	} else { \n"
	"		color = gl_BackLightModelProduct.sceneColor + gl_BackLightProduct[0].ambient + gl_BackMaterial.emission; \n"
	"		diffuse = gl_BackLightProduct[0].diffuse; \n"
	"		specular = gl_BackLightProduct[0].specular; \n"
	"		shininess = gl_BackMaterial.shininess;  \n"
	"		if(enabledAttribsSBTS.w != 0) \n"
	"			N = -N; } \n"
	"	float NdotL = max(dot(N, L), 0.0); \n"
	"	color += NdotL * diffuse; \n"
	"	if(NdotL > 0.0) { \n"
	"		float NdotH = max(dot(N, normalize(L-normalize(ecPosition))), 0.0); \n"
	"		color += pow(NdotH, shininess) * specular; } } \n"
	"if(enabledTextures.x != 0) \n"
	"	color *= texture2D(textureImage[0], gl_TexCoord[1].st); \n"
	"if(enabledTextures.y != 0) \n"
	"	color *= texture2D(textureImage[1], gl_TexCoord[1].pq); \n"
	"if(enabledTextures.z != 0) \n"
	"	color *= texture2D(textureImage[2], gl_TexCoord[2].st); \n"
	"if(enabledTextures.w != 0) \n"
	"	color *= texture2D(textureImage[3], gl_TexCoord[2].pq); \n"
	"gl_FragColor = color; }";

#endif

/** Set boolean render parameter.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup opengl
 */
void GIAPIENTRY giGLRenderParameterb(GIenum pname, GIboolean param)
{
	GIRenderer *pRenderer = &(GIContext_current()->renderer);

	/* select state and set value */
	switch(pname)
	{
	case GI_GL_USE_VERTEX_TEXTURE:
	case GI_GL_USE_GEOMETRY_SHADER:
		if(param)
			pRenderer->gim_flags |= pname;
		else
			pRenderer->gim_flags &= ~pname;
		break;
	default:
		GIContext_error(pRenderer->context, GI_INVALID_ENUM);
	}
}

/** Set integer render parameter.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup opengl
 */
void GIAPIENTRY giGLRenderParameteri(GIenum pname, GIint param)
{
	GIRenderer *pRenderer = &(GIContext_current()->renderer);

	/* select state and set value */
	switch(pname)
	{
	case GI_RENDER_RESOLUTION_U:
	case GI_RENDER_RESOLUTION_V:
		if(param != 0 && param < 2)
			GIContext_error(pRenderer->context, GI_INVALID_VALUE);
		else
			pRenderer->render_res[pname-GI_RENDER_RESOLUTION_U] = param;
		break;
	case GI_RENDER_CACHE_SIZE:
		if(param >= 0)
			pRenderer->gim_cache_size = param;
		else
			GIContext_error(pRenderer->context, GI_INVALID_VALUE);
		break;
	default:
		GIContext_error(pRenderer->context, GI_INVALID_ENUM);
	}
}

/** Set integer attribute render parameter.
 *  \param attrib attribute channel
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup opengl
 */
void GIAPIENTRY giGLAttribRenderParameteri(GIuint attrib, GIenum pname, GIint param)
{
	GIRenderer *pRenderer = &(GIContext_current()->renderer);
	if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pRenderer->context, GI_INVALID_VALUE);
		return;
	}

	/* select state and set value */
	switch(pname)
	{
	case GI_GL_RENDER_SEMANTIC:
		if(param == GI_GL_VERTEX || param == GI_GL_NORMAL || 
			param == GI_GL_COLOR || param == GI_GL_SECONDARY_COLOR || 
			param == GI_GL_FOG_COORD || param == GI_GL_EVAL_COORD || 
			param == GI_GL_TEXTURE_COORD || param == GI_GL_VERTEX_ATTRIB || 
			param == GI_NONE)
			pRenderer->attrib_semantic[attrib] = param;
		else
			GIContext_error(pRenderer->context, GI_INVALID_VALUE);
		break;
	case GI_GL_RENDER_CHANNEL:
		if(param >= 0)
			pRenderer->attrib_channel[attrib] = param;
		else
			GIContext_error(pRenderer->context, GI_INVALID_VALUE);
		break;
	case GI_TEXTURE_COORD_DOMAIN:
		if(param == GI_UNIT_SQUARE || param == GI_HALF_TEXEL_INDENT)
			pRenderer->image_domain[attrib] = param;
		else
			GIContext_error(pRenderer->context, GI_INVALID_VALUE);
		break;
	default:
		GIContext_error(pRenderer->context, GI_INVALID_ENUM);
	}
}

/** Render current mesh immediately.
 *  This functions renders the current bound mesh by \c glBegin()/glEnd() 
 *  using the current OpenGL context. This function should be built into a 
 *  display list when needed more than once.
 *  \ingroup opengl
 */
void GIAPIENTRY giGLDrawMesh()
{
#if 0
	GIRenderer *pRenderer = &(GIContext_current()->renderer);
	GIMesh *pMesh = pRenderer->context->mesh;
	GIPatch *pPatch, *pPEnd;
	GIFace *pFace, *pFEnd;
	GIHalfEdge *pHalfEdge;
	GIRenderAttrib arrAttribs[GI_ATTRIB_COUNT];
	GIuint i, j, a, uiNumAttribs;

	/* setup everything */
	if(!GIRenderer_setup_mesh_rendering(pRenderer, 
		arrAttribs, &uiNumAttribs, &pPatch, &pPEnd))
		return;
	if(pPatch)
	{
		pFace = pPatch->faces;
		pFEnd = pPEnd->faces;
	}
	else
		pFace = pFEnd = pMesh->faces;

	/* render faces */
	glBegin(GL_TRIANGLES);
		do
		{
			for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
			{
				for(j=0; j<uiNumAttribs; ++j)
				{
					a = arrAttribs[j].attrib;
					switch(pMesh->asemantic[a])
					{
					case GI_POSITION_ATTRIB:
						if(arrAttribs[j].index < 0)
						{
							//printf("vertex %d ", pHalfEdge->vstart->id);
							(*(GIGLattribdfunc)arrAttribs[j].func)(pHalfEdge->vstart->coords);
						}else
						{
							//printf("vertex %d ", pHalfEdge->vstart->id);
							(*(GIGLidxattribdfunc)arrAttribs[j].func)(
								arrAttribs[j].index, pHalfEdge->vstart->coords);
						}
						break;
					case GI_PARAM_ATTRIB:
						if(arrAttribs[j].index < 0)
						{
							//printf("UV %d ", pHalfEdge->pstart->vertex->id);
							(*(GIGLattribdfunc)arrAttribs[j].func)(pHalfEdge->pstart->params);
						}
						else
						{
							//printf("UV %d ", pHalfEdge->pstart->vertex->id);
							(*(GIGLidxattribdfunc)arrAttribs[j].func)(
								arrAttribs[j].index, pHalfEdge->pstart->params);
						}
						break;
					//case GI_PARAM_STRETCH_ATTRIB:
					//	if(arrAttribs[j].index < 0)
					//		(*(GIGLattribdfunc)arrAttribs[j].func)(&pHalfEdge->pstart->stretch);
					//	else
					//		(*(GIGLidxattribdfunc)arrAttribs[j].func)(
					//			arrAttribs[j].index, &pHalfEdge->pstart->stretch);
					//	break;
					default:
						if(arrAttribs[j].index < 0)
							(*(GIGLattribffunc)arrAttribs[j].func)(
								(GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]));
						else
							(*(GIGLidxattribffunc)arrAttribs[j].func)(arrAttribs[j].index, 
								(GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]));
					}
				}
			}
			pFace = pFace->next;
		}while(pFace != pFEnd);
	glEnd();
#endif
}


/** Render current cut of current mesh immediately.
 *  This functions renders the current cut of the current bound mesh as lines 
 *  by \c glBegin()/glEnd() using the current OpenGL context. This function should 
 *  be built into a display list when needed more than once.
 *  \ingroup opengl
 */
void GIAPIENTRY giGLDrawCut()
{
#if 0
	GIRenderer *pRenderer = &(GIContext_current()->renderer);
	GIMesh *pMesh = pRenderer->context->mesh;
	GIPatch *pPatch, *pPEnd;
	GIHalfEdge *pHalfEdge;
	GIParam *pParam;
	GIRenderAttrib arrAttribs[GI_ATTRIB_COUNT];
	GIuint j, a, uiNumAttribs;

	/* setup everything */
	if(!GIRenderer_setup_mesh_rendering(pRenderer, 
		arrAttribs, &uiNumAttribs, &pPatch, &pPEnd) || !pMesh->patches)
		return;

	/* render half edges */
	do
	{
		pParam = pPatch->params;
		glBegin(GL_LINE_LOOP);
			do
			{
				pHalfEdge = pParam->cut_hedge;
				for(j=0; j<uiNumAttribs; ++j)
				{
					a = arrAttribs[j].attrib;
					switch(pMesh->asemantic[a])
					{
					case GI_POSITION_ATTRIB:
						if(arrAttribs[j].index < 0)
							(*(GIGLattribdfunc)arrAttribs[j].func)(pHalfEdge->vstart->coords);
						else
							(*(GIGLidxattribdfunc)arrAttribs[j].func)(
								arrAttribs[j].index, pHalfEdge->vstart->coords);
						break;
					case GI_PARAM_ATTRIB:
						if(arrAttribs[j].index < 0)
							(*(GIGLattribdfunc)arrAttribs[j].func)(pHalfEdge->pstart->params);
						else
							(*(GIGLidxattribdfunc)arrAttribs[j].func)(
								arrAttribs[j].index, pHalfEdge->pstart->params);
						break;
					case GI_PARAM_STRETCH_ATTRIB:
						if(arrAttribs[j].index < 0)
							(*(GIGLattribdfunc)arrAttribs[j].func)(&pHalfEdge->pstart->stretch);
						else
							(*(GIGLidxattribdfunc)arrAttribs[j].func)(
								arrAttribs[j].index, &pHalfEdge->pstart->stretch);
						break;
					default:
						if(arrAttribs[j].index < 0)
							(*(GIGLattribffunc)arrAttribs[j].func)(
								(GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]));
						else
							(*(GIGLidxattribffunc)arrAttribs[j].func)(arrAttribs[j].index, 
								(GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]));
					}
				}
				pParam = pHalfEdge->next->pstart;
			}while(pParam != pPatch->params);
		glEnd();
		pPatch = pPatch->next;
	}while(pPatch != pPEnd);
#endif
}

/** Draw Geometry Images as triangular mesh.
 *  \ingroup opengl
 */
void GIAPIENTRY giGLDrawGIM()
{
#if 0
	GIRenderer *pRenderer = &(GIContext_current()->renderer);
	GIGLManager *pGL = pRenderer->context->gl_manager;
	GIImage *pImage, *pGIM = NULL;
	GIGLShader *pGIMRenderer;
	GIGLGIMCache *pTCache = NULL, *pICache = NULL;
	GIfloat *pTexCoords = NULL;
	GIuint *pIndices = NULL;
	GIuint uiResU = pRenderer->render_res[0];
	GIuint uiResV = pRenderer->render_res[1];
	GIuint uiWidth, uiHeight, uiPackedWH, uiPackedRes, uiTSize, uiISize;
	GIfloat fTexScale[16] = { 1.0f, 1.0f, 1.0f, 1.0f, 
		1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 
		1.0f, 1.0f, 1.0f, 1.0f };
	GIfloat fTexBias[16] = { 0.0f, 0.0f, 0.0f, 0.0f, 
		0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
		0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
	GIfloat fOffset[16];
	GIfloat fResRatioU, fResRatioV;
	GIuint i, j, uiUnit = 0;
	GLint iPackAlign, iUnpackAlign, iMinFilter, iMagFilter;
	GLint iNormals = 0, iColors = 0, iScaleBiasGeometry = 0, iScaleBiasNormal = 0, iTwoSide = 0;
	GLint iGImage = 0, iNImage = 0, iCImage = 0;
	GLint iTextures[4] = { 0, 0, 0, 0 }, iTImage[4] = { 0, 0, 0, 0 };
	GIboolean bRange, bHalf, bVTF, bGS, bTCreate = GI_TRUE, bICreate = GI_TRUE;
	GIbitfield flags = pRenderer->gim_flags;
	GIboolean bRender[GI_ATTRIB_COUNT];
	GIint iVertex = -1, iColor = -1, iNormal = -1;
	GIint iTexture[4] = { -1, -1, -1, -1 };

	//find attribute images
	memset(bRender, 0, GI_ATTRIB_COUNT*sizeof(GIboolean));
	for(i=0; i<GI_ATTRIB_COUNT; ++i)
	{
		bRender[i] = GI_TRUE;
		switch(pRenderer->attrib_semantic[i])
		{
		case GI_GL_VERTEX:
			if(iVertex == -1)
			{
				pGIM = pRenderer->context->attrib_image[i];
				iVertex = i;
			}
			else
				bRender[i] = GI_FALSE;
			break;
		case GI_GL_NORMAL:
			if(iNormal == -1)
				iNormal = i;
			else
				bRender[i] = GI_FALSE;
			break;
		case GI_GL_COLOR:
			if(iColor == -1)
				iColor = i;
			else
				bRender[i] = GI_FALSE;
			break;
		case GI_GL_TEXTURE_COORD:
			if(pRenderer->attrib_channel[i] < 4)
			{
				if(iTexture[pRenderer->attrib_channel[i]] == -1)
					iTexture[pRenderer->attrib_channel[i]] = i;
				else
					bRender[i] = GI_FALSE;
			}
			else
				bRender[i] = GI_FALSE;
			break;
		default:
			bRender[i] = GI_FALSE;
		}
	}
	if(!pGIM || !(pGIM->data || glIsTexture(pGIM->texture) || pGIM->buffer))
	{
		GIContext_error(pRenderer->context, GI_INVALID_OPERATION);
		return;
	}

	/* error checking */
	if(!pGL->gl_version)
		GIGLManager_init(pGL);
	if(!uiResU)
		uiResU = pGIM->sub_width;
	if(!uiResV)
		uiResV = pGIM->sub_height;
	bRange = pGL->gl_version >= 0x0102;
	bHalf = pGL->half_float_pixel && pGL->half_float_vertex;
	if((flags & GI_GL_USE_VERTEX_TEXTURE) && !pGL->vertex_texture)
	{
		flags &= ~GI_GL_USE_VERTEX_TEXTURE;
		GIContext_error(pRenderer->context, GI_UNSUPPORTED_OPERATION);
	}
	if((flags & GI_GL_USE_GEOMETRY_SHADER) && !pGL->geometry_shader)
	{
		flags &= ~GI_GL_USE_GEOMETRY_SHADER;
		GIContext_error(pRenderer->context, GI_UNSUPPORTED_OPERATION);
	}
	pGIMRenderer = pGL->gim_renderer[flags&0x03];
	if(!pGIMRenderer)
	{
		GIContext_error(pRenderer->context, GI_UNSUPPORTED_OPERATION);
		return;
	}
	bVTF = (flags & GI_GL_USE_VERTEX_TEXTURE);
	bGS = (flags & GI_GL_USE_GEOMETRY_SHADER);
	if(!pGL->glsl || (pGIM->buffer && !pGL->vbo) || (!bVTF && 
		!pGIM->texture && pGIM->type == GI_HALF_FLOAT && 
		!pGL->half_float_vertex) || (!bVTF && (pGIM->sub_height != 
		pGIM->width || pGIM->sub_height != pGIM->height)))
	{
		GIContext_error(pRenderer->context, GI_UNSUPPORTED_OPERATION);
		return;
	}
	if(pGIM->buffer && !pGL->_glIsBuffer(pGIM->buffer))
	{
		GIContext_error(pRenderer->context, GI_INVALID_OPERATION);
		return;
	}

	/* matching resolutions? */
	fResRatioU = (GIfloat)(pGIM->sub_width-1) / (GIfloat)(uiResU-1);
	fResRatioV = (GIfloat)(pGIM->sub_height-1) / (GIfloat)(uiResV-1);
	if(fabs(fResRatioU-floor(fResRatioU)) > 1e-4)
	{
		uiResU = pGIM->sub_width;
		fResRatioU = 1.0f;
		GIContext_error(pRenderer->context, GI_INVALID_VALUE);
	}
	if(fabs(fResRatioV-floor(fResRatioV)) > 1e-4)
	{
		uiResV = pGIM->sub_height;
		fResRatioV = 1.0f;
		GIContext_error(pRenderer->context, GI_INVALID_VALUE);
	}
	if(bVTF)
	{
		uiWidth = bGS ? (uiResU-1) : uiResU;
		uiHeight = bGS ? (uiResV-1) : uiResV;
		fResRatioU = fResRatioV = 1.0f;
	}
	else
	{
		uiWidth = pGIM->sub_width;
		uiHeight = pGIM->sub_height;
	}
	uiPackedWH = uiWidth | (uiHeight<<16);
	uiPackedRes = (bGS ? GI_QUADS : 0) | uiResU | (uiResV<<16);
	uiTSize = ((uiWidth*uiHeight)<<1) * sizeof(GIfloat);
	if(bGS)
		uiISize = ((uiResU-1)*(uiResV-1)) << 2;
	else
#ifdef GI_SINGLE_TRISTRIP
		uiISize = ((uiResU+1)*(uiResV-1)-1) << 1;
#else
		uiISize = (uiResU*(uiResV-1)) << 1;
#endif

	/* collect information about images to be rendered */
	glGetIntegerv(GL_PACK_ALIGNMENT, &iPackAlign);
	glGetIntegerv(GL_UNPACK_ALIGNMENT, &iUnpackAlign);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	for(i=0; i<GI_ATTRIB_COUNT; ++i)
	{
		pImage = pRenderer->context->attrib_image[i];
		if(!pImage)
			bRender[i] = GI_FALSE;
		if(!bRender[i] || (pImage == pGIM && !bVTF))
			continue;

		/* error checking */
		if(!pImage || !(pImage->data || glIsTexture(pImage->texture) || pImage->buffer))
		{
			bRender[i] = GI_FALSE;
			GIContext_error(pRenderer->context, GI_INVALID_OPERATION);
			continue;
		}
		if(uiUnit >= pGL->max_texunits || (pImage->buffer && !pGL->pbo) || 
			(pImage->texture && !pGL->non_power_of_2 && 
			(GI_NON_POWER_OF_2(pImage->width) || GI_NON_POWER_OF_2(pImage->width))) || 
			(pImage->type != GI_UNSIGNED_BYTE && !pGL->texture_float))
		{
			bRender[i] = GI_FALSE;
			GIContext_error(pRenderer->context, GI_UNSUPPORTED_OPERATION);
			continue;
		}
		if(pImage->buffer && !pGL->_glIsBuffer(pImage->buffer))
		{
			bRender[i] = GI_FALSE;
			GIContext_error(pRenderer->context, GI_INVALID_OPERATION);
			continue;
		}

		/* set texture */
		pGL->_glActiveTexture(GL_TEXTURE0+uiUnit);
		if(pImage->texture)
		{
			glBindTexture(GL_TEXTURE_2D, pImage->texture);
			if(pImage == pGIM)
			{
				glGetTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, &iMinFilter);
				glGetTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, &iMagFilter);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			}
		}
		else
		{
			glGenTextures(1, &pImage->texture);
			glBindTexture(GL_TEXTURE_2D, pImage->texture);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
				(pImage==pGIM) ? GL_NEAREST : GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, 
				(pImage==pGIM) ? GL_NEAREST : GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			if(pImage->buffer)
				pGL->_glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pImage->buffer);
			glTexImage2D(GL_TEXTURE_2D, 0, pImage->gl_internal, pImage->width, 
				pImage->height, 0, pImage->gl_format, pImage->type, pImage->data);
			if(pImage->buffer)
				pGL->_glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
		}

		/* get image specific shader parameters */
		if(bVTF && pImage == pGIM)
		{
			iGImage = uiUnit;
			j = 12;
		}
		else if(pRenderer->attrib_semantic[i] == GI_GL_NORMAL)
		{
			iNormals = 1;
			iNImage = uiUnit;
			iScaleBiasNormal = (pImage->type == GI_UNSIGNED_BYTE ? 1 : 0);
			j = 0;
		}
		else if(pRenderer->attrib_semantic[i] == GI_GL_COLOR)
		{
			iColors = 1;
			iCImage = uiUnit;
			j = 2;
		}
		else if(pRenderer->attrib_semantic[i] == GI_GL_TEXTURE_COORD)
		{
			GIuint uiTex = pRenderer->attrib_channel[i];
			iTextures[uiTex] = 1;
			iTImage[uiTex] = uiUnit;
			j = 4 + (uiTex<<1);
		}
		else
			continue;

		/* adjust texCoords */
		if(pRenderer->image_domain[i] == GI_HALF_TEXEL_INDENT)
		{
			fTexScale[j] = (GIfloat)(pImage->sub_width-1) / (GIfloat)pImage->width;
			fTexScale[j+1] = (GIfloat)(pImage->sub_height-1) / (GIfloat)pImage->height;
			fTexBias[j] = ((GIfloat)pImage->offset_x+0.5f) / (GIfloat)pImage->width;
			fTexBias[j+1] = ((GIfloat)pImage->offset_y+0.5f) / (GIfloat)pImage->height;
		}
		if(bGS && bVTF)
		{
			GIfloat f1ResU = fTexScale[j] / (GIfloat)(uiResU-1);
			GIfloat f1ResV = fTexScale[j+1] / (GIfloat)(uiResV-1);
			fOffset[j] = 0.5f * f1ResU;
			fOffset[j+1] = 0.5f * f1ResV;
			fTexScale[j] -= f1ResU;
			fTexScale[j+1] -= f1ResV;
			fTexBias[j] += fOffset[j];
			fTexBias[j+1] += fOffset[j+1];
		}
		++uiUnit;
	}
	pGL->_glActiveTexture(GL_TEXTURE0);

	/* set shader parameters */
	glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE, &iTwoSide);
	pGL->_glUseProgram(pGIMRenderer->program);
	if(bVTF)
	{
		pGL->_glUniform1i((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
			"geometryImage")-1, iGImage);
		if(bGS)
		{
			static const GIfloat fSign[8] = { -1.0f, 1.0f, -1.0f, -1.0f, 
				1.0f, -1.0f, 1.0f, 1.0f };
			pGL->_glUniform4fv((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
				"offset")-1, 4, fOffset);
			pGL->_glUniform2fv((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
				"sign")-1, 4, fSign);
		}
	}
	pGL->_glUniform4fv((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
		"texScale")-1, 4, fTexScale);
	pGL->_glUniform4fv((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
		"texBias")-1, 4, fTexBias);
	pGL->_glUniform4i((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
		"enabledAttribsSBTS")-1, iNormals, iColors, iScaleBiasNormal, iTwoSide);
	pGL->_glUniform4iv((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
		"enabledTextures")-1, 1, iTextures);
	pGL->_glUniform1i((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
		"normalImage")-1, iNImage);
	pGL->_glUniform1i((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
		"colorImage")-1, iCImage);
	pGL->_glUniform1iv((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
		"textureImage")-1, 4, iTImage);
	if(!iColors)
	{
		GLfloat color[4];
		glGetFloatv(GL_CURRENT_COLOR, color);
		pGL->_glUniform4fv((GLint)GIHash_find(&pGIMRenderer->uniform_locs, 
			"baseColor")-1, 1, color);
	}

	/* cache size changed? */
	if(pRenderer->gim_cache_size < pGL->gim_cache_size)
	{
		GIGLGIMCache *pCache, *pCache2;
		GIuint c = 0;
		GI_LIST_FOREACH(pGL->texCoord_cache, pCache)
			if((++c) > pRenderer->gim_cache_size)
				break;
		GI_LIST_NEXT(pGL->texCoord_cache, pCache)
		if(c > pRenderer->gim_cache_size)
		{
			pCache->prev->next = pGL->texCoord_cache;
			pGL->texCoord_cache->prev = pCache->prev;
			do
			{
				pCache2 = pCache;
				pCache = pCache->next;
				GIGLGIMCache_destruct(pCache2, pGL);
				GI_FREE_SINGLE(pCache2, sizeof(GIGLGIMCache));
			}while(pCache != pGL->texCoord_cache);
			if(!pRenderer->gim_cache_size)
				pGL->texCoord_cache = NULL;
		}
		c = 0;
		GI_LIST_FOREACH(pGL->index_cache, pCache)
			if((++c) > pRenderer->gim_cache_size)
				break;
		GI_LIST_NEXT(pGL->index_cache, pCache)
		if(c > pRenderer->gim_cache_size)
		{
			pCache->prev->next = pGL->index_cache;
			pGL->index_cache->prev = pCache->prev;
			do
			{
				pCache2 = pCache;
				pCache = pCache->next;
				GIGLGIMCache_destruct(pCache2, pGL);
				GI_FREE_SINGLE(pCache2, sizeof(GIGLGIMCache));
			}while(pCache != pGL->index_cache);
			if(!pRenderer->gim_cache_size)
				pGL->index_cache = NULL;
		}
	}
	pGL->gim_cache_size = pRenderer->gim_cache_size;

	/* search cache for match or LRU element */
	if(pGL->gim_cache_size && pGL->vbo)
	{
		/* search texCoord cache */
		GIuint c = 0;
		GI_LIST_FOREACH(pGL->texCoord_cache, pTCache)
			if(pTCache->width_height == uiPackedWH)
			{
				bTCreate = GI_FALSE;
				c = pGL->gim_cache_size;
				break;
			}
			++c;
		GI_LIST_NEXT(pGL->texCoord_cache, pTCache)
		if(c < pGL->gim_cache_size)
			pTCache = NULL;
		if(pTCache)
		{
			if(pTCache != pGL->texCoord_cache)
			{
				GI_LIST_REMOVE(pGL->texCoord_cache, pTCache);
				GI_LIST_ADD(pGL->texCoord_cache, pTCache);
				pGL->texCoord_cache = pTCache;
			}
		}
		else
		{
			pTCache = (GIGLGIMCache*)GI_CALLOC_SINGLE(sizeof(GIGLGIMCache));
			GI_LIST_ADD(pGL->texCoord_cache, pTCache);
			pGL->texCoord_cache = pTCache;
		}

		/* bind or create buffer and data if neccessary */
		if(pTCache->buffer)
			pGL->_glBindBuffer(GL_ARRAY_BUFFER, pTCache->buffer);
		else
		{
			pGL->_glGenBuffers(1, &pTCache->buffer);
			pGL->_glBindBuffer(GL_ARRAY_BUFFER, pTCache->buffer);
		}
		if(bTCreate)
		{
			pTCache->width_height = uiPackedWH;
			pGL->_glBufferData(GL_ARRAY_BUFFER, uiTSize, NULL, GL_STATIC_DRAW);
			pTexCoords = (GIfloat*)pGL->_glMapBuffer(
				GL_ARRAY_BUFFER, GL_WRITE_ONLY);
		}

		/* search index cache if needed */
		if(!bVTF || !bGS)
		{
			c = 0;
			GI_LIST_FOREACH(pGL->index_cache, pICache)
				if(pICache->width_height == uiPackedWH && 
					pICache->res_uv == uiPackedRes)
				{
					bICreate = GI_FALSE;
					c = pGL->gim_cache_size;
					break;
				}
				++c;
			GI_LIST_NEXT(pGL->index_cache, pICache)
			if(c < pGL->gim_cache_size)
				pICache = NULL;
			if(pICache)
			{
				if(pICache != pGL->index_cache)
				{
					GI_LIST_REMOVE(pGL->index_cache, pICache);
					GI_LIST_ADD(pGL->index_cache, pICache);
					pGL->index_cache = pICache;
				}
			}
			else
			{
				pICache = (GIGLGIMCache*)GI_CALLOC_SINGLE(sizeof(GIGLGIMCache));
				GI_LIST_ADD(pGL->index_cache, pICache);
				pGL->index_cache = pICache;
			}

			/* bind or create buffer and data if neccessary */
			if(pICache->buffer)
				pGL->_glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pICache->buffer);
			else
			{
				pGL->_glGenBuffers(1, &pICache->buffer);
				pGL->_glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pICache->buffer);
			}
			if(bICreate)
			{
				pICache->width_height = uiPackedWH;
				pICache->res_uv = uiPackedRes;
				pGL->_glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
					uiISize*sizeof(GIuint), NULL, GL_STATIC_DRAW);
				pIndices = (GIuint*)pGL->_glMapBuffer(
					GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
			}
		}
		else
			bICreate = GI_FALSE;
	}

	/* create texCoord data if neccessary */
	if(bTCreate)
	{
		GIuint idx;
		GIfloat fScaleX = 1.0f / (uiWidth-1), fScaleY = 1.0f / (uiHeight-1), fX;

		/* create texCoord data */
		if(!pTexCoords)
			pTexCoords = (GIfloat*)GI_MALLOC_ARRAY(uiTSize, 1);
		for(i=0; i<uiWidth; ++i)
		{
			fX = (GIfloat)i * fScaleX;
			for(j=0; j<uiHeight; ++j)
			{
				idx = (i+j*uiWidth) << 1;
				pTexCoords[idx] = fX;
				pTexCoords[idx+1] = (GIfloat)j * fScaleY;
			}
		}

		/* unmap buffer if neccessary */
		if(pTCache)
		{
			pGL->_glUnmapBuffer(GL_ARRAY_BUFFER);
			pTexCoords = NULL;
		}
	}

	/* set geometry data */
	if(bVTF)
	{
		/* texCoords as vertices */
		glVertexPointer(2, GL_FLOAT, 0, pTexCoords);
		glEnableClientState(GL_VERTEX_ARRAY);
	}
	else
	{
		GIuint uiGSize = pGIM->size;
		if(pGIM->type == GI_HALF_FLOAT && !bHalf)
			uiGSize <<= 1;

		/* set texCoord data */
		pGL->_glClientActiveTexture(GL_TEXTURE0);
		glTexCoordPointer(2, GL_FLOAT, 0, pTexCoords);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);

		/* get geometry data */
		if(pGIM->texture)
		{
			GLint iBoundT;
			glGetIntegerv(GL_TEXTURE_BINDING_2D, &iBoundT);
			glBindTexture(GL_TEXTURE_2D, pGIM->texture);
			if(!pGL->pbo || (pGIM->type == GI_HALF_FLOAT && !bHalf))
			{
				pGIM->data = GI_MALLOC_ARRAY(uiGSize, 1);
				glGetTexImage(GL_TEXTURE_2D, 0, pGIM->gl_format, 
					(pGIM->type==GI_HALF_FLOAT && !bHalf) ? 
					GL_FLOAT : pGIM->type, pGIM->data);
			}
			else
			{
				pGIM->buffer = pGL->gim_buffer;
				pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, pGIM->buffer);
				pGL->_glBufferData(GL_PIXEL_PACK_BUFFER, uiGSize, 0, GL_STATIC_DRAW);
				glGetTexImage(GL_TEXTURE_2D, 0, pGIM->gl_format, pGIM->type, NULL);
				pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
			}
			glBindTexture(GL_TEXTURE_2D, iBoundT);
		}
		else if(pGIM->data && pTCache)
		{
			pGIM->buffer = pGL->gim_buffer;
			pGL->_glBindBuffer(GL_ARRAY_BUFFER, pGIM->buffer);
			pGL->_glBufferData(GL_ARRAY_BUFFER, uiGSize, pGIM->data, GL_STATIC_DRAW);
		}
		if(pGL->vbo)
			pGL->_glBindBuffer(GL_ARRAY_BUFFER, pGIM->buffer);

		/* set arrays */
		pGL->_glVertexAttribPointer(0, pGIM->comp, (pGIM->type==GI_HALF_FLOAT && 
			!bHalf) ? GL_FLOAT : pGIM->type, pGIM->type==GI_UNSIGNED_BYTE ? 
			GL_TRUE : GL_FALSE, 0, pGIM->buffer ? NULL : pGIM->data);
		pGL->_glEnableVertexAttribArray(0);
	}

	/* create index data if neccessary */
	if(bICreate)
	{
		GIuint *pIndex;
		GIuint uiRowStart, idx, uiOffsetU = (GIuint)GI_ROUND(fResRatioU), 
			uiOffsetV = (GIuint)GI_ROUND(fResRatioV)*uiWidth;
		if(!pIndices)
			pIndices = (GIuint*)GI_MALLOC_ARRAY(uiISize, sizeof(GIuint));
		pIndex = pIndices;

		/* create index data */
		if(bGS)
		{
			/* quad set */
			for(i=0,uiRowStart=0; i<uiResV-1; ++i,uiRowStart+=uiOffsetV)
			{
				for(j=0,idx=uiRowStart; j<uiResU-1; ++j,idx+=uiOffsetU)
				{
					*(pIndex++) = idx + uiOffsetV;
					*(pIndex++) = idx;
					*(pIndex++) = idx + uiOffsetU;
					*(pIndex++) = idx + uiOffsetU + uiOffsetV;
				}
			}
		}
		else
		{
			/* triangle strip */
			for(i=0,uiRowStart=0; i<uiResV-1; ++i,uiRowStart+=uiOffsetV)
			{
#ifdef GI_SINGLE_TRISTRIP
				if(i)
				{
					*(pIndex++) = idx - uiOffsetU;
					*(pIndex++) = uiRowStart + uiOffsetV;
				}
#endif
				for(j=0,idx=uiRowStart; j<uiResU; ++j,idx+=uiOffsetU)
				{
					*(pIndex++) = idx + uiOffsetV;
					*(pIndex++) = idx;
				}
			}
		}

		/* unmap buffer if neccessary */
		if(pICache)
		{
			pGL->_glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
			pIndices = NULL;
		}
	}

	/* render arrays */
#ifdef GI_SINGLE_TRISTRIP
	if(bGS && bVTF)
		glDrawArrays(GL_POINTS, 0, (uiResU-1)*(uiResV-1));
	else if(bRange && uiISize <= pGL->max_indices && 
			(uiWidth*uiHeight) <= pGL->max_vertices)
		pGL->_glDrawRangeElements(bGS ? GL_LINES_ADJACENCY_ARB : GL_TRIANGLE_STRIP, 
			0, uiWidth*uiHeight-1, uiISize, GL_UNSIGNED_INT, pIndices);
	else
		glDrawElements(bGS ? GL_LINES_ADJACENCY_ARB : GL_TRIANGLE_STRIP, 
			uiISize, GL_UNSIGNED_INT, pIndices);
#else
	if(bGS)
	{
		if(bVTF)
			glDrawArrays(GL_POINTS, 0, (uiResU-1)*(uiResV-1));
		else if(bRange && uiISize <= pGL->max_indices && 
				(uiWidth*uiHeight) <= pGL->max_vertices)
			pGL->_glDrawRangeElements(GL_LINES_ADJACENCY_ARB, 0, 
				uiWidth*uiHeight-1, uiISize, GL_UNSIGNED_INT, pIndices);
		else
			glDrawElements(GL_LINES_ADJACENCY_ARB, 
				uiISize, GL_UNSIGNED_INT, pIndices);
	}
	else
	{
		GLsizei iStripSize = uiResU << 1;
		if(bRange && iStripSize <= pGL->max_indices && iStripSize <= pGL->max_vertices)
			for(i=0,j=0; i<uiResV-1; ++i,j+=uiResU)
				pGL->_glDrawRangeElements(GL_TRIANGLE_STRIP, j, j+iStripSize-1, 
					iStripSize, GL_UNSIGNED_INT, pIndices+(j<<1));
		else
			for(i=0,j=0; i<uiResV-1; ++i,j+=iStripSize)
				glDrawElements(GL_TRIANGLE_STRIP, iStripSize, 
					GL_UNSIGNED_INT, pIndices+j);
	}
#endif

	/* restore GL state */
	if(bVTF)
		glDisableClientState(GL_VERTEX_ARRAY);
	else
	{
		pGL->_glDisableVertexAttribArray(0);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	}
	if(pGIM->type == GI_UNSIGNED_BYTE)
		glPopMatrix();
	if(pGL->vbo)
	{
		pGL->_glBindBuffer(GL_ARRAY_BUFFER, 0);
		if(pICache)
			pGL->_glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		if(pGL->pbo)
			pGL->_glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
	}
	pGL->_glUseProgram(0);
	glPixelStorei(GL_PACK_ALIGNMENT, iPackAlign);
	glPixelStorei(GL_UNPACK_ALIGNMENT, iUnpackAlign);
	pGL->_glActiveTexture(GL_TEXTURE0);

	/* delete temporary data */
	if(!bVTF)
	{
		if(pGIM->buffer && (pGIM->texture || pGIM->data))
			pGIM->buffer = 0;
		else if(pGIM->texture)
		{
			GI_FREE_ARRAY(pGIM->data);
			pGIM->data = 0;
		}
	}
	if(pTexCoords)
		GI_FREE_ARRAY(pTexCoords);
	if(pIndices)
		GI_FREE_ARRAY(pIndices);

	/* delete textures if neccessary */
	for(i=0; i<GI_ATTRIB_COUNT; ++i)
	{
		pImage = pRenderer->context->attrib_image[i];
		if(!bRender[i] || (pImage == pGIM && !bVTF))
			continue;
		if(pImage->texture)
		{
			if(pImage->data || pImage->buffer)
			{
				glDeleteTextures(1, &pImage->texture);
				pImage->texture = 0;
			}
			else if(pImage == pGIM)
			{
				glBindTexture(GL_TEXTURE_2D, pImage->texture);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, iMinFilter);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, iMagFilter);
			}
		}
	}
#endif
}

/** Clean up OpenGL context specific data.
 *  This function should be called before the GL context used by OpenGI is destroyed 
 *  or changes during the lifetime of the active OpenGI context, as OpenGI may create 
 *  GL context local data as e.g. shaders. It need not be called when the GI context 
 *  is destroyed before the GL context, but it has to be called between two calls to 
 *  giSample() or giGLDrawGIM() within two different GL contexts.
 *  \ingroup opengl
 */
void GIAPIENTRY giGLCleanUp()
{
	/* just call gl manager's destructor */
	GIGLManager_destruct(GIContext_current()->gl_manager);
}

/** \internal
 *  \brief Renderer constructor.
 *  \param renderer renderer to construct
 *  \param context context to construct in
 *  \ingroup opengl
 */
void GIRenderer_construct(GIRenderer *renderer, struct _GIContext *context)
{
	GIuint a;

	/* initialize state */
	renderer->context = context;
	renderer->render_res[0] = renderer->render_res[1] = 0;
	renderer->gim_cache_size = 8;
	renderer->gim_flags = 0;
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
	{
		renderer->attrib_semantic[a] = GI_NONE;
		renderer->attrib_channel[a] = 0;
		renderer->image_domain[a] = GI_HALF_TEXEL_INDENT;
	}
	renderer->attrib_semantic[0] = GI_GL_VERTEX;
}

/** \internal
 *  \brief Gather information for mesh rendering.
 *  \param renderer renderer to work on
 *  \param attribs array to take attribute information
 *  \param num_attribs address to store number of rendered attribs at
 *  \param pstart address to store starting patch at
 *  \param pend address to store end patch at
 *  \retval GI_TRUE if set up successfully
 *  \retval GI_FALSE on error
 *  \ingroup opengl
 */
GIboolean GIRenderer_setup_mesh_rendering(GIRenderer *renderer, 
										  GIRenderAttrib *attribs, 
										  GIuint *num_attribs, 
										  GIPatch **pstart, GIPatch **pend)
{
	GIGLManager *pGL = renderer->context->gl_manager;
	GIMesh *pMesh = renderer->context->mesh;
	GIint iSemanticAttrib[GI_GL_SEMANTIC_COUNT];
	GIuint i, a;
	GIboolean bFloat, bParams, bStretch;

	/* error checking */
	if(!pMesh)
	{
		GIContext_error(renderer->context, GI_INVALID_OPERATION);
		return GI_FALSE;
	}
#if 0
	if(!pGL->gl_version)
		GIGLManager_init(pGL);
#endif
	if(pMesh->active_patch)
	{
		*pstart = pMesh->active_patch;
		*pend = pMesh->active_patch->next;
		bParams = pMesh->param_patches && pMesh->active_patch->parameterized;
		bStretch = pMesh->active_patch->param_metric != 0;
	}
	else
	{
		*pstart = *pend = pMesh->patches;
		bParams = pMesh->param_patches > 0;
		bStretch = pMesh->param_metric != 0;
	}

	/* gather attribute information */
	for(i=0; i<GI_GL_SEMANTIC_COUNT; ++i)
		iSemanticAttrib[i] = (i==(GI_GL_TEXTURE_COORD-GI_GL_SEMANTIC_BASE)||
			i==(GI_GL_VERTEX_ATTRIB-GI_GL_SEMANTIC_BASE)) ? 0 : -1;
	for(a=0,i=0; a<GI_ATTRIB_COUNT; ++a)
	{
		GIuint uiSemantic = renderer->attrib_semantic[a];
		if(uiSemantic != GI_NONE)
		{
			/* error checking */
			GIboolean bUsed = (uiSemantic==GI_GL_TEXTURE_COORD||
				uiSemantic==GI_GL_VERTEX_ATTRIB) ? (iSemanticAttrib[uiSemantic-
				GI_GL_SEMANTIC_BASE]&(1<<renderer->attrib_channel[a])) : 
				(iSemanticAttrib[uiSemantic-GI_GL_SEMANTIC_BASE]>=0);
			if(bUsed || (pMesh->asemantic[a] == GI_PARAM_ATTRIB && !bParams) || 
				(pMesh->asemantic[a] == GI_PARAM_STRETCH_ATTRIB && !bStretch) || 
				(pMesh->asemantic[a] == GI_NONE && pMesh->aoffset[a]<0))
			{
				continue;
			}
			if((uiSemantic==GI_GL_VERTEX && pMesh->asize[a]<2) || 
				((uiSemantic==GI_GL_NORMAL || uiSemantic==GI_GL_COLOR || 
				uiSemantic==GI_GL_SECONDARY_COLOR) && pMesh->asize[a]<3) || 
				((uiSemantic==GI_GL_SECONDARY_COLOR || 
				uiSemantic==GI_GL_FOG_COORD) && pGL->gl_version<0x104) || 
				(uiSemantic==GI_GL_TEXTURE_COORD && renderer->attrib_channel[a]!=0 && 
				!pGL->multitexture) || (uiSemantic==GI_GL_VERTEX_ATTRIB && !pGL->glsl))
			{
				GIContext_error(renderer->context, GI_UNSUPPORTED_OPERATION);
				continue;
			}

			/* set GL attribute functions */
			attribs[i].attrib = a;
			attribs[i].index = -1;

			bFloat = (pMesh->asemantic[a]==GI_NONE) ? GI_TRUE : GI_FALSE;
			switch(uiSemantic)
			{
#if 10
			case GI_GL_VERTEX:
				//if(pMesh->asize[a] == 4)
				//	attribs[i].func = bFloat ? (GIvoid*)&glVertex4fv : (GIvoid*)&glVertex4dv;
				//else if(pMesh->asize[a] == 3)
				//	attribs[i].func = bFloat ? (GIvoid*)my_glVertex3fv : (GIvoid*)my_glVertex3dv;
				//else
				//	attribs[i].func = bFloat ? (GIvoid*)&glVertex2fv : (GIvoid*)&glVertex2dv;
				iSemanticAttrib[GI_GL_VERTEX_ATTRIB-GI_GL_SEMANTIC_BASE] |= 1;
				break;
			case GI_GL_NORMAL:
				//attribs[i].func = bFloat ? (GIvoid*)my_glNormal3fv : (GIvoid*)my_glNormal3dv;
				break;
			case GI_GL_COLOR:
				//if(pMesh->asize[a] == 4)
				//	attribs[i].func = bFloat ? (GIvoid*)&glColor4fv : (GIvoid*)&glColor4dv;
				//else
				//	attribs[i].func = bFloat ? (GIvoid*)&glColor3fv : (GIvoid*)&glColor3dv;
				break;
			case GI_GL_SECONDARY_COLOR:
				//attribs[i].func = bFloat ? (GIvoid*)pGL->_glSecondaryColor3fv : (GIvoid*)pGL->_glSecondaryColor3dv;
				break;
			case GI_GL_FOG_COORD:
				//attribs[i].func = bFloat ? (GIvoid*)pGL->_glFogCoordfv : (GIvoid*)pGL->_glFogCoorddv;
				break;
			case GI_GL_EVAL_COORD:
				//if(pMesh->asize[a] > 1)
				//	attribs[i].func = bFloat ? (GIvoid*)&glEvalCoord2fv : (GIvoid*)&glEvalCoord2dv;
				//else
				//	attribs[i].func = bFloat ? (GIvoid*)&glEvalCoord1fv : (GIvoid*)&glEvalCoord1dv;
				break;
			case GI_GL_TEXTURE_COORD:
				if(pGL->multitexture)
				{
					//if(pMesh->asize[a] == 4)
					//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glMultiTexCoord4fv : (GIvoid*)pGL->_glMultiTexCoord4dv;
					//else if(pMesh->asize[a] == 3)
					//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glMultiTexCoord3fv : (GIvoid*)pGL->_glMultiTexCoord3dv;
					//else if(pMesh->asize[a] == 2)
					//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glMultiTexCoord2fv : (GIvoid*)pGL->_glMultiTexCoord2dv;
					//else
					//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glMultiTexCoord1fv : (GIvoid*)pGL->_glMultiTexCoord1dv;
					attribs[i].index = renderer->attrib_channel[a] + GL_TEXTURE0;
				}
				else
				{
					//if(pMesh->asize[a] == 4)
					//	attribs[i].func = bFloat ? (GIvoid*)&glTexCoord4fv : (GIvoid*)&glTexCoord4dv;
					//else if(pMesh->asize[a] == 3)
					//	attribs[i].func = bFloat ? (GIvoid*)&glTexCoord3fv : (GIvoid*)&glTexCoord3dv;
					//else if(pMesh->asize[a] == 2)
					//	attribs[i].func = bFloat ? (GIvoid*)&glTexCoord2fv : (GIvoid*)&glTexCoord2dv;
					//else
					//	attribs[i].func = bFloat ? (GIvoid*)&glTexCoord1fv : (GIvoid*)&glTexCoord1dv;
				}
				break;
#endif
			case GI_GL_VERTEX_ATTRIB:
				//if(pMesh->asize[a] == 4)
				//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glVertexAttrib4fv : (GIvoid*)pGL->_glVertexAttrib4dv;
				//else if(pMesh->asize[a] == 3)
				//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glVertexAttrib3fv : (GIvoid*)pGL->_glVertexAttrib3dv;
				//else if(pMesh->asize[a] == 2)
				//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glVertexAttrib2fv : (GIvoid*)pGL->_glVertexAttrib2dv;
				//else
				//	attribs[i].func = bFloat ? (GIvoid*)pGL->_glVertexAttrib1fv : (GIvoid*)pGL->_glVertexAttrib1dv;

				attribs[i].index = renderer->attrib_channel[a];

				if(renderer->attrib_channel[a] == 0)
					iSemanticAttrib[GI_GL_VERTEX-GI_GL_SEMANTIC_BASE] = i;
			}
			if(uiSemantic == GI_GL_TEXTURE_COORD || uiSemantic == GI_GL_VERTEX_ATTRIB)
				iSemanticAttrib[uiSemantic-GI_GL_SEMANTIC_BASE] |= 
					1 << renderer->attrib_channel[a];
			else
				iSemanticAttrib[uiSemantic-GI_GL_SEMANTIC_BASE] = i;
			++i;
		}
	}

	/* move glVertex to end */
	if(!(*num_attribs = i) || iSemanticAttrib[GI_GL_VERTEX-GI_GL_SEMANTIC_BASE] < 0)
		return GI_FALSE;
	i = iSemanticAttrib[GI_GL_VERTEX-GI_GL_SEMANTIC_BASE];
	if(i != *num_attribs-1)
	{
		GIRenderAttrib temp;
		GI_SWAP(attribs[i], attribs[*num_attribs-1], temp);
	}
	return GI_TRUE;
}

/** \internal
 *  \brief Shader constructor.
 *  \param shader shader to construct
 *  \param gl valid GL manager
 *  \param vsource source code for vertex shader or NULL if not used
 *  \param gsource source code for geometry shader or NULL if not used
 *  \param fsource source code for fragment shader or NULL if not used
 *  \param g_in input primitive type for geometry shader
 *  \param g_out output primitive type for geometry shader
 *  \param num_g_out number of vertices output by geometry shader
 *  \param num_defines number of preprocessor defines
 *  \param ... preprocessor defines as strings
 *  \ingroup opengl
 */
void GIGLShader_construct(GIGLShader *shader, GIGLManager *gl, 
						  const GIchar *vsource, const GIchar *gsource, 
						  const GIchar *fsource, GIenum g_in, GIenum g_out, 
						  GIuint num_g_out, GIuint num_defines, ...)
{
#if 0
	GIchar *szDefines;
	GLuint uiFailed = 0;
	GLint iStatus;
	GLint i, iNum, iMaxLen, iLength, iSize, iNumStrings = (num_defines ? 2 : 1);
	GLenum uiType;
	GLchar *szName, *szSource[2];

	/* clear shader */
	GIGLShader_destruct(shader, gl);

	/* fetch defines */
	if(num_defines)
	{
		GIchar szDefine[64];
		va_list args;
		va_start(args, num_defines);
		szDefines = (GIchar*)GI_MALLOC_ARRAY(num_defines, 64*sizeof(GIchar));
		szDefines[0] = 0;
		for(i=0; i<num_defines; ++i)
		{
			sprintf(szDefine, "#define %s 1 \n", va_arg(args, GIchar*));
			strcat(szDefines, szDefine);
		}
		va_end(args);
		szSource[0] = szDefines;
	}

	/* create and compile shaders */
	if(vsource)
	{
		szSource[iNumStrings-1] = vsource;
		shader->v_shader = gl->_glCreateShader(GL_VERTEX_SHADER);
		gl->_glShaderSource(shader->v_shader, iNumStrings, szSource, NULL);
		gl->_glCompileShader(shader->v_shader);
		gl->_glGetShaderiv(shader->v_shader, GL_COMPILE_STATUS, &iStatus);
		if(!iStatus)
			uiFailed = shader->v_shader;
	}
	if(gsource && gl->geometry_shader)
	{
		szSource[iNumStrings-1] = gsource;
		//shader->g_shader = gl->_glCreateShader(GL_GEOMETRY_SHADER_ARB);
		//gl->_glShaderSource(shader->g_shader, iNumStrings, szSource, NULL);
		//gl->_glCompileShader(shader->g_shader);
		//gl->_glGetShaderiv(shader->g_shader, GL_COMPILE_STATUS, &iStatus);
		if(!iStatus)
			uiFailed = shader->g_shader;
	}
	if(fsource)
	{
		szSource[iNumStrings-1] = fsource;
		//shader->f_shader = gl->_glCreateShader(GL_FRAGMENT_SHADER);
		//gl->_glShaderSource(shader->f_shader, iNumStrings, szSource, NULL);
		//gl->_glCompileShader(shader->f_shader);
		//gl->_glGetShaderiv(shader->f_shader, GL_COMPILE_STATUS, &iStatus);
		if(!iStatus)
			uiFailed = shader->f_shader;
	}
	if(num_defines)
		GI_FREE_ARRAY(szDefines);

	/* error occured? */
	if(uiFailed)
	{
		/* print shader info log */
		GLint iLength;
		gl->_glGetShaderiv(uiFailed, GL_INFO_LOG_LENGTH, &iLength);
		if(iLength > 1)
		{
			char *pInfoLog = (char*)GI_MALLOC_ARRAY(iLength, sizeof(char));
			gl->_glGetShaderInfoLog(uiFailed, iLength, NULL, pInfoLog);
			GIDebug(printf("OpenGI could not compile shader:\n%s\n", pInfoLog));
			GI_FREE_ARRAY(pInfoLog);
		}
		if(shader->v_shader)
			gl->_glDeleteShader(shader->v_shader);
		if(shader->g_shader)
			gl->_glDeleteShader(shader->g_shader);
		if(shader->f_shader)
			gl->_glDeleteShader(shader->f_shader);
		return;
	}

	/* create and link program */
	shader->program = gl->_glCreateProgram();
	if(shader->v_shader)
		gl->_glAttachShader(shader->program, shader->v_shader);
	if(shader->g_shader)
	{
		gl->_glAttachShader(shader->program, shader->g_shader);
		gl->_glProgramParameteri(shader->program, GL_GEOMETRY_INPUT_TYPE_EXT, g_in);
		gl->_glProgramParameteri(shader->program, GL_GEOMETRY_OUTPUT_TYPE_EXT, g_out);
		gl->_glProgramParameteri(shader->program, GL_GEOMETRY_VERTICES_OUT_EXT, num_g_out);
	}
	if(shader->f_shader)
		gl->_glAttachShader(shader->program, shader->f_shader);
	gl->_glLinkProgram(shader->program);

	/* error occured? */
	gl->_glGetProgramiv(shader->program, GL_LINK_STATUS, &iStatus);
	if(!iStatus)
	{
		/* print program info log */
		GLint iLength;
		gl->_glGetProgramiv(shader->program, GL_INFO_LOG_LENGTH, &iLength);
		if(iLength > 1)
		{
			char *pInfoLog = (char*)GI_MALLOC_ARRAY(iLength, sizeof(char));
			gl->_glGetProgramInfoLog(shader->program, iLength, NULL, pInfoLog);
			GIDebug(printf("OpenGI could not link shader program:\n%s\n", pInfoLog));
			GI_FREE_ARRAY(pInfoLog);
		}
		gl->_glDeleteProgram(shader->program);
		if(shader->v_shader)
			gl->_glDeleteShader(shader->v_shader);
		if(shader->g_shader)
			gl->_glDeleteShader(shader->g_shader);
		if(shader->f_shader)
			gl->_glDeleteShader(shader->f_shader);
		return;
	}

	/* get attribute locations */
	GIHash_construct(&shader->attrib_locs, 32, 0.0f, 32*sizeof(GIchar), hash_string, compare_string, copy_string);
	gl->_glGetProgramiv(shader->program, GL_ACTIVE_ATTRIBUTES, &iNum);
	if(iNum > 0)
	{
		gl->_glGetProgramiv(shader->program, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH, &iMaxLen);
		szName = (GIchar*)GI_MALLOC_ARRAY(iMaxLen, sizeof(GIchar));
		for(i=0; i<iNum; ++i)
		{
			gl->_glGetActiveAttrib(shader->program, i, iMaxLen, &iLength, &iSize, &uiType, szName);
			if(strncmp(szName, "gl_", 3))
				GIHash_insert(&shader->attrib_locs, szName, 
					(GIvoid*)(gl->_glGetAttribLocation(shader->program, szName)+1));
//			printf("attribute: %s : %d\n", szName, gl->_glGetAttribLocation(shader->program, szName));
		}
		GI_FREE_ARRAY(szName);
	}

	/* get uniform locations */
	GIHash_construct(&shader->uniform_locs, 32, 0.0f, 32*sizeof(GIchar), hash_string, compare_string, copy_string);
	gl->_glGetProgramiv(shader->program, GL_ACTIVE_UNIFORMS, &iNum);
	if(iNum > 0)
	{
		gl->_glGetProgramiv(shader->program, GL_ACTIVE_UNIFORM_MAX_LENGTH, &iMaxLen);
		szName = (GIchar*)GI_MALLOC_ARRAY(iMaxLen, sizeof(GIchar));
		for(i=0; i<iNum; ++i)
		{
			gl->_glGetActiveUniform(shader->program, i, iMaxLen, &iLength, &iSize, &uiType, szName);
			if(strncmp(szName, "gl_", 3))
			{
				szName[strcspn(szName, "[")] = 0;
				GIHash_insert(&shader->uniform_locs, szName, 
					(GIvoid*)(gl->_glGetUniformLocation(shader->program, szName)+1));
			}
//			printf("uniform: %s : %d\n", szName, gl->_glGetUniformLocation(shader->program, szName));
		}
		GI_FREE_ARRAY(szName);
	}
#endif
}

/** \internal
 *  \brief Shader destructor.
 *  \param shader shader to destruct
 *  \param gl valid GL manager
 *  \ingroup opengl
 */
void GIGLShader_destruct(GIGLShader *shader, GIGLManager *gl)
{
#if 0
	/* delete GL objects */
	if(shader->program)
	{
		gl->_glDeleteProgram(shader->program);
		shader->program = 0;
		if(shader->v_shader)
		{
			gl->_glDeleteShader(shader->v_shader);
			shader->v_shader = 0;
		}
		if(shader->g_shader)
		{
			gl->_glDeleteShader(shader->g_shader);
			shader->g_shader = 0;
		}
		if(shader->f_shader)
		{
			gl->_glDeleteShader(shader->f_shader);
			shader->f_shader = 0;
		}
	}
#endif
	/* clean up tables */
	GIHash_destruct(&shader->attrib_locs, 0);
	GIHash_destruct(&shader->uniform_locs, 0);
}

/** \internal
 *  \brief GIM cache destructor.
 *  \param cache cache to destruct
 *  \param gl OpenGL manager
 *  \ingroup opengl
 */
void GIGLGIMCache_destruct(GIGLGIMCache *cache, GIGLManager *gl)
{
	if(cache->buffer)
	{
		//gl->_glDeleteBuffers(1, &cache->buffer);
		cache->width_height = cache->res_uv = 0;
	}
}

/** \internal
 *  \brief Initialize OpenGL manager.
 *  \param mgr OpenGL manager to initialize
 *  \ingroup opengl
 */
void GIGLManager_init(GIGLManager *mgr)
{
#if 0
	const unsigned char* szExtensions = glGetString(GL_EXTENSIONS);
	GIuint i, uiMajor, uiMinor;
	GIboolean v12, v13, v14, v15, v20, v21, v30;

	/* check version and extensions */
	memset(mgr, 0, sizeof(GIGLManager));
	sscanf(glGetString(GL_VERSION), "%d.%d", &uiMajor, &uiMinor);
	mgr->gl_version = (uiMajor<<8) | uiMinor;
	v12 = mgr->gl_version >= 0x102;
	v13 = mgr->gl_version >= 0x103;
	v14 = mgr->gl_version >= 0x104;
	v15 = mgr->gl_version >= 0x105;
	v20 = mgr->gl_version >= 0x200;
	v21 = mgr->gl_version >= 0x201;
	v30 = mgr->gl_version >= 0x300;
	mgr->texture_3d = (mgr->gl_version >= 0x102);
	mgr->texture_cube_map = (v13 || strstr(szExtensions, "GL_ARB_texture_cube_map") != NULL);
	mgr->non_power_of_2 = (v20 || strstr(szExtensions, "GL_ARB_texture_non_power_of_two") != NULL);

	/* check if range elements supported */
	if(v12)
	{
		mgr->_glDrawRangeElements = (PFNGLDRAWRANGEELEMENTSPROC)
			GI_GL_PROC_ADDRESS("glDrawRangeElements");
		glGetIntegerv(GL_MAX_ELEMENTS_VERTICES, &mgr->max_vertices);
		glGetIntegerv(GL_MAX_ELEMENTS_INDICES, &mgr->max_indices);
	}

	/* check if multitexturing supported */
	mgr->multitexture = (v13 || strstr(szExtensions, "GL_ARB_multitexture") != NULL);
	if(mgr->multitexture)
	{
		/* get function pointers */
		mgr->_glActiveTexture = (PFNGLACTIVETEXTUREPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glActiveTexture" : "glActiveTextureARB");
		mgr->_glClientActiveTexture = (PFNGLCLIENTACTIVETEXTUREPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glClientActiveTexture" : "glClientActiveTextureARB");
		mgr->_glMultiTexCoord1i = (PFNGLMULTITEXCOORD1IPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord1i" : "glMultiTexCoord1iARB");
		mgr->_glMultiTexCoord1f = (PFNGLMULTITEXCOORD1FPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord1f" : "glMultiTexCoord1fARB");
		mgr->_glMultiTexCoord1fv = (PFNGLMULTITEXCOORD1FVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord1fv" : "glMultiTexCoord1fvARB");

#if 0
		mgr->_glMultiTexCoord2fv = (PFNGLMULTITEXCOORD2FVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord2fv" : "glMultiTexCoord2fvARB");
#else
		mgr->_glMultiTexCoord2fv = (PFNGLMULTITEXCOORD2FVPROC)my_glMultiTexCoord2fv;
#endif		
		mgr->_glMultiTexCoord3fv = (PFNGLMULTITEXCOORD3FVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord3fv" : "glMultiTexCoord3fvARB");
		mgr->_glMultiTexCoord4fv = (PFNGLMULTITEXCOORD4FVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord4fv" : "glMultiTexCoord4fvARB");
		mgr->_glMultiTexCoord1d = (PFNGLMULTITEXCOORD1DPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord1d" : "glMultiTexCoord1dARB");
		mgr->_glMultiTexCoord1dv = (PFNGLMULTITEXCOORD1DVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord1dv" : "glMultiTexCoord1dvARB");

#if 0
		mgr->_glMultiTexCoord2dv = (PFNGLMULTITEXCOORD2DVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord2dv" : "glMultiTexCoord2dvARB");
#else
		mgr->_glMultiTexCoord2dv = (PFNGLMULTITEXCOORD2DVPROC)my_glMultiTexCoord2dv;
#endif
		mgr->_glMultiTexCoord3dv = (PFNGLMULTITEXCOORD3DVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord3dv" : "glMultiTexCoord3dvARB");
		mgr->_glMultiTexCoord4dv = (PFNGLMULTITEXCOORD4DVPROC)
			GI_GL_PROC_ADDRESS(v13 ? "glMultiTexCoord4dv" : "glMultiTexCoord4dvARB");
		glGetIntegerv(GL_MAX_TEXTURE_UNITS, &mgr->max_texunits);
	}
	else
		mgr->max_texunits = 1;

	/* check for 1.4 */
	if(v14)
	{
		mgr->_glFogCoordfv = (PFNGLFOGCOORDFVPROC)GI_GL_PROC_ADDRESS("glFogCoordfv");
		mgr->_glFogCoorddv = (PFNGLFOGCOORDDVPROC)GI_GL_PROC_ADDRESS("glFogCoorddv");
		mgr->_glSecondaryColor3fv = (PFNGLSECONDARYCOLOR3FVPROC)
			GI_GL_PROC_ADDRESS("glSecondaryColor3fv");
		mgr->_glSecondaryColor3dv = (PFNGLSECONDARYCOLOR3DVPROC)
			GI_GL_PROC_ADDRESS("glSecondaryColor3dv");
	}

	/* check if VBOs supported */
	mgr->vbo = (v15 || strstr(szExtensions, "GL_ARB_vertex_buffer_object") != NULL);
	mgr->pbo = (v21 || strstr(szExtensions, "GL_ARB_pixel_buffer_object") != NULL);
	if(mgr->vbo)
	{
		mgr->_glBindBuffer = (PFNGLBINDBUFFERPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glBindBuffer" : "glBindBufferARB");
		mgr->_glBufferData = (PFNGLBUFFERDATAPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glBufferData" : "glBufferDataARB");
		mgr->_glBufferSubData = (PFNGLBUFFERSUBDATAPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glBufferSubData" : "glBufferSubDataARB");
		mgr->_glDeleteBuffers = (PFNGLDELETEBUFFERSPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glDeleteBuffers" : "glDeleteBuffersARB");
		mgr->_glGenBuffers = (PFNGLGENBUFFERSPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glGenBuffers" : "glGenBuffersARB");
		mgr->_glGetBufferSubData = (PFNGLGETBUFFERSUBDATAPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glGetBufferSubData" : "glGetBufferSubDataARB");
		mgr->_glIsBuffer = (PFNGLISBUFFERPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glIsBuffer" : "glIsBufferARB");
		mgr->_glMapBuffer = (PFNGLMAPBUFFERPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glMapBuffer" : "glMapBufferARB");
		mgr->_glUnmapBuffer = (PFNGLUNMAPBUFFERPROC)
			GI_GL_PROC_ADDRESS(v15 ? "glUnmapBuffer" : "glUnmapBufferARB");
		mgr->_glGenBuffers(1, &mgr->gim_buffer);
	}

	/* check if shader supported */
	mgr->glsl = (mgr->max_texunits >= 4 && (v20 || strstr(szExtensions, "GL_ARB_shading_language_100") != NULL));
	if(mgr->glsl)
	{
		static const char *szDefine[GI_NUM_SAMPLERS] = { 
			"DEFAULT", "TEXTURE1D", "TEXTURE2D", "TEXTURE3D", 
			"TEXTURECUBE", "NORMALIZED", "PACKED" };
		GLint iVTUnits;

		/* get function pointers */
		mgr->_glAttachShader = (PFNGLATTACHSHADERPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glAttachShader" : "glAttachObjectARB");
		mgr->_glCompileShader = (PFNGLCOMPILESHADERPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glCompileShader" : "glCompileShaderARB");
		mgr->_glCreateProgram = (PFNGLCREATEPROGRAMPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glCreateProgram" : "glCreateProgramObjectARB");
		mgr->_glCreateShader = (PFNGLCREATESHADERPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glCreateShader" : "glCreateShaderObjectARB");
		mgr->_glDeleteProgram = (PFNGLDELETEPROGRAMPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glDeleteProgram" : "glDeleteObjectARB");
		mgr->_glDeleteShader = (PFNGLDELETESHADERPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glDeleteShader" : "glDeleteObjectARB");
		mgr->_glDetachShader = (PFNGLDETACHSHADERPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glDetachShader" : "glDetachObjectARB");
		mgr->_glDisableVertexAttribArray = (PFNGLDISABLEVERTEXATTRIBARRAYPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glDisableVertexAttribArray" : "glDisableVertexAttribArrayARB");
		mgr->_glEnableVertexAttribArray = (PFNGLENABLEVERTEXATTRIBARRAYPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glEnableVertexAttribArray" : "glEnableVertexAttribArrayARB");
		mgr->_glGetActiveAttrib = (PFNGLGETACTIVEATTRIBPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetActiveAttrib" : "glGetActiveAttribARB");
		mgr->_glGetActiveUniform = (PFNGLGETACTIVEUNIFORMPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetActiveUniform" : "glGetActiveUniformARB");
		mgr->_glGetAttribLocation = (PFNGLGETATTRIBLOCATIONPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetAttribLocation" : "glGetAttribLocationARB");
		mgr->_glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetProgramInfoLog" : "glGetInfoLogARB");
		mgr->_glGetProgramiv = (PFNGLGETPROGRAMIVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetProgramiv" : "glGetObjectParameterivARB");
		mgr->_glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetShaderInfoLog" : "glGetInfoLogARB");
		mgr->_glGetShaderiv = (PFNGLGETSHADERIVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetShaderiv" : "glGetObjectParameterivARB");
		mgr->_glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glGetUniformLocation" : "glGetUniformLocationARB");
		mgr->_glLinkProgram = (PFNGLLINKPROGRAMPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glLinkProgram" : "glLinkProgramARB");
		mgr->_glShaderSource = (PFNGLSHADERSOURCEPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glShaderSource" : "glShaderSourceARB");
		mgr->_glUniform1i = (PFNGLUNIFORM1IPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform1i" : "glUniform1iARB");
		mgr->_glUniform2i = (PFNGLUNIFORM2IPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform2i" : "glUniform2iARB");
		mgr->_glUniform3i = (PFNGLUNIFORM3IPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform3i" : "glUniform3iARB");
		mgr->_glUniform4i = (PFNGLUNIFORM4IPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform4i" : "glUniform4iARB");
		mgr->_glUniform1iv = (PFNGLUNIFORM1IVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform1iv" : "glUniform1ivARB");
		mgr->_glUniform4iv = (PFNGLUNIFORM4IVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform4iv" : "glUniform4ivARB");
		mgr->_glUniform2f = (PFNGLUNIFORM2FPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform2f" : "glUniform2fARB");
		mgr->_glUniform2fv = (PFNGLUNIFORM2FVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform2fv" : "glUniform2fvARB");
		mgr->_glUniform4fv = (PFNGLUNIFORM4FVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUniform4fv" : "glUniform4fvARB");
		mgr->_glUseProgram = (PFNGLUSEPROGRAMPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glUseProgram" : "glUseProgramObjectARB");
		mgr->_glValidateProgram = (PFNGLVALIDATEPROGRAMPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glValidateProgram" : "glValidateProgramARB");
		mgr->_glVertexAttrib1fv = (PFNGLVERTEXATTRIB1FVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib1fv" : "glVertexAttrib1fvARB");
		mgr->_glVertexAttrib2fv = (PFNGLVERTEXATTRIB2FVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib2fv" : "glVertexAttrib2fvARB");
		mgr->_glVertexAttrib3fv = (PFNGLVERTEXATTRIB3FVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib3fv" : "glVertexAttrib3fvARB");
		mgr->_glVertexAttrib4fv = (PFNGLVERTEXATTRIB4FVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib4fv" : "glVertexAttrib4fvARB");
		mgr->_glVertexAttrib1dv = (PFNGLVERTEXATTRIB1DVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib1dv" : "glVertexAttrib1dvARB");
		mgr->_glVertexAttrib2dv = (PFNGLVERTEXATTRIB2DVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib2dv" : "glVertexAttrib2dvARB");
		mgr->_glVertexAttrib3dv = (PFNGLVERTEXATTRIB3DVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib3dv" : "glVertexAttrib3dvARB");
		mgr->_glVertexAttrib4dv = (PFNGLVERTEXATTRIB4DVPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttrib4dv" : "glVertexAttrib4dvARB");
		mgr->_glVertexAttribPointer = (PFNGLVERTEXATTRIBPOINTERPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glVertexAttribPointer" : "glVertexAttribPointerARB");
		if(v20)
		{
			sscanf(glGetString(GL_SHADING_LANGUAGE_VERSION), "%d.%d", &uiMajor, &uiMinor);
			mgr->glsl_version = (uiMajor<<8) | uiMinor;
		}
		else
		{
			mgr->glsl_version = 0x100;
			mgr->_glGetHandle = (PFNGLGETHANDLEARBPROC)GI_GL_PROC_ADDRESS("glGetHandleARB");
		}

		/* check if vertex texturing supported */
		glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &iVTUnits);
		mgr->vertex_texture = (iVTUnits > 0);

		/* check if geometry shader supported */
		mgr->geometry_shader = (strstr(szExtensions, "GL_EXT_geometry_shader4") != NULL);
		if(mgr->geometry_shader)
		{
			/* get function pointers */
			mgr->_glProgramParameteri = (PFNGLPROGRAMPARAMETERIEXTPROC)
				GI_GL_PROC_ADDRESS("glProgramParameteriEXT");
		}

		/* create sampling shaders */
		//for(i=0; i<GI_NUM_SAMPLERS; ++i)
		//{
		//	mgr->gim_sampler[i] = (GIGLShader*)GI_CALLOC_SINGLE(sizeof(GIGLShader));
		//	GIGLShader_construct(mgr->gim_sampler[i], mgr, 
		//		g_sample_vs, NULL, g_sample_fs, 0, 0, 0, 1, szDefine[i]);
		//	if(!mgr->gim_sampler[i]->program)
		//	{
		//		GI_FREE_SINGLE(mgr->gim_sampler[i], sizeof(GIGLShader));
		//		mgr->gim_sampler[i] = NULL;
		//	}
		//}

		/* create rendering shaders */
		for(i=0; i<4; ++i)
		{
			if(((i & GI_GL_USE_VERTEX_TEXTURE) && !mgr->vertex_texture) || 
			   ((i & GI_GL_USE_GEOMETRY_SHADER) && !mgr->geometry_shader))
				mgr->gim_renderer[i] = NULL;
			else
			{
				mgr->gim_renderer[i] = (GIGLShader*)GI_CALLOC_SINGLE(sizeof(GIGLShader));
				if(i & GI_GL_USE_GEOMETRY_SHADER)
					GIGLShader_construct(mgr->gim_renderer[i], mgr, g_render_vs, g_render_gs, 
						g_render_fs, (i&GI_GL_USE_VERTEX_TEXTURE) ? GL_POINTS : GL_LINES_ADJACENCY_ARB, 
						GL_TRIANGLE_STRIP, 4, (i&GI_GL_USE_VERTEX_TEXTURE) ? 2 : 1, "GSHADER", "TEXTURE");
				else
					GIGLShader_construct(mgr->gim_renderer[i], mgr, g_render_vs, NULL, 
						g_render_fs, 0, 0, 0, (i&GI_GL_USE_VERTEX_TEXTURE) ? 1 : 0, "TEXTURE");
				if(!mgr->gim_renderer[i]->program)
				{
					GI_FREE_SINGLE(mgr->gim_renderer[i], sizeof(GIGLShader));
					mgr->gim_renderer[i] = NULL;
				}
			}
		}
	}

	/* check if FBOs supported */
	mgr->fbo = (strstr(szExtensions, "GL_EXT_framebuffer_object") != NULL);
	if(mgr->fbo)
	{
		/* get function pointers */
		mgr->_glBindFramebuffer = (PFNGLBINDFRAMEBUFFEREXTPROC)
			GI_GL_PROC_ADDRESS("glBindFramebufferEXT");
		mgr->_glBindRenderbuffer = (PFNGLBINDRENDERBUFFEREXTPROC)
			GI_GL_PROC_ADDRESS("glBindRenderbufferEXT");
		mgr->_glCheckFramebufferStatus = (PFNGLCHECKFRAMEBUFFERSTATUSEXTPROC)
			GI_GL_PROC_ADDRESS("glCheckFramebufferStatusEXT");
		mgr->_glDeleteFramebuffers = (PFNGLDELETEFRAMEBUFFERSEXTPROC)
			GI_GL_PROC_ADDRESS("glDeleteFramebuffersEXT");
		mgr->_glDeleteRenderbuffers = (PFNGLDELETERENDERBUFFERSEXTPROC)
			GI_GL_PROC_ADDRESS("glDeleteRenderbuffersEXT");
		mgr->_glFramebufferRenderbuffer = (PFNGLFRAMEBUFFERRENDERBUFFEREXTPROC)
			GI_GL_PROC_ADDRESS("glFramebufferRenderbufferEXT");
		mgr->_glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DEXTPROC)
			GI_GL_PROC_ADDRESS("glFramebufferTexture2DEXT");
		mgr->_glGenFramebuffers = (PFNGLGENFRAMEBUFFERSEXTPROC)
			GI_GL_PROC_ADDRESS("glGenFramebuffersEXT");
		mgr->_glGenRenderbuffers = (PFNGLGENRENDERBUFFERSEXTPROC)
			GI_GL_PROC_ADDRESS("glGenRenderbuffersEXT");
		mgr->_glRenderbufferStorage = (PFNGLRENDERBUFFERSTORAGEEXTPROC)
			GI_GL_PROC_ADDRESS("glRenderbufferStorageEXT");
		glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS_EXT, &mgr->max_attachments);
	}

	/* check if floating point color buffers, pixel and vertex data supported */
	mgr->color_buffer_float = (v30 || strstr(szExtensions, "GL_ARB_color_buffer_float") != NULL);
	if(mgr->color_buffer_float)
		mgr->_glClampColor = (PFNGLCLAMPCOLORPROC)GI_GL_PROC_ADDRESS(
			v30 ? "glClampColor" : "glClampColorARB");
	mgr->half_float_pixel = (v30 || strstr(szExtensions, "GL_ARB_half_float_pixel") != NULL);
	mgr->half_float_vertex = (v30 || strstr(szExtensions, "GL_ARB_half_float_vertex") != NULL);

	/* check if floating point and RG textures supported */
	mgr->texture_float = (v30 || strstr(szExtensions, "GL_ARB_texture_float") != NULL);
	mgr->texture_rg = (v30 || strstr(szExtensions, "GL_ARB_texture_rg") != NULL);

	/* check if MRTs supported */
	mgr->draw_buffers = (v20 || strstr(szExtensions, "GL_ARB_draw_buffers") != NULL);
	if(mgr->draw_buffers)
	{
		/* get function pointers */
		mgr->_glDrawBuffers = (PFNGLDRAWBUFFERSPROC)
			GI_GL_PROC_ADDRESS(v20 ? "glDrawBuffers" : "glDrawBuffersARB");
		glGetIntegerv(GL_MAX_DRAW_BUFFERS, &mgr->max_draw_buffers);
	}
	else
		mgr->max_draw_buffers = 1;
#endif
}

/** \internal
 *  \brief GL manager destructor.
 *  \param mgr OpenGL manager to initialize
 *  \ingroup opengl
 */
void GIGLManager_destruct(GIGLManager *mgr)
{
	GIuint i;
#if 0
	/* delete textures */
	if(mgr->norm_cubemap[0])
		glDeleteTextures(1, mgr->norm_cubemap);
	if(mgr->norm_cubemap[1])
		glDeleteTextures(1, mgr->norm_cubemap+1);
	if(mgr->norm_cubemap[2])
		glDeleteTextures(1, mgr->norm_cubemap+2);

	/* delete shaders */
	for(i=0; i<GI_NUM_SAMPLERS; ++i)
	{
		if(mgr->gim_sampler[i])
		{
			GIGLShader_destruct(mgr->gim_sampler[i], mgr);
			GI_FREE_SINGLE(mgr->gim_sampler[i], sizeof(GIGLShader));
		}
	}
	for(i=0; i<4; ++i)
	{
		if(mgr->gim_renderer[i])
		{
			GIGLShader_destruct(mgr->gim_renderer[i], mgr);
			GI_FREE_SINGLE(mgr->gim_renderer[i], sizeof(GIGLShader));
		}
	}
#endif
	/* delete GIM cache */
	if(mgr->gim_cache_size)
	{
		GIGLGIMCache *pCache;
		GI_LIST_FOREACH(mgr->texCoord_cache, pCache)
			GIGLGIMCache_destruct(pCache, mgr);
		GI_LIST_NEXT(mgr->texCoord_cache, pCache)
		GI_LIST_FOREACH(mgr->index_cache, pCache)
			GIGLGIMCache_destruct(pCache, mgr);
		GI_LIST_NEXT(mgr->index_cache, pCache)
		GI_LIST_CLEAR(mgr->texCoord_cache, sizeof(GIGLGIMCache));
		GI_LIST_CLEAR(mgr->index_cache, sizeof(GIGLGIMCache));
	}
#if 0
	if(mgr->gim_buffer)
		mgr->_glDeleteBuffers(1, &mgr->gim_buffer);
#endif
	memset(mgr, 0, sizeof(GIGLManager));
}

/** \internal
 *  \brief Get normalization cube map texture.
 *  \param mgr OpenGL manager to work on
 *  \param type type to create for
 *  \return texture object of normalization cube map or 0 if not supported
 *  \ingroup opengl
 */
GIuint GIGLManager_normalization_cubemap(GIGLManager *mgr, GIenum type)
{
#if 0
	static const GIfloat half[3] = { 0.5f, 0.5f, 0.5f };
	GIfloat pData[GI_NORMALIZATION_CUBE_SIZE*GI_NORMALIZATION_CUBE_SIZE*3];
	GIint iBoundTex, iSize = GI_NORMALIZATION_CUBE_SIZE, i, j, f;
	GIfloat fHS = GI_NORMALIZATION_CUBE_SIZE >> 1;
	GIfloat *vec;
	GIuint t = (type==GI_UNSIGNED_BYTE) ? 0 : ((type==GI_HALF_FLOAT) ? 1 : 2);
	if(mgr->norm_cubemap[t])
		return mgr->norm_cubemap[t];
	else if(!mgr->texture_cube_map || 
		(type != GI_UNSIGNED_BYTE && !mgr->texture_float))
		return 0;

	/* init texture */
	glGetIntegerv(GL_TEXTURE_BINDING_CUBE_MAP, &iBoundTex);
	glGenTextures(1, mgr->norm_cubemap+t);
	glBindTexture(GL_TEXTURE_CUBE_MAP, mgr->norm_cubemap[t]);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	/* set faces */
	for(f=0; f<6; ++f)
	{
		for(j=0,vec=pData; j<iSize; ++j)
		{
			for(i=0; i<iSize; ++i,vec+=3)
			{
				switch(f)
				{
				case 0:
					GI_VEC3_SET(vec, fHS, -(j+0.5f-fHS), -(i+0.5f-fHS));
					break;
				case 1:
					GI_VEC3_SET(vec, -fHS, -(j+0.5f-fHS), i+0.5f-fHS);
					break;
				case 2:
					GI_VEC3_SET(vec, i+0.5f-fHS, fHS, j+0.5f-fHS);
					break;
				case 3:
					GI_VEC3_SET(vec, i+0.5f-fHS, -fHS, -(j+0.5f-fHS));
					break;
				case 4:
					GI_VEC3_SET(vec, i+0.5f-fHS, -(j+0.5f-fHS), fHS);
					break;
				case 5:
					GI_VEC3_SET(vec, -(i+0.5f-fHS), -(j+0.5f-fHS), -fHS);
				}
				GIvec3f_normalize(vec);
				if(type == GI_UNSIGNED_BYTE)
				{
					GI_VEC3_ADD_SCALED(vec, half, vec, 0.5f);
				}
			}
		}
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X+f, 0, 
			(type==GI_UNSIGNED_BYTE) ? GL_RGB8 : 
			((type==GI_HALF_FLOAT) ? GL_RGBA16F : GL_RGBA32F), 
			iSize, iSize, 0, GL_RGB, GL_FLOAT, pData);
	}

	/* restore state */
	glBindTexture(GL_TEXTURE_CUBE_MAP, iBoundTex);
	return mgr->norm_cubemap[t];
#endif
	return 0;
}


void fileConv(char* wrkfile, char* filename, int* iNumIndices, int* iNumVertices, unsigned int **pIndices, GLfloat **pVertices, GLfloat **pUV, GLfloat **pNormals);
void GIAPIENTRY putMesh(char* filename, int* iNumIndices, int* iNumVertices, unsigned int **pIndices, GLfloat **pVertices, GLfloat **pUV, GLfloat **pNormals, double scale)
{
	int t;
	FILE* fp;
	GIRenderer *pRenderer = &(GIContext_current()->renderer);
	GIMesh *pMesh = pRenderer->context->mesh;
	GIPatch *pPatch, *pPEnd;
	GIFace *pFace, *pFEnd;
	GIHalfEdge *pHalfEdge;
	GIRenderAttrib arrAttribs[GI_ATTRIB_COUNT];
	GIuint i, j, a, uiNumAttribs;

	char tmpfile[256];

	if (iNumIndices) *iNumIndices = 0;
	if (iNumVertices) *iNumVertices = 0;
	if (pIndices) *pIndices = NULL;
	if (pVertices) *pVertices = NULL;
	if (pUV) *pUV = NULL;
	if (pNormals) *pNormals = NULL;

	/* setup everything */
	if(!GIRenderer_setup_mesh_rendering(pRenderer, 
		arrAttribs, &uiNumAttribs, &pPatch, &pPEnd))
		return;
	if(pPatch)
	{
		pFace = pPatch->faces;
		pFEnd = pPEnd->faces;
	}
	else
		pFace = pFEnd = pMesh->faces;

	t = time(NULL);
	sprintf(tmpfile, "tmp$$$_%04d.tmp", t);

	fp = fopen(tmpfile, "w");
		do
		{
			GIHalfEdge *wrk = pFace->hedges;
#if 0
			bool boundary = false;
			for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
			{
				if ( pHalfEdge->pstart->cut_hedge ) boundary = true;
			}
			if ( boundary )
			{
				pHalfEdge = wrk;
				for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
				{
					for(j=0; j<uiNumAttribs; ++j)
					{
						a = arrAttribs[j].attrib;
						switch(pMesh->asemantic[a])
						{
						case GI_POSITION_ATTRIB:
							//printf("v %d %.16f %.16f %.16f\n", pHalfEdge->vstart->id, (double)pHalfEdge->vstart->coords[0], (double)pHalfEdge->vstart->coords[1], (double)pHalfEdge->vstart->coords[2]);
							break;
						case GI_PARAM_ATTRIB:
							printf("cut %d\n", (int)(pHalfEdge->pstart->cut_hedge));
							printf("UV %d %.16f %.16f\n", pHalfEdge->pstart->vertex->id, (double)pHalfEdge->pstart->params[0], (double)pHalfEdge->pstart->params[1]);
							break;
						}
					}
				}
				printf("--------------\n\n");
				//pFace = pFace->next;
				//continue;
			}
#endif
			pHalfEdge = wrk;
			for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
			{
				for(j=0; j<uiNumAttribs; ++j)
				{
					int cut = 0;
					a = arrAttribs[j].attrib;
					switch(pMesh->asemantic[a])
					{
					case GI_POSITION_ATTRIB:
						if(arrAttribs[j].index < 0)
						{
							fprintf(fp, "v %d %.16f %.16f %.16f\n", pHalfEdge->vstart->id, (double)pHalfEdge->vstart->coords[0] / scale, (double)pHalfEdge->vstart->coords[1] / scale, (double)pHalfEdge->vstart->coords[2] / scale);
						}else
						{
							fprintf(fp, "vertex %d \n", pHalfEdge->vstart->id);
						}
						break;
					case GI_PARAM_ATTRIB:
						cut = (pHalfEdge->pstart->cut_hedge != NULL )? 1: 0;
						if(arrAttribs[j].index < 0)
						{
							fprintf(fp, "UV %d %.16f %.16f %d\n", pHalfEdge->pstart->vertex->id, (double)pHalfEdge->pstart->params[0], (double)pHalfEdge->pstart->params[1], cut);
						}
						else
						{
							fprintf(fp, "UV %d %.16f %.16f %d\n", pHalfEdge->pstart->vertex->id, (double)pHalfEdge->pstart->params[0], (double)pHalfEdge->pstart->params[1], cut);
						}
						break;
					default:
						if(arrAttribs[j].index < 0)
						{
							GIfloat* n = (GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]);
							fprintf(fp, "NORMAL %d %.16f %.16f %.16f\n", pHalfEdge->pstart->vertex->id, n[0], n[1], n[2]);
						}
						else
						{
							GIfloat* n = (GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]);
							fprintf(fp, "NORMAL %d %.16f %.16f %.16f\n", pHalfEdge->pstart->vertex->id, n[0], n[1], n[2]);
						}
					}
				}
			}
			pFace = pFace->next;
		}while(pFace != pFEnd);
	fclose(fp);

	fileConv(tmpfile, filename, iNumIndices, iNumVertices, pIndices, pVertices, pUV, pNormals);
	remove(tmpfile);
}

struct vertex_type
{
	int id;
	double v[3];
	double uv[2];
	double normal[3];
	int cut;
};

struct face_type
{
	unsigned int vertex_id[3];
};

void fileConv(char* wrkfile, char* filename, int* iNumIndices, int* iNumVertices, unsigned int **pIndices, GLfloat **pVertices, GLfloat **pUV, GLfloat **pNormals)
{
	char buf[256];
	FILE* fp;
	int index;
	int vrtnum, vrtnum0;
	int facenum;
	int cutnum;
	struct vertex_type* vertex;
	struct face_type* face;
	int i;
	int j;

	cutnum = 0;
	facenum = 0;
	vrtnum = 0;
	fp = fopen(wrkfile, "r");
	do
	{
		struct vertex_type v;
		if ( fgets(buf, 256, fp) == NULL )
		{
			break;
		}
		sscanf(buf, "UV %d %lf %lf %d", &index, v.uv, v.uv+1, &(v.cut));
		fgets(buf, 256, fp);
		cutnum += v.cut;
		fgets(buf, 256, fp);
		sscanf(buf, "v %d", &index);
		if ( vrtnum < index )
		{
			vrtnum = index;
		}
		facenum++;
	}while(1);
	fclose( fp );

	printf("vrtnum %d cutnum %d\n", vrtnum, cutnum);
	vrtnum += 1;
	vrtnum0 = vrtnum;
	vrtnum += cutnum;
	printf("vrtnum0 %d vrtnum %d cutnum %d\n", vrtnum0, vrtnum, cutnum);
	facenum /= 3;
	vertex = (struct vertex_type*)malloc( vrtnum*sizeof(struct vertex_type));
	face = (struct face_type*)malloc( facenum*sizeof(struct face_type));

	memset(vertex, '\0', vrtnum*sizeof(struct vertex_type));
	memset(face, '\0', facenum*sizeof(struct face_type));

	for ( i = 0; i < vrtnum; i++ )
	{
		vertex[i].id = -1;
	}

	fp = fopen(wrkfile, "r");

	int cutnum2 = 0;
	int cutindex = vrtnum0;
	for ( i = 0; i < facenum; i++ )
	{
		struct vertex_type v;
		struct face_type f;

		for ( j = 0; j < 3; j++ )
		{
			fgets(buf, 256, fp);
			sscanf(buf, "UV %d %lf %lf %d", &index, v.uv, v.uv+1, &(v.cut));
			fgets(buf, 256, fp);
			sscanf(buf, "NORMAL %d %lf %lf %lf", &index, v.normal, v.normal+1, v.normal+2);
			fgets(buf, 256, fp);
			sscanf(buf, "v %d %lf %lf %lf", &index, v.v, v.v+1, v.v+2);
			v.id = index;

			if ( vertex[index].id == -1 && vertex[index].cut == 0 )
			{
				vertex[index] = v;
			}else
			{
				double du = (vertex[index].uv[0] - v.uv[0]);
				double dv = (vertex[index].uv[1] - v.uv[1]);
				if ( (v.cut || vertex[index].cut) && (fabs(du) > 1.0e-8 || fabs(dv) > 1.0e-8) )
				{
					index = -1;
					for ( int k = 0; k < cutnum2; k++ )
					{
						if ( vertex[vrtnum0+k].id != v.id )
						{
							continue;
						}

						double du = (vertex[vrtnum0+k].uv[0] - v.uv[0]);
						double dv = (vertex[vrtnum0+k].uv[1] - v.uv[1]);
						if ( (fabs(du) > 1.0e-8 || fabs(dv) > 1.0e-8) )
						{
							continue;
						}
						index = vrtnum0+k;
						break;
					}
					if ( index == -1 )
					{
						vertex[cutindex] = v;
						index = cutindex;
						cutindex++;
						cutnum2++;
					}
				}else
				{
					if ( fabs(du) > 0.0001 || fabs(dv) > 0.0001 )
					{
						printf("boundary error %d u:%f u:%f   v:%f  v:%f =>du %.16f  dv %.16f\n", index, vertex[index].uv[0], v.uv[0], vertex[index].uv[1], v.uv[1], du, dv);
					}
				}
			}
			face[i].vertex_id[j] = index;
		}
	}
	fclose( fp );
	printf("cutnum2 %d\n", cutnum2);

	vrtnum = vrtnum0 + cutnum2;
	printf("vrtnum %d %d+%d\n", vrtnum, vrtnum0, cutnum2);
	if ( pIndices )
	{
		*iNumIndices = facenum;
		*iNumVertices = vrtnum;
		*pIndices = (unsigned int*)malloc(facenum*sizeof(unsigned int)*3);
		*pVertices = (GLfloat*)malloc(vrtnum*sizeof(GLfloat)*3);
		*pNormals = (GLfloat*)malloc(vrtnum*sizeof(GLfloat)*3);
		*pUV = (GLfloat*)malloc(vrtnum*sizeof(GLfloat)*2);
	}

	if (filename != NULL )
	{
		fp = fopen( filename, "w");
	}else
	{
		fp = NULL;
	}

	if ( fp ) fprintf(fp, "mtllib %s.mtl\n", filename);
	for ( i = 0; i < vrtnum; i++ )
	{
		if ( fp ) fprintf(fp, "v %f %f %f\n", vertex[i].v[0], vertex[i].v[1], vertex[i].v[2]);

		if ( pIndices && *pVertices )
		{
			(*pVertices)[3*i] = vertex[i].v[0];
			(*pVertices)[3*i+1] = vertex[i].v[1];
			(*pVertices)[3*i+2] = vertex[i].v[2];
		}
	}
	for ( i = 0; i < vrtnum; i++ )
	{
		if ( fp ) fprintf(fp, "vt %f %f\n", vertex[i].uv[0], vertex[i].uv[1]);
		if ( pIndices && *pUV )
		{
			(*pUV)[2*i] = vertex[i].uv[0];
			(*pUV)[2*i+1] = vertex[i].uv[1];
		}
	}
	for ( i = 0; i < vrtnum; i++ )
	{
		if ( fp ) fprintf(fp, "vn %f %f %f\n", vertex[i].normal[0], vertex[i].normal[1], vertex[i].normal[2]);
		if ( pIndices && *pNormals )
		{
			(*pNormals)[3*i] = vertex[i].normal[0];
			(*pNormals)[3*i+1] = vertex[i].normal[1];
			(*pNormals)[3*i+2] = vertex[i].normal[2];
		}
	}

	if ( fp ) fprintf(fp, "usemtl mat\n");
	for ( i = 0; i < facenum; i++ )
	{
		if ( fp ) fprintf( fp, "f %d/%d/%d",   face[i].vertex_id[0]+1, face[i].vertex_id[0]+1, face[i].vertex_id[0]+1);
		if ( fp ) fprintf( fp, "  %d/%d/%d",   face[i].vertex_id[1]+1, face[i].vertex_id[1]+1, face[i].vertex_id[1]+1);
		if ( fp ) fprintf( fp, "  %d/%d/%d\n", face[i].vertex_id[2]+1, face[i].vertex_id[2]+1, face[i].vertex_id[2]+1);
		//fprintf( fp, "f %d %d %d\n",   face[i].vertex_id[0]+1, face[i].vertex_id[1]+1, face[i].vertex_id[2]+1);
		if ( pIndices && *pIndices )
		{
			(*pIndices)[3*i] = face[i].vertex_id[0];
			(*pIndices)[3*i+1] = face[i].vertex_id[1];
			(*pIndices)[3*i+2] = face[i].vertex_id[2];
		}
	}
	if ( fp ) fclose(fp);
	
	if ( fp ) 
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
	free(vertex);
	free(face);
}


void  putCut(double r, N_Cylinder& boundaryLine, double scale)
{
	r /= scale;
	GIRenderer *pRenderer = &(GIContext_current()->renderer);
	GIMesh *pMesh = pRenderer->context->mesh;
	GIPatch *pPatch, *pPEnd;
	GIHalfEdge *pHalfEdge;
	GIParam *pParam;
	GIRenderAttrib arrAttribs[GI_ATTRIB_COUNT];
	GIuint j, a, uiNumAttribs;
	FILE* fp;
	char tmpfile[256];

	/* setup everything */
	if(!GIRenderer_setup_mesh_rendering(pRenderer, 
		arrAttribs, &uiNumAttribs, &pPatch, &pPEnd) || !pMesh->patches)
		return;

	int t = time(NULL);
	sprintf(tmpfile, "tmp$$$_%04d.tmp", t);
	
	fp = fopen( tmpfile, "w");
	std::vector<int> vcount;

	int vnum = 0;

	do
	{
		pParam = pPatch->params;
		do
		{
			pHalfEdge = pParam->cut_hedge;
			for(j=0; j<uiNumAttribs; ++j)
			{
				a = arrAttribs[j].attrib;
				switch(pMesh->asemantic[a])
				{
				case GI_POSITION_ATTRIB:
					fprintf(fp, "vertex %d %f %f %f\n", pHalfEdge->vstart->id, (double)pHalfEdge->vstart->coords[0] / scale, (double)pHalfEdge->vstart->coords[1] / scale, (double)pHalfEdge->vstart->coords[2] / scale);
					vnum++;
					break;
				case GI_PARAM_ATTRIB:
					fprintf(fp, "UV %d %.16f %.16f\n", pHalfEdge->pstart->vertex->id, (double)pHalfEdge->pstart->params[0], (double)pHalfEdge->pstart->params[1]);
				}
			}
			pParam = pHalfEdge->next->pstart;
		}while(pParam != pPatch->params);
		vcount.push_back(vnum);
		vnum = 0;
		pPatch = pPatch->next;
	}while(pPatch != pPEnd);

	fclose( fp );



	fp = fopen( tmpfile, "r");
	char buf[256];
	double p1[3];
	double p2[3];
	int id[2];
	double uv[2][2];

	printf("boundary line %d\n", vcount.size());
	for ( int j = 0; j < vcount.size(); j++ )
	{
		fgets(buf, 256, fp);
		sscanf(buf, "UV %d %lf %lf", id, uv[0], uv[0]+1);
		fgets(buf, 256, fp);
		sscanf(buf, "vertex %d %lf %lf %lf", &id, p1, p1+1, p1+2);
		for ( int i = 1; i < vcount[j]; i++ )
		{
			fgets(buf, 256, fp);
			sscanf(buf, "UV %d %lf %lf", id+1, uv[1], uv[1]+1);
			fgets(buf, 256, fp);
			sscanf(buf, "vertex %d %lf %lf %lf", &id, p2, p2+1, p2+2);
			Cylinder cyl(p1, p2, r);
			cyl.setVertexID(id);
			cyl.setUV(uv);

			boundaryLine.Add( cyl);

			p1[0] = p2[0];
			p1[1] = p2[1];
			p1[2] = p2[2];
			id[0] = id[1];
			uv[0][0] = uv[1][0];
			uv[0][1] = uv[1][1];
		}
		boundaryLine.close_loop();
	}
	fclose( fp );
	remove(tmpfile);
}

void calcBoundaryBox(double min[3], double max[3])
{
	min[0] = min[1] = min[2] = 0.0;
	max[0] = max[1] = max[2] = 0.0;
	GIContext *pContext = GIContext_current();
	if(pContext->mesh && pContext->mesh->id)
	{
		min[0] = pContext->mesh->aabb_min[0];
		min[1] = pContext->mesh->aabb_min[1];
		min[2] = pContext->mesh->aabb_min[2];
		max[0] = pContext->mesh->aabb_max[0];
		max[1] = pContext->mesh->aabb_max[1];
		max[2] = pContext->mesh->aabb_max[2];
		return;
	}
}
void calcBoundaryBox2(GIdouble min[3], GIdouble max[3])
{
	min[0] = min[1] = min[2] = 0.0;
	max[0] = max[1] = max[2] = 0.0;
	GIContext *pContext = GIContext_current();
	if(pContext->mesh && pContext->mesh->id)
	{
		min[0] = pContext->mesh->aabb_min[0];
		min[1] = pContext->mesh->aabb_min[1];
		min[2] = pContext->mesh->aabb_min[2];
		max[0] = pContext->mesh->aabb_max[0];
		max[1] = pContext->mesh->aabb_max[1];
		max[2] = pContext->mesh->aabb_max[2];
		return;
	}
}
