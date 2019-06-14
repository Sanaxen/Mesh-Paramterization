#if 0
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
 *  \brief Implementation of functions for geometry image sampling.
 */

#include "gi_sampler.h"
#include "gi_context.h"
#include "gi_memory.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#if OPENGI_SSE >= 1
	#if defined(_MSC_VER) && _MSC_VER >= 1400
		#include <intrin.h>
	#else
		#include <xmmintrin.h>
		#if OPENGI_SSE >= 2
			#include <emmintrin.h>
			#if OPENGI_SSE >= 3
				#include <pmmintrin.h>
			#endif
		#endif
	#endif
#endif

#define GI_GROUP_WITH_TEXTURE			0x80000000


/** internal
 *  \brief Identity matrix.
 *  \ingroup sampling
 */
static const GIfloat g_identity[16] = 
{ 
	1.0f, 0.0f, 0.0f, 0.0f, 
	0.0f, 1.0f, 0.0f, 0.0f, 
	0.0f, 0.0f, 1.0f, 0.0f, 
	0.0f, 0.0f, 0.0f, 1.0f
};

/** Set boolean configuration parameter for sampling.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup sampling
 */
void GIAPIENTRY giSamplerParameterb(GIenum pname, GIboolean param)
{
	GISampler *pSampler = &(GIContext_current()->sampler);

	/* select state and set value */
	switch(pname)
	{
	case GI_SAMPLER_USE_SHADER:
		pSampler->use_shader = param;
		break;
	case GI_SAMPLER_USE_RENDER_TO_TEXTURE:
		pSampler->use_fbo = param;
		break;
	default:
		GIContext_error(pSampler->context, GI_INVALID_ENUM);
	}
}

/** Set integer configuration parameter for sampling.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup sampling
 */
void GIAPIENTRY giSamplerParameteri(GIenum pname, GIint param)
{
	GISampler *pSampler = &(GIContext_current()->sampler);

	/* select state and set value */
	switch(pname)
	{
	case GI_SAMPLER:
		switch(param)
		{
		case GI_SAMPLER_SOFTWARE:
		case GI_SAMPLER_OPENGL:
			pSampler->sampler = param;
			break;
		default:
			GIContext_error(pSampler->context, GI_INVALID_ENUM);
		}
		break;
	default:
		GIContext_error(pSampler->context, GI_INVALID_ENUM);
	}
}

/** Set integer attribute sampling parameter.
 *  \param attrib attribute channel
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup sampling
 */
void GIAPIENTRY giAttribSamplerParameteri(GIuint attrib, 
										  GIenum pname, GIint param)
{
	GISampler *pSampler = &(GIContext_current()->sampler);
	if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pSampler->context, GI_INVALID_VALUE);
		return;
	}

	/* select state and set value */
	switch(pname)
	{
	case GI_SAMPLING_MODE:
		if(param == GI_SAMPLE_DEFAULT || param == GI_SAMPLE_NORMALIZED || 
			param == GI_SAMPLE_TEXTURED)
			pSampler->attrib_mode[attrib] = param;
		else
			GIContext_error(pSampler->context, GI_INVALID_VALUE);
		break;
	case GI_TEXTURE_DIMENSION:
		if(param >= 0 && param <= 4)
			pSampler->texture_dim[attrib] = param;
		else
			GIContext_error(pSampler->context, GI_INVALID_VALUE);
		break;
	case GI_GL_SAMPLE_TEXTURE:
		if(param >= 0)
			pSampler->attrib_texture[attrib] = param;
		else
			GIContext_error(pSampler->context, GI_INVALID_VALUE);
		break;
	default:
		GIContext_error(pSampler->context, GI_INVALID_ENUM);
	}
}

/** Set floating point array attribute sampling parameter.
 *  \param attrib attribute channel
 *  \param pname state to set
 *  \param params values to set
 *  \ingroup sampling
 */
void GIAPIENTRY giAttribSamplerParameterfv(GIuint attrib, GIenum pname, 
										   const GIfloat *params)
{
	GISampler *pSampler = &(GIContext_current()->sampler);
	if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pSampler->context, GI_INVALID_VALUE);
		return;
	}

	/* select state and set value */
	switch(pname)
	{
	case GI_SAMPLING_TRANSFORM:
		memcpy(pSampler->attrib_matrix[attrib], params, 16*sizeof(GIfloat));
		if(memcmp(params, g_identity, 16*sizeof(GIfloat)))
			pSampler->identity_matrix &= ~(1<<attrib);
		else
			pSampler->identity_matrix |= 1 << attrib;
		break;
	default:
		GIContext_error(pSampler->context, GI_INVALID_ENUM);
	}
}

/** Sample attributes of current mesh (with UVs) to according images.
 *  This function samples the specified attributes of the current bound mesh to 
 *  their according attribute images. The mesh has to have a valid parameterization.
 *  When using OpenGL acceleration an OpenGL context has to be current. 
 *  The function should preserve every state it changes of the current OpenGL context, so it may be 
 *  called from anywhere in the program (except between glBegin()/glEnd() of course).
 *  \ingroup sampling
 */
void GIAPIENTRY giSample()
{
	GISampler *pSampler = &(GIContext_current()->sampler);
	GIMesh *pMesh = pSampler->context->mesh;
	GIPatch *pPatch;

	/* error checking */
	pSampler->sampled_attribs = 0;
	if(!pMesh || (!pMesh->active_patch && pMesh->patch_count!=1))
	{
		GIContext_error(pSampler->context, GI_INVALID_OPERATION);
		return;
	}
	pPatch = pMesh->active_patch ? pMesh->active_patch : pMesh->patches;
	if(!pPatch->parameterized || !pPatch->corners[0])
	{
		GIContext_error(pSampler->context, GI_INVALID_PARAMETERIZATION);
		return;
	}

	/* sample using specified method */
	switch(pSampler->sampler)
	{
	case GI_SAMPLER_SOFTWARE:
		GISampler_sample_software(pSampler, pPatch);
		break;
	case GI_SAMPLER_OPENGL:
		GISampler_sample_opengl(pSampler, pPatch);
	}
}

/** \internal
 *  \brief Sampler constructor.
 *  \param sampler sampler to construct
 *  \param context context to construct in
 *  \ingroup sampling
 */
void GISampler_construct(GISampler *sampler, GIContext *context)
{
	GIuint a;

	/* initialize state */
	sampler->context = context;
	sampler->sampler = GI_SAMPLER_SOFTWARE;
	sampler->use_shader = GI_TRUE;
	sampler->use_fbo = GI_TRUE;
	sampler->sampled_attribs = 0;
	sampler->identity_matrix = ~0;
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
	{
		sampler->attrib_mode[a] = GI_SAMPLE_DEFAULT;
		sampler->attrib_texture[a] = 0;
		sampler->texture_dim[a] = 2;
		sampler->callback[a] = 0;
		sampler->cdata[a] = 0;
		memcpy(sampler->attrib_matrix[a], g_identity, 16*sizeof(GIfloat));
	}
}

/** \internal
 *  \brief Sample mesh using software rasterization.
 *  \param sampler sampler to use
 *  \param patch patch to sample
 *  \ingroup sampling
 */
void GISampler_sample_software(GISampler *sampler, GIPatch *patch)
{
	GIGLManager *pGL = sampler->context->gl_manager;
	GIMesh *pMesh = patch->mesh;
	GIImage *pImage, *pImage2;
	GIRasterizerData data;
	GIdouble *pCoords;
	GIImage arrRestore[GI_ATTRIB_COUNT];
	GIbitfield arrGroups[GI_ATTRIB_COUNT];
	GIuint arrGroupSize[GI_ATTRIB_COUNT];
	GIImage *arrImages[GI_ATTRIB_COUNT];
	GIpfunc pfnSetPixel[GI_ATTRIB_COUNT];
	GITexture arrTextures[GI_ATTRIB_COUNT];
	GIuint i, k, idx, uiAttribs = 0;
	GIint j, iNumGroups = 0;
	GIboolean bMatch;
	GLint iUnpackAlign, iBoundTex, iBoundBO;

	/* save GL state */
	if(!pGL->gl_version)
		GIGLManager_init(pGL);
	glGetIntegerv(GL_UNPACK_ALIGNMENT, &iUnpackAlign);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &iBoundTex);
	if(pGL->vbo)
		glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &iBoundBO);

	/* collect information about images to be sampled */
	memset(arrTextures, 0, GI_ATTRIB_COUNT*sizeof(GITexture));
	for(i=0; i<GI_ATTRIB_COUNT; ++i)
	{
		pfnSetPixel[i] = NULL;
		pImage = sampler->context->attrib_image[i];
		if(!pImage || (pMesh->asemantic[i] == GI_NONE && pMesh->aoffset[i]<0) || 
			(pMesh->asemantic[i] == GI_PARAM_STRETCH_ATTRIB && !patch->param_metric))
			continue;

		/* error checking */
		if(!(pImage->data || glIsTexture(pImage->texture) || pImage->buffer) || 
			(sampler->attrib_mode[i] == GI_SAMPLE_TEXTURED && 
			!glIsTexture(sampler->attrib_texture[i])))
		{
			GIContext_error(sampler->context, GI_INVALID_OPERATION);
			continue;
		}
		if((pImage->buffer && !pGL->vbo) || (sampler->attrib_mode[i] == 
			GI_SAMPLE_TEXTURED && ((sampler->texture_dim[i]==3 && 
			!pGL->texture_3d) || (sampler->texture_dim[i]==4 && 
			!pGL->texture_cube_map))))
		{
			GIContext_error(sampler->context, GI_UNSUPPORTED_OPERATION);
			continue;
		}
		if(pImage->buffer && !pGL->_glIsBuffer(pImage->buffer))
		{
			GIContext_error(sampler->context, GI_INVALID_OPERATION);
			continue;
		}

		/* create texture if neccessary */
		if(sampler->attrib_mode[i] == GI_SAMPLE_TEXTURED)
			GITexture_construct(arrTextures+i, 
				sampler->attrib_texture[i], sampler->texture_dim[i]);

		/* images with same size into same group */
		bMatch = GI_FALSE;
		for(j=iNumGroups-1; j>=0 && !bMatch; --j)
		{
			pImage2 = arrImages[j];
			if(pImage->sub_width == pImage2->sub_width && 
				pImage->sub_height == pImage2->sub_height)
			{
				arrGroups[j] |= 1 << i;
				++arrGroupSize[j];
				bMatch = GI_TRUE;
			}
		}
		if(!bMatch)
		{
			arrGroups[iNumGroups] = 1 << i;
			arrGroupSize[iNumGroups] = 1;
			arrImages[iNumGroups++] = pImage;
		}

		/* create data storage if needed */
		pfnSetPixel[i] = NULL;
		if(pImage->texture)
		{
			arrRestore[i] = *pImage;
			if(pImage->width != pImage->sub_width || 
				pImage->height != pImage->sub_height)
			{
				pImage->size = pImage->pixel_size * 
					pImage->sub_width * pImage->sub_height;
				pImage->width = pImage->sub_width;
				pImage->height = pImage->sub_height;
				pImage->offset_x = pImage->offset_y = 0;
			}
			if(pImage->type == GI_HALF_FLOAT && !pGL->half_float_pixel)
			{
				pImage->type = GI_FLOAT;
				pImage->size <<= 1;
				pImage->pixel_size <<= 1;
			}
			if(pImage->comp == 2 && !pGL->texture_rg)
			{
				pImage->comp <<= 1;
				pImage->size <<= 1;
				pImage->pixel_size <<= 1;
				pImage->gl_format = GL_RGBA;
			}
			else if(pImage->type == GI_FLOAT && pImage->comp == 4)
				pfnSetPixel[i] = &GIImage_setpixel_float4;
			pImage->data = pImage->sub_data = GI_MALLOC_ALIGNED(
				GI_SSE_SIZE(pImage->size), GI_SSE_ALIGN_FLOAT);
		}
		else if(pImage->buffer)
		{
			pGL->_glBindBuffer(GL_ARRAY_BUFFER, pImage->buffer);
			pImage->data = pGL->_glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			pImage->sub_data = (GIbyte*)pImage->data + pImage->pixel_size*
				(pImage->offset_y*pImage->width+pImage->offset_x);
		}

		/* set function pointers */
		if(pImage->type == GI_FLOAT)
			pfnSetPixel[i] = &GIImage_setpixel_float;
		else if(pImage->type == GI_HALF_FLOAT)
			pfnSetPixel[i] = &GIImage_setpixel_half;
		else
			pfnSetPixel[i] = &GIImage_setpixel_ubyte;
		uiAttribs |= 1 << i;
	}
	if(!iNumGroups)
		return;

	/* init rasterizer configuration */
	data.sampler = sampler;
	data.patch = patch;
	data.params = pCoords = (GIdouble*)GI_MALLOC_ARRAY(patch->pcount, 2*sizeof(GIdouble));
	data.textures = arrTextures;
	data.setpixel_fn = pfnSetPixel;

	/* process image groups */
	for(i=0; i<iNumGroups; ++i)
	{
		GIImage *pImage = arrImages[i];
		GIParam *pParam, **pCorners = patch->corners;
		GIdouble dWidth1 = (GIdouble)(pImage->sub_width-1);
		GIdouble dHeight1 = (GIdouble)(pImage->sub_height-1);
		GIdouble dEPS = DBL_EPSILON * GI_MAX(dWidth1, dHeight1);
		GIfloat fResRatio = dWidth1 / (GIfloat)(patch->resolution-1);
		GIboolean bRound = pImage->sub_width == pImage->sub_height && 
			fabs(fResRatio-floor(fResRatio)) < 1e-4;
		GIdouble dWidthEPS = dWidth1 /*+ dEPS*/, dHeightEPS = dHeight1 /*+ dEPS*/;

		/* scale params and readjust to integers */
		GI_LIST_FOREACH(patch->params, pParam)
			idx = pParam->id << 1;
			pCoords[idx] = dWidth1 * pParam->params[0];
			pCoords[idx+1] = dHeight1 * pParam->params[1];
		GI_LIST_NEXT(patch->params, pParam)
		pParam = pCorners[0];
		idx = pParam->id << 1;
		pCoords[idx+1] = 0.0;
		while(pParam != pCorners[1])
		{
			pParam = pParam->cut_hedge->next->pstart;
			idx = pParam->id << 1;
			pCoords[idx+1] = 0.0;
			if(bRound && (pParam->vertex->cut_degree != 2 || 
			   pParam->vertex->flags & GI_VERTEX_EXACT_BIT))
				pCoords[idx] = GI_ROUND(pCoords[idx]);
		}
		pCoords[idx] = dWidthEPS;
		while(pParam != pCorners[2])
		{
			pParam = pParam->cut_hedge->next->pstart;
			idx = pParam->id << 1;
			pCoords[idx] = dWidthEPS;
			if(bRound && (pParam->vertex->cut_degree != 2 || 
			   pParam->vertex->flags & GI_VERTEX_EXACT_BIT))
				pCoords[idx+1] = GI_ROUND(pCoords[idx+1]);
		}
		pCoords[idx+1] = dHeightEPS;
		while(pParam != pCorners[3])
		{
			pParam = pParam->cut_hedge->next->pstart;
			idx = pParam->id << 1;
			pCoords[idx+1] = dHeightEPS;
			if(bRound && (pParam->vertex->cut_degree != 2 || 
			   pParam->vertex->flags & GI_VERTEX_EXACT_BIT))
				pCoords[idx] = GI_ROUND(pCoords[idx]);
		}
		pCoords[idx] = 0.0;
		while(pParam != pCorners[0])
		{
			pParam = pParam->cut_hedge->next->pstart;
			idx = pParam->id << 1;
			pCoords[idx] = 0.0;
			if(bRound && (pParam->vertex->cut_degree != 2 || 
			   pParam->vertex->flags & GI_VERTEX_EXACT_BIT))
				pCoords[idx+1] = GI_ROUND(pCoords[idx+1]);
		}

		/* rasterizer configuration */
		data.num_attribs = arrGroupSize[i];
		for(j=0,k=0; j<GI_ATTRIB_COUNT; ++j)
		{
			if(arrGroups[i] & (1<<j))
			{
				data.attribs[k++] = j;
				if(sampler->attrib_mode[i] == GI_SAMPLE_TEXTURED)
				{
					GIfloat fRatioW = (GIfloat)pImage->sub_width / 
						(GIfloat)arrTextures[j].width[0];
					GIfloat fRatioH = (GIfloat)pImage->sub_height / 
						(GIfloat)arrTextures[j].height[0];
					if(fabs(fRatioW-1.0f) < fabs(fRatioH-1.0f))
						fRatioW = fRatioH;
					arrTextures[j].filter_mode = fRatioW>1.0f ? 
						arrTextures[j].mag_filter : arrTextures[j].min_filter;
					break;
				}
			}
		}

		/* rasterize faces */
#if OPENGI_NUM_THREADS > 1
		if(sampler->context->use_threads)
		{
			GIthread threads[OPENGI_NUM_THREADS-1];
			GIuint uiFaces;
			GIfloat *pTemp = GI_MALLOC_ALIGNED(GI_SSE_SIZE(
				data.num_attribs*12*sizeof(GIfloat)), GI_SSE_ALIGN_FLOAT);

			/* rasterize first face and create threads */
			GISampler_rasterize_triangle(patch->faces, &data, pTemp);
			data.next_face = patch->faces->next;
			data.end_face = patch->next->faces;
			GIMutex_construct(&data.mutex);
			for(j=0; j<OPENGI_NUM_THREADS-1; ++j)
				threads[j] = GIthread_create(GISampler_rasterize_thread, &data);

			/* print number of rasterized faces */
			uiFaces = (GIuint)GISampler_rasterize_thread(&data) + 1;
			GIDebug(printf("rasterized triangles: %d", uiFaces));
			for(j=0; j<OPENGI_NUM_THREADS-1; ++j)
			{
				uiFaces = (GIuint)GIthread_join(threads[j]);
				GIDebug(printf(", %d", uiFaces));
			}
			GIDebug(printf("\n"));
			GIMutex_destruct(&data.mutex);
			GI_FREE_ALIGNED(pTemp);
		}
		else
#endif
		{
			GIFace *pFace = patch->faces, *pFEnd = patch->next->faces;
			GIfloat *pTemp = GI_MALLOC_ALIGNED(GI_SSE_SIZE(
				data.num_attribs*12*sizeof(GIfloat)), GI_SSE_ALIGN_FLOAT);
			do
			{
				GISampler_rasterize_triangle(pFace, &data, pTemp);
				pFace = pFace->next;
			}while(pFace != pFEnd);
			GI_FREE_ALIGNED(pTemp);
		}

		/* fuse cut */
		if(bRound && patch->mesh->patch_count == 1 && patch->path_count != patch->groups)
		{
			GICutPath *pPath;
			GIuint *pPosMap = (GIuint*)GI_CALLOC_ARRAY(patch->groups, sizeof(GIuint));
			GIuint p1 = 0, p2, uiGLength = 0, uiResRatio = GI_ROUND(fResRatio);
			GI_LIST_FOREACH(patch->paths, pPath)
				uiGLength += uiResRatio * pPath->glength;
				if(pPath->twin && pPosMap[pPath->group])
				{
					p2 = pPosMap[pPath->group];
					for(; p1<uiGLength; ++p1,--p2)
					{
						for(j=0; j<data.num_attribs; ++j)
						{
							k = data.attribs[j];
							pImage2 = sampler->context->attrib_image[k];
							memcpy(GIImage_border_pixel(pImage2, p1), 
								GIImage_border_pixel(pImage2, p2), 
								pImage2->pixel_size);
						}
					}
				}
				else
				{
					p1 = uiGLength;
					pPosMap[pPath->group] = p1;
				}
			GI_LIST_NEXT(patch->paths, pPath)
			GI_FREE_ARRAY(pPosMap);
		}
	}
	GI_FREE_ARRAY(pCoords);

	/* copy data to texture or buffer if neccessary */
	for(i=0; i<GI_ATTRIB_COUNT; ++i)
	{
		if(uiAttribs & (1<<i))
		{
			/* copy and release data */
			pImage = sampler->context->attrib_image[i];
			if(pImage->texture)
			{
				glBindTexture(GL_TEXTURE_2D, pImage->texture);
				glTexSubImage2D(GL_TEXTURE_2D, 0, arrRestore[i].offset_x, 
					arrRestore[i].offset_y, pImage->sub_width, pImage->sub_height, 
					pImage->gl_format, pImage->type, pImage->data);
				GI_FREE_ALIGNED(pImage->data);
				*pImage = arrRestore[i];
			}
			else if(pImage->buffer)
			{
				pGL->_glBindBuffer(GL_ARRAY_BUFFER, pImage->buffer);
				pGL->_glUnmapBuffer(GL_ARRAY_BUFFER);
				pImage->data = pImage->sub_data = 0;
			}

			/* release texture data */
			if(sampler->attrib_mode[i] == GI_SAMPLE_TEXTURED)
				GITexture_destruct(arrTextures+i);
		}
	}
	sampler->sampled_attribs = uiAttribs;

	/* restore GL state */
	glBindTexture(GL_TEXTURE_2D, iBoundTex);
	if(pGL->vbo)
		pGL->_glBindBuffer(GL_ARRAY_BUFFER, iBoundBO);
	glPixelStorei(GL_UNPACK_ALIGNMENT, iUnpackAlign);
}

/** \internal
 *  \brief Sample mesh using hardware accelerated OpenGL rasterization.
 *  \param sampler sampler to use
 *  \param patch patch to sample
 *  \ingroup sampling
 */
void GISampler_sample_opengl(GISampler *sampler, GIPatch *patch)
{
	static const GIenum target[4] = { 
		GL_TEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_3D, GL_TEXTURE_CUBE_MAP };
	static const GIenum binding[4] = { GL_TEXTURE_BINDING_1D, 
		GL_TEXTURE_BINDING_2D, GL_TEXTURE_BINDING_3D, 
		GL_TEXTURE_BINDING_CUBE_MAP };
	GIGLManager *pGL = sampler->context->gl_manager;
	GIMesh *pMesh = patch->mesh;
	GIboolean bGLSL, bFBO, bTexFloat, bBufFloat;
	GLboolean bDepth, bTexture, bLight;
	GLint iPackAlign, iClamp, iTexUnit, 
		iBoundProg, iBoundFBO, iBoundRB, iBoundBO;
	GLuint i, uiFBO, uiRB = 0;

	/* collect information about supported and used GL features */
	if(!pGL->gl_version)
		GIGLManager_init(pGL);
	bGLSL = (pGL->glsl && sampler->use_shader && pGL->gim_sampler);
	bFBO = (pGL->fbo && sampler->use_fbo);
	bTexFloat = (pGL->texture_float && bFBO);
	bBufFloat = (pGL->color_buffer_float && bFBO);

	/* create FBO and RBs */
	if(bFBO)
	{
		glGetIntegerv(GL_FRAMEBUFFER_BINDING, &iBoundFBO);
		pGL->_glGenFramebuffers(1, &uiFBO);
		pGL->_glBindFramebuffer(GL_FRAMEBUFFER, uiFBO);
	}
	else
	{
		bDepth = glIsEnabled(GL_DEPTH_TEST);
		if(bDepth)
			glDisable(GL_DEPTH_TEST);
	}

	/* save and set GL state */
	glPushAttrib(GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT | GL_VIEWPORT_BIT);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_TEXTURE);
	glPushMatrix();
	glGetIntegerv(GL_PACK_ALIGNMENT, &iPackAlign);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDisable(GL_CULL_FACE);
	if(pGL->multitexture)
	{
		glGetIntegerv(GL_ACTIVE_TEXTURE, &iTexUnit);
		if(iTexUnit != GL_TEXTURE0)
			pGL->_glActiveTexture(GL_TEXTURE0);
	}
	if(pGL->pbo)
	{
		glGetIntegerv(GL_PIXEL_PACK_BUFFER_BINDING, &iBoundBO);
		pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
	}
	if(bGLSL)
	{
		if(pGL->gl_version >= 0x200)
			glGetIntegerv(GL_CURRENT_PROGRAM, &iBoundProg);
		else
			iBoundProg = pGL->_glGetHandle(GL_PROGRAM_OBJECT_ARB);
	}
	else
	{
		if(pGL->color_buffer_float)
		{
			glGetIntegerv(GL_CLAMP_VERTEX_COLOR, &iClamp);
			if(iClamp != GL_FIXED_ONLY)
				pGL->_glClampColor(GL_CLAMP_VERTEX_COLOR, GL_FIXED_ONLY);
		}
		bLight = glIsEnabled(GL_LIGHTING);
		glDisable(GL_LIGHTING);
		glMatrixMode(GL_TEXTURE);
		glPushMatrix();
		glLoadIdentity();
	}

	/* sample images */
	for(i=0; i<GI_ATTRIB_COUNT; ++i)
	{
		GIImage *pImage = sampler->context->attrib_image[i];
		GLuint uiFBStatus = GL_FRAMEBUFFER_COMPLETE;
		GLint iMinFilter, iMagFilter, iTexEnv = -1, iBoundTex, iDim = 0, iRowLength;
		if(!pImage || (pMesh->asemantic[i] == GI_NONE && pMesh->aoffset[i]<0) || 
			(pMesh->asemantic[i] == GI_PARAM_STRETCH_ATTRIB && !patch->param_metric))
			continue;

		/* error checking */
		if(!(pImage->data || glIsTexture(pImage->texture) || pImage->buffer) || 
			(sampler->attrib_mode[i] == GI_SAMPLE_TEXTURED && 
			!glIsTexture(sampler->attrib_texture[i])))
		{
			GIContext_error(sampler->context, GI_INVALID_OPERATION);
			continue;
		}
		if((pImage->buffer && !pGL->pbo) || 
			((pImage->data || pImage->buffer) && 
			((pImage->type==GI_HALF_FLOAT && !pGL->half_float_pixel) || 
			(pImage->type!=GI_UNSIGNED_BYTE && bFBO && !bBufFloat))) || 
			(pImage->texture && pImage->type!=GI_UNSIGNED_BYTE && !bTexFloat) || 
			(pImage->comp==2 && !pGL->texture_rg) || 
			(pImage->texture && pImage->comp==1 && !pGL->texture_rg && bFBO) || 
			(sampler->attrib_mode[i]==GI_SAMPLE_TEXTURED && 
			((sampler->texture_dim[i]==3 && !pGL->texture_3d) ||
			(sampler->texture_dim[i]==4 && !pGL->texture_cube_map))) || 
			(sampler->attrib_mode[i]==GI_SAMPLE_NORMALIZED && 
			!bGLSL && !GIGLManager_normalization_cubemap(pGL, pImage->type)))
		{
			GIContext_error(sampler->context, GI_UNSUPPORTED_OPERATION);
			continue;
		}
		if(pImage->buffer && !pGL->_glIsBuffer(pImage->buffer))
		{
			GIContext_error(sampler->context, GI_INVALID_OPERATION);
			continue;
		}

		/* FBOs supported? */
		if(bFBO)
		{
			if(pImage->texture)
			{
				/* set destination texture */
				if(pImage->type != GI_UNSIGNED_BYTE)
				{
					glGetIntegerv(GL_TEXTURE_BINDING_2D, &iBoundTex);
					glBindTexture(GL_TEXTURE_2D, pImage->texture);
					glGetTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, &iMinFilter);
					glGetTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, &iMagFilter);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
					glBindTexture(GL_TEXTURE_2D, iBoundTex);
				}
				pGL->_glFramebufferTexture2D(GL_FRAMEBUFFER, 
					GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pImage->texture, 0);
				glViewport(pImage->offset_x, pImage->offset_y, 
					pImage->sub_width, pImage->sub_height);
			}
			else
			{
				/* initialize and set destination RB */
				GIenum uiInternal = pImage->gl_internal;
				if(pImage->comp < 3 && !pGL->texture_rg)
				{
					if(pImage->type == GI_FLOAT)
						uiInternal = GL_RGBA32F;
					if(pImage->type == GI_HALF_FLOAT)
						uiInternal = GL_RGBA16F;
					else
						uiInternal = GL_RGBA8;
				}
				if(!uiRB)
				{
					glGetIntegerv(GL_RENDERBUFFER_BINDING, &iBoundRB);
					pGL->_glGenRenderbuffers(1, &uiRB);
					pGL->_glBindRenderbuffer(GL_RENDERBUFFER, uiRB);
				}
				pGL->_glRenderbufferStorage(GL_RENDERBUFFER, 
					uiInternal, pImage->sub_width, pImage->sub_height);
				pGL->_glFramebufferRenderbuffer(GL_FRAMEBUFFER, 
					GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, uiRB);
				glViewport(0, 0, pImage->sub_width, pImage->sub_height);
			}

			/* set draw buffer and check framebuffer status */
			glDrawBuffer(GL_COLOR_ATTACHMENT0);
			uiFBStatus = pGL->_glCheckFramebufferStatus(GL_FRAMEBUFFER);
		}
		else
			glViewport(0, 0, pImage->sub_width, pImage->sub_height);

		/* framebuffer complete (or not used)? */
		if(uiFBStatus != GL_FRAMEBUFFER_COMPLETE)
		{
			GIContext_error(sampler->context, GI_UNSUPPORTED_OPERATION);
			continue;
		}

		/* set render state */
//		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
//		glClear(GL_COLOR_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslated(0.5, 0.5, 0.0);
		glScaled(pImage->sub_width-1, pImage->sub_height-1, 1.0);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0.0, pImage->sub_width, 0.0, pImage->sub_height, -1.0, 1.0);
		glMatrixMode(GL_TEXTURE);
		glLoadMatrixf(sampler->attrib_matrix[i]);
		if(sampler->attrib_mode[i] == GI_SAMPLE_TEXTURED || 
			(sampler->attrib_mode[i] == GI_SAMPLE_NORMALIZED && !bGLSL))
		{
			iDim = (sampler->attrib_mode[i]==GI_SAMPLE_NORMALIZED) ? 
				4 : sampler->texture_dim[i];
			glGetIntegerv(binding[iDim-1], &iBoundTex);
			glBindTexture(target[iDim-1], 
				(sampler->attrib_mode[i]==GI_SAMPLE_NORMALIZED) ? 
				GIGLManager_normalization_cubemap(pGL, pImage->type) : 
				sampler->attrib_texture[i]);
		}

		/* sample mesh */
		if(bGLSL)
		{
			/* set shader */
			if(sampler->attrib_mode[i] == GI_SAMPLE_TEXTURED)
			{
				pGL->_glUseProgram(pGL->gim_sampler[sampler->texture_dim[i]]->program);
				pGL->_glUniform1i((GIint)GIHash_find(&pGL->gim_sampler[
					sampler->texture_dim[i]]->uniform_locs, "tex")-1, 0);
			}
			else if(sampler->attrib_mode[i] == GI_SAMPLE_NORMALIZED)
				pGL->_glUseProgram(pGL->gim_sampler[(pImage->type==GI_UNSIGNED_BYTE) ? 
					GI_SAMPLE_SHADER_PACKED : GI_SAMPLE_SHADER_NORMALIZED]->program);
			else
				pGL->_glUseProgram(pGL->gim_sampler[GI_SAMPLE_SHADER_DEFAULT]->program);
			GISampler_sample_gl_shader(sampler, patch, pGL, i);
		}
		else
		{
			/* set texture state */
			if(iDim)
			{
				bTexture = glIsEnabled(target[iDim-1]);
				glEnable(target[iDim-1]);
				glGetTexEnviv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, &iTexEnv);
				glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
			}
			GISampler_sample_gl_fixed_function(sampler, patch, pGL, i);
			if(iDim)
			{
				if(!bTexture)
					glDisable(target[iDim-1]);
				glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, iTexEnv);
			}
		}

		/* restore state */
		if(iDim)
			glBindTexture(target[iDim-1], iBoundTex);

		/* retrieve data */
		if(!pImage->texture)
		{
			glGetIntegerv(GL_PACK_ROW_LENGTH, &iRowLength);
			glPixelStorei(GL_PACK_ROW_LENGTH, pImage->width);
		}
		if(bFBO)
		{
			if(pImage->data)
			{
				glReadBuffer(GL_COLOR_ATTACHMENT0);
				glReadPixels(0, 0, pImage->sub_width, pImage->sub_height, 
					pImage->gl_format, pImage->type, pImage->sub_data);
			}
			else if(pImage->buffer)
			{
				glReadBuffer(GL_COLOR_ATTACHMENT0);
				pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, pImage->buffer);
				glReadPixels(0, 0, pImage->sub_width, pImage->sub_height, 
					pImage->gl_format, pImage->type, (GIbyte*)(
					(GIbyte*)pImage->sub_data-(GIbyte*)pImage->data));
				pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
			}
			else if(pImage->type != GI_UNSIGNED_BYTE)
			{
				glGetIntegerv(GL_TEXTURE_BINDING_2D, &iBoundTex);
				glBindTexture(GL_TEXTURE_2D, pImage->texture);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, iMinFilter);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, iMagFilter);
				glBindTexture(GL_TEXTURE_2D, iBoundTex);
			}
		}
		else if(pImage->texture)
		{
			glGetIntegerv(GL_TEXTURE_BINDING_2D, &iBoundTex);
			glBindTexture(GL_TEXTURE_2D, pImage->texture);
			glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 
				pImage->offset_x, pImage->offset_y, 
				0, 0, pImage->sub_width, pImage->sub_height);
			glBindTexture(GL_TEXTURE_2D, iBoundTex);
		}
		else if(pImage->buffer)
		{
			pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, pImage->buffer);
			glReadPixels(0, 0, pImage->sub_width, pImage->sub_height, 
				pImage->gl_format, pImage->type, (GIbyte*)(
				(GIbyte*)pImage->sub_data-(GIbyte*)pImage->data));
			pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
		}
		else
			glReadPixels(0, 0, pImage->sub_width, pImage->sub_height, 
				pImage->gl_format, pImage->type, pImage->sub_data);
		if(!pImage->texture)
			glPixelStorei(GL_PACK_ROW_LENGTH, iRowLength);
		sampler->sampled_attribs |= 1 << i;
	}

	/* delete FBO and RBs */
	if(bFBO)
	{
		pGL->_glBindFramebuffer(GL_FRAMEBUFFER, iBoundFBO);
		pGL->_glDeleteFramebuffers(1, &uiFBO);
		if(uiRB)
		{
			pGL->_glBindRenderbuffer(GL_RENDERBUFFER, iBoundRB);
			pGL->_glDeleteRenderbuffers(1, &uiRB);
		}
	}
	else if(bDepth)
		glEnable(GL_DEPTH_TEST);

	/* restore GL state */
	glMatrixMode(GL_TEXTURE);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glPixelStorei(GL_PACK_ALIGNMENT, iPackAlign);
	if(pGL->multitexture && iTexUnit != GL_TEXTURE0)
		pGL->_glActiveTexture(iTexUnit);
	if(pGL->pbo)
		pGL->_glBindBuffer(GL_PIXEL_PACK_BUFFER, iBoundBO);
	if(bGLSL)
		pGL->_glUseProgram(iBoundProg);
	else
	{
		if(pGL->color_buffer_float && iClamp != GL_FIXED_ONLY)
			pGL->_glClampColor(GL_CLAMP_VERTEX_COLOR, iClamp);
		if(bLight)
			glEnable(GL_LIGHTING);
		glMatrixMode(GL_TEXTURE);
		glPopMatrix();
	}
	glPopAttrib();
}

/** \internal
 *  \brief Sample mesh for use with fixed function pipeline.
 *  \param sampler sampler to use
 *  \param patch patch to sample
 *  \param gl valid OpenGL manager
 *  \param attrib attribute to sample
 *  \ingroup sampling
 */
void GISampler_sample_gl_fixed_function(GISampler *sampler, GIPatch *patch, 
										GIGLManager *gl, GIuint attrib)
{
	GIMesh *pMesh = patch->mesh;
	GIFace *pFace = patch->faces, *pFEnd = patch->next->faces;
	GIHalfEdge *pHalfEdge;
	GIParam *pParam;
	GIGLattribffunc glAttribXfv;
	GIGLattribdfunc glAttribXdv;
	GIfloat *fvec, *m = sampler->attrib_matrix[attrib];
	GIfloat v[4];
	GIuint i, j, uiSemantic = pMesh->asemantic[attrib];
	GIboolean b2, bTex = sampler->attrib_mode[attrib] == GI_SAMPLE_TEXTURED || 
		sampler->attrib_mode[attrib] == GI_SAMPLE_NORMALIZED, 
		bTrans = !bTex && !(sampler->identity_matrix&(1<<attrib));

	/* set attribute functions */
	switch(uiSemantic)
	{
	case GI_POSITION_ATTRIB:
		glAttribXdv = bTex ? &glTexCoord3dv : &glColor3dv;
		break;
	case GI_PARAM_ATTRIB:
		glAttribXdv = bTex ? &glTexCoord2dv : NULL;
		break;
	case GI_PARAM_STRETCH_ATTRIB:
		break;
	default:
		if(bTrans)
			glAttribXfv = NULL;
		else if(bTex)
		{
			switch(pMesh->asize[attrib])
			{
			case 1: glAttribXfv = &glTexCoord1fv; break;
			case 2: glAttribXfv = &glTexCoord2fv; break;
			case 3: glAttribXfv = &glTexCoord3fv; break;
			default: glAttribXfv = &glTexCoord4fv;
			}
		}
		else if(pMesh->asize[attrib] == 4)
			glAttribXfv = &glColor4fv;
		else if(pMesh->asize[attrib] == 3)
			glAttribXfv = &glColor3fv;
		else
		{
			glAttribXfv = NULL;
			b2 = pMesh->asize[attrib] == 2;
		}
	}

	/* render mesh */
	glBegin(GL_TRIANGLES);
		do
		{
			for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
			{
				switch(uiSemantic)
				{
				case GI_POSITION_ATTRIB:
					if(bTrans)
					{
						GIdouble *v = pHalfEdge->vstart->coords;
						glColor4f(m[0]*v[0]+m[4]*v[1]+m[8]*v[2]+m[12], 
								  m[1]*v[0]+m[5]*v[1]+m[9]*v[2]+m[13], 
								  m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14], 
								  m[3]*v[0]+m[7]*v[1]+m[11]*v[2]+m[15]);
					}
					else
						glAttribXdv(pHalfEdge->vstart->coords);
					break;
				case GI_PARAM_ATTRIB:
					if(glAttribXdv)
						glAttribXdv(pHalfEdge->pstart->params);
					else if(bTrans)
					{
						GIdouble *v = pHalfEdge->pstart->params;
						glColor4f(m[0]*v[0]+m[4]*v[1]+m[12], 
								  m[1]*v[0]+m[5]*v[1]+m[13], 
								  m[2]*v[0]+m[6]*v[1]+m[14], 
								  m[3]*v[0]+m[7]*v[1]+m[15]);
					}
					else
						glColor3d(pHalfEdge->pstart->params[0], 
							pHalfEdge->pstart->params[1], 0.0);
					break;
				case GI_PARAM_STRETCH_ATTRIB:
					if(bTex)
						glTexCoord1d(pHalfEdge->pstart->stretch);
					else if(bTrans)
						glColor4f(m[0]*pHalfEdge->pstart->stretch+m[12], 
								  m[1]*pHalfEdge->pstart->stretch+m[13], 
								  m[2]*pHalfEdge->pstart->stretch+m[14], 
								  m[3]*pHalfEdge->pstart->stretch+m[15]);
					else
						glColor3d(pHalfEdge->pstart->stretch, 0.0, 0.0);
					break;
				default:
					if(glAttribXfv)
						glAttribXfv((GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[attrib]));
					else if(bTrans)
					{
						v[1] = v[2] = 0.0f;
						v[3] = 1.0f;
						for(j=0; j<pMesh->asize[attrib]; ++j)
							v[j] = ((GIfloat*)((GIbyte*)pHalfEdge->astart+
								pMesh->aoffset[attrib]))[j];
						glColor4f(m[0]*v[0]+m[4]*v[1]+m[8]*v[2]+m[12]*v[3], 
								  m[1]*v[0]+m[5]*v[1]+m[9]*v[2]+m[13]*v[3], 
								  m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14]*v[3], 
								  m[3]*v[0]+m[7]*v[1]+m[11]*v[2]+m[15]*v[3]);
					}
					else
					{
						fvec = (GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[attrib]);
						glColor3f(fvec[0], b2 ? fvec[1] : 0.0f, 0.0f);
					}
				}
				glVertex2dv(pHalfEdge->pstart->params);
			}
			pFace = pFace->next;
		}while(pFace != pFEnd);
	glEnd();
	glBegin(GL_LINE_LOOP);
		pParam = patch->params;
		do
		{
			pHalfEdge = pParam->cut_hedge;
			switch(uiSemantic)
			{
			case GI_POSITION_ATTRIB:
				if(bTrans)
				{
					GIdouble *v = pHalfEdge->vstart->coords;
					glColor4f(m[0]*v[0]+m[4]*v[1]+m[8]*v[2]+m[12], 
							  m[1]*v[0]+m[5]*v[1]+m[9]*v[2]+m[13], 
							  m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14], 
							  m[3]*v[0]+m[7]*v[1]+m[11]*v[2]+m[15]);
				}
				else
					glAttribXdv(pHalfEdge->vstart->coords);
				break;
			case GI_PARAM_ATTRIB:
				if(glAttribXdv)
					glAttribXdv(pHalfEdge->pstart->params);
				else if(bTrans)
				{
					GIdouble *v = pHalfEdge->pstart->params;
					glColor4f(m[0]*v[0]+m[4]*v[1]+m[12], 
							  m[1]*v[0]+m[5]*v[1]+m[13], 
							  m[2]*v[0]+m[6]*v[1]+m[14], 
							  m[3]*v[0]+m[7]*v[1]+m[15]);
				}
				else
					glColor3d(pHalfEdge->pstart->params[0], 
						pHalfEdge->pstart->params[1], 0.0);
				break;
			case GI_PARAM_STRETCH_ATTRIB:
				if(bTex)
					glTexCoord1d(pHalfEdge->pstart->stretch);
				else if(bTrans)
					glColor4f(m[0]*pHalfEdge->pstart->stretch+m[12], 
							  m[1]*pHalfEdge->pstart->stretch+m[13], 
							  m[2]*pHalfEdge->pstart->stretch+m[14], 
							  m[3]*pHalfEdge->pstart->stretch+m[15]);
				else
					glColor3d(pHalfEdge->pstart->stretch, 0.0, 0.0);
				break;
			default:
				if(glAttribXfv)
					glAttribXfv((GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[attrib]));
				else if(bTrans)
				{
					v[1] = v[2] = 0.0f;
					v[3] = 1.0f;
					for(j=0; j<pMesh->asize[attrib]; ++j)
						v[j] = ((GIfloat*)((GIbyte*)pHalfEdge->astart+
							pMesh->aoffset[attrib]))[j];
					glColor4f(m[0]*v[0]+m[4]*v[1]+m[8]*v[2]+m[12]*v[3], 
							  m[1]*v[0]+m[5]*v[1]+m[9]*v[2]+m[13]*v[3], 
							  m[2]*v[0]+m[6]*v[1]+m[10]*v[2]+m[14]*v[3], 
							  m[3]*v[0]+m[7]*v[1]+m[11]*v[2]+m[15]*v[3]);
				}
				else
				{
					fvec = (GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[attrib]);
					glColor3f(fvec[0], b2 ? fvec[1] : 0.0f, 0.0f);
				}
			}
			glVertex2dv(pHalfEdge->pstart->params);
			pParam = pHalfEdge->next->pstart;
		}while(pParam != patch->params);
	glEnd();
}

/** \internal
 *  \brief Sample mesh for use with sampling shaders.
 *  \param sampler sampler to use
 *  \param patch patch to sample
 *  \param gl valid OpenGL manager
 *  \param attrib attribute to sample
 *  \ingroup sampling
 */
void GISampler_sample_gl_shader(GISampler *sampler, GIPatch *patch, 
								GIGLManager *gl, GIuint attrib)
{
	GIMesh *pMesh = patch->mesh;
	GIFace *pFace = patch->faces, *pFEnd = patch->next->faces;
	GIHalfEdge *pHalfEdge;
	GIParam *pParam;
	GIGLattribffunc glAttribXfv;
	GIuint i, uiSemantic = pMesh->asemantic[attrib];

	/* check attribute sizes */
	if(uiSemantic == GI_NONE)
	{
		switch(pMesh->asize[attrib])
		{
		case 1: glAttribXfv = &glTexCoord1fv; break;
		case 2: glAttribXfv = &glTexCoord2fv; break;
		case 3: glAttribXfv = &glTexCoord3fv; break;
		default: glAttribXfv = &glTexCoord4fv;
		}
	}

	/* render mesh */
	glBegin(GL_TRIANGLES);
		do
		{
			for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
			{
				switch(uiSemantic)
				{
				case GI_POSITION_ATTRIB:
					glTexCoord3dv(pHalfEdge->vstart->coords);
					break;
				case GI_PARAM_ATTRIB:
					glTexCoord2dv(pHalfEdge->pstart->params);
					break;
				case GI_PARAM_STRETCH_ATTRIB:
					glTexCoord1d(pHalfEdge->pstart->stretch);
					break;
				default:
					glAttribXfv((GIfloat*)((GIbyte*)pHalfEdge->astart+
						pMesh->aoffset[attrib]));
				}
				glVertex2dv(pHalfEdge->pstart->params);
			}
			pFace = pFace->next;
		}while(pFace != pFEnd);
	glEnd();
	glBegin(GL_LINE_LOOP);
		pParam = patch->params;
		do
		{
			pHalfEdge = pParam->cut_hedge;
			switch(uiSemantic)
			{
			case GI_POSITION_ATTRIB:
				glTexCoord3dv(pHalfEdge->vstart->coords);
				break;
			case GI_PARAM_ATTRIB:
				glTexCoord2dv(pHalfEdge->pstart->params);
				break;
			case GI_PARAM_STRETCH_ATTRIB:
				glTexCoord1d(pHalfEdge->pstart->stretch);
				break;
			default:
				glAttribXfv((GIfloat*)((GIbyte*)pHalfEdge->astart+
					pMesh->aoffset[attrib]));
			}
			glVertex2dv(pParam->params);
			pParam = pHalfEdge->next->pstart;
		}while(pParam != patch->params);
	glEnd();
}

/** \internal
 *  \brief Rasterize triangle in software sampling.
 *  \param face face to rasterize
 *  \param data rasterizer data
 *  \param temp temporary storage for vertex attributes
 *  \ingroup sampling
 */
void GISampler_rasterize_triangle(GIFace *face, GIRasterizerData *data, 
								  GIfloat *temp)
{
	static GIfloat vcHalf[4] = { 0.5f, 0.5f, 0.5f, 0.5f };
	GISampler *pSampler = data->sampler;
	GIContext *pContext = pSampler->context;
	GIPatch *pPatch = data->patch;
	GIMesh *pMesh = pPatch->mesh;
	GIHalfEdge *pHalfEdge;
	GIdouble *param[3];
	GIdouble *p1, *p2, *pTemp;
	GIdouble dMinY, dY, dX1, dX2;
	GIuint i, j, x, y, a, e1, e2, uiTemp, yMax, xMax;
	GIdouble incX[3], minX[3];
	GIuint lines[3], minY[3];
	GIdouble xBorder = pContext->attrib_image[data->attribs[0]]->sub_width - 1;
	GIdouble yBorder = pContext->attrib_image[data->attribs[0]]->sub_height - 1;

	/* gather face data */
	dMinY = DBL_MAX;
	for(i=0,pHalfEdge=face->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
	{
		/* compute extents and edge data */
		p1 = param[i] = data->params + (pHalfEdge->pstart->id<<1);
		p2 = data->params + (pHalfEdge->next->pstart->id<<1);
		dMinY = GI_MIN(dMinY, p1[1]);
		if(p1[1] > p2[1])
		{
			GI_SWAP(p1, p2, pTemp);
		}
		dY = ceil(p1[1]);
		lines[i] = ((p2[1]==yBorder) ? ((GIuint)yBorder+1) : (GIuint)ceil(p2[1])) - 
			((p1[1]==yBorder) ? ((GIuint)yBorder+1) : (GIuint)ceil(p1[1]));
		if(!lines[i])
		{
			incX[i] = minX[i] = 0.0;
			minY[i] = UINT_MAX;
		}
		else
		{
			incX[i] = (p2[0]-p1[0]) / (p2[1]-p1[1]);
			minX[i] = p1[0] + incX[i]*(dY-p1[1]);
			minY[i] = (GIuint)dY;
		}

		/* gather (and transform) attributes */
		for(j=0; j<data->num_attribs; ++j)
		{
			GIfloat *vec = temp + (j*12) + (i<<2);
#if OPENGI_SSE >= 1
			__m128 XMM7;
#else
			GIuint c, offset;
#endif
			a = data->attribs[j];
			switch(pMesh->asemantic[a])
			{
#if OPENGI_SSE >= 2
			case GI_POSITION_ATTRIB:
				{
					__m128d XMM0 = _mm_loadu_pd(pHalfEdge->vstart->coords);
					__m128d XMM1 = _mm_load_sd(pHalfEdge->vstart->coords+2);
					__m128 XMM2 = _mm_set1_ps(1.0f);
					XMM7 = _mm_cvtpd_ps(XMM0);
					XMM2 = _mm_cvtsd_ss(XMM2, XMM1);
					XMM7 = _mm_shuffle_ps(XMM7, XMM2, _MM_SHUFFLE(1, 0, 1, 0));
				}
				break;
			case GI_PARAM_ATTRIB:
				{
					__m128d XMM0 = _mm_loadu_pd(pHalfEdge->pstart->params);
					__m128 XMM1 = _mm_set_ss(1.0f);
					XMM7 = _mm_cvtpd_ps(XMM0);
					XMM7 = _mm_shuffle_ps(XMM7, XMM1, _MM_SHUFFLE(0, 1, 1, 0));
				}
				break;
			case GI_PARAM_STRETCH_ATTRIB:
				{
					__m128d XMM0 = _mm_load_sd(&pHalfEdge->pstart->stretch);
					__m128 XMM1 = _mm_set_ss(1.0f);
					XMM7 = _mm_cvtpd_ps(XMM0);
					XMM7 = _mm_shuffle_ps(XMM7, XMM1, _MM_SHUFFLE(0, 1, 1, 0));
				}
				break;
#else
			case GI_POSITION_ATTRIB:
				GI_VEC3_COPY(vec, pHalfEdge->vstart->coords);
				vec[3] = 1.0f;
#if OPENGI_SSE >= 1
				if(!(pSampler->identity_matrix & (1<<a)))
					XMM7 = _mm_load_ps(vec);
#endif
				break;
			case GI_PARAM_ATTRIB:
				GI_VEC2_COPY(vec, pHalfEdge->pstart->params);
				vec[2] = 0.0f;
				vec[3] = 1.0f;
#if OPENGI_SSE >= 1
				if(!(pSampler->identity_matrix & (1<<a)))
					XMM7 = _mm_load_ps(vec);
#endif
				break;
			case GI_PARAM_STRETCH_ATTRIB:
				vec[0] = pHalfEdge->pstart->stretch;
				vec[1] = vec[2] = 0.0f;
				vec[3] = 1.0f;
#if OPENGI_SSE >= 1
				if(!(pSampler->identity_matrix & (1<<a)))
					XMM7 = _mm_load_ps(vec);
#endif
				break;
#endif
#if OPENGI_SSE >= 1
			default:
				if(pMesh->asize[a] == 4)
					XMM7 = _mm_loadu_ps((GIfloat*)((GIbyte*)pHalfEdge
						->astart+pMesh->aoffset[a]));
				else
				{
					XMM7 = _mm_load_ss((GIfloat*)(
						(GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]));
					if(pMesh->asize[a] > 1)
					{
						__m128 XMM0 = _mm_load_ss((GIfloat*)(
							(GIbyte*)pHalfEdge->astart+pMesh->aoffset[a])+1);
						XMM7 = _mm_unpacklo_ps(XMM7, XMM0);
						if(pMesh->asize[a] > 2)
						{
							__m128 XMM1 = _mm_load_ss((GIfloat*)(
								(GIbyte*)pHalfEdge->astart+pMesh->aoffset[a])+2);
							XMM0 = _mm_set1_ps(1.0f);
							XMM0 = _mm_move_ss(XMM0, XMM1);
							XMM7 = _mm_shuffle_ps(XMM7, XMM0, _MM_SHUFFLE(1, 0, 1, 0));
						}
					}
					if(pMesh->asize[a] < 3)
					{
						__m128 XMM0 = _mm_set_ss(1.0f);
						XMM7 = _mm_shuffle_ps(XMM7, XMM0, _MM_SHUFFLE(0, 1, 1, 0));
					}
				}
			}
			if(!(pSampler->identity_matrix & (1<<a)))
			{
				__m128 XMM0 = _mm_loadu_ps(pSampler->attrib_matrix[a]);
				__m128 XMM1 = _mm_loadu_ps(pSampler->attrib_matrix[a]+4);
				__m128 XMM2 = _mm_loadu_ps(pSampler->attrib_matrix[a]+8);
				__m128 XMM3 = _mm_loadu_ps(pSampler->attrib_matrix[a]+12);
				__m128 XMM4 = _mm_shuffle_ps(XMM7, XMM7, _MM_SHUFFLE(0, 0, 0, 0));
				__m128 XMM5 = _mm_shuffle_ps(XMM7, XMM7, _MM_SHUFFLE(1, 1, 1, 1));
				__m128 XMM6 = _mm_shuffle_ps(XMM7, XMM7, _MM_SHUFFLE(2, 2, 2, 2));
				XMM7 = _mm_shuffle_ps(XMM7, XMM7, _MM_SHUFFLE(3, 3, 3, 3));
				XMM0 = _mm_mul_ps(XMM0, XMM4);
				XMM1 = _mm_mul_ps(XMM1, XMM5);
				XMM2 = _mm_mul_ps(XMM2, XMM6);
				XMM3 = _mm_mul_ps(XMM3, XMM7);
				XMM0 = _mm_add_ps(XMM0, XMM1);
				XMM0 = _mm_add_ps(XMM0, XMM2);
				XMM0 = _mm_add_ps(XMM0, XMM3);
				_mm_store_ps(vec, XMM0);
			}
#if OPENGI_SSE < 2
			else if(pMesh->asemantic[a] == GI_NONE)
#else
			else
#endif
				_mm_store_ps(vec, XMM7);
#else
			default:
				vec[1] = vec[2] = 0.0f; vec[3] = 1.0f;
				for(c=0,offset=pMesh->aoffset[a]; c<pMesh->asize[a]; 
					++c,offset+=sizeof(GIfloat))
					vec[c] = *(GIfloat*)((GIbyte*)pHalfEdge->astart+offset);
			}
			if(!(pSampler->identity_matrix & (1<<a)))
			{
				GIfloat *m = pSampler->attrib_matrix[a];
				GIfloat v[4] = { vec[0], vec[1], vec[2], vec[3] };
				vec[0] = m[0]*v[0] + m[4]*v[1] + m[8]*v[2] + m[12]*v[3];
				vec[1] = m[1]*v[0] + m[5]*v[1] + m[9]*v[2] + m[13]*v[3];
				vec[2] = m[2]*v[0] + m[6]*v[1] + m[10]*v[2] + m[14]*v[3];
				vec[3] = m[3]*v[0] + m[7]*v[1] + m[11]*v[2] + m[15]*v[3];
			}
#endif
		}
	}

	/* initialize active edge data */
	if(minY[0] < minY[1])
	{
		e1 = 0;
		e2 = (minY[1]<minY[2] ? 1 : 2);
	}
	else
	{
		e1 = 1;
		e2 = (minY[0]<minY[2] ? 0 : 2);
	}
	if(minX[e1] > minX[e2] || (minX[e1] == minX[e2] && incX[e1] > incX[e2]))
	{
		GI_SWAP(e1, e2, uiTemp);
	}
	dX1 = minX[e1];
	dX2 = minX[e2];

	/* process scanlines */
	for(y=ceil(dMinY),yMax=y+GI_MAX(lines[e1],lines[e2]); y<yMax; )
	{
#if OPENGI_SSE >= 2
		__m128d XMM7 = _mm_setzero_pd();
#if OPENGI_SSE < 3
		__m128d XMM6 = _mm_set_sd(-0.0);
		XMM6 = _mm_shuffle_pd(XMM6, XMM6, _MM_SHUFFLE2(0, 1));
#endif
		XMM7 = _mm_cvtsi32_sd(XMM7, y);
		XMM7 = _mm_shuffle_pd(XMM7, XMM7, _MM_SHUFFLE2(0, 1));
#else
		GIdouble p[2];
		p[1] = (GIdouble)y;
#endif
		for(x=ceil(dX1),xMax=(dX2==xBorder /*&& dX2!=dX1*/) ? 
			(dX2+1.0) : ceil(dX2); x<xMax; ++x)
		{
			/* compute barycentric coordinates */
#if OPENGI_SSE >= 1
			__m128 XMM11, XMM12, XMM13, XOUT;
#endif
#if OPENGI_SSE >= 2
			__m128d XMM0, XMM1, XMM2, XMM3, XMM4;
			XMM7 = _mm_cvtsi32_sd(XMM7, x);
			XMM1 = _mm_loadu_pd(param[1]);
			XMM2 = _mm_loadu_pd(param[2]);
			XMM3 = _mm_sub_pd(XMM1, XMM7);
			XMM4 = _mm_sub_pd(XMM2, XMM7);
			XMM4 = _mm_shuffle_pd(XMM4, XMM4, _MM_SHUFFLE2(0, 1));
			XMM3 = _mm_mul_pd(XMM3, XMM4);
#if OPENGI_SSE >= 3
			XMM3 = _mm_hsub_pd(XMM3, XMM3);
#else
			XMM4 = _mm_shuffle_pd(XMM3, XMM3, _MM_SHUFFLE2(0, 1));
			XMM3 = _mm_sub_pd(XMM3, XMM4);
			XMM3 = _mm_xor_pd(XMM3, XMM6);
#endif
			XMM0 = _mm_loadu_pd(param[0]);
			XMM2 = _mm_sub_pd(XMM2, XMM7);
			XMM4 = _mm_sub_pd(XMM0, XMM7);
			XMM4 = _mm_shuffle_pd(XMM4, XMM4, _MM_SHUFFLE2(0, 1));
			XMM2 = _mm_mul_pd(XMM2, XMM4);
#if OPENGI_SSE >= 3
			XMM2 = _mm_hsub_pd(XMM2, XMM2);
#else
			XMM4 = _mm_shuffle_pd(XMM2, XMM2, _MM_SHUFFLE2(0, 1));
			XMM2 = _mm_sub_pd(XMM2, XMM4);
			XMM2 = _mm_xor_pd(XMM2, XMM6);
#endif
			XMM0 = _mm_sub_pd(XMM0, XMM7);
			XMM1 = _mm_sub_pd(XMM1, XMM7);
			XMM1 = _mm_shuffle_pd(XMM1, XMM1, _MM_SHUFFLE2(0, 1));
			XMM1 = _mm_mul_pd(XMM1, XMM0);
#if OPENGI_SSE >= 3
			XMM1 = _mm_hsub_pd(XMM1, XMM1);
#else
			XMM4 = _mm_shuffle_pd(XMM1, XMM1, _MM_SHUFFLE2(0, 1));
			XMM1 = _mm_sub_pd(XMM1, XMM4);
			XMM1 = _mm_xor_pd(XMM1, XMM6);
#endif
			XMM0 = _mm_set1_pd(1.0);
			XMM4 = _mm_add_pd(XMM1, XMM2);
			XMM4 = _mm_add_pd(XMM4, XMM3);
			XMM0 = _mm_div_pd(XMM0, XMM4);
			XMM1 = _mm_mul_pd(XMM1, XMM0);
			XMM2 = _mm_mul_pd(XMM2, XMM0);
			XMM3 = _mm_mul_pd(XMM3, XMM0);
			XMM11 = _mm_cvtpd_ps(XMM3);
			XMM12 = _mm_cvtpd_ps(XMM2);
			XMM13 = _mm_cvtpd_ps(XMM1);
			XMM11 = _mm_shuffle_ps(XMM11, XMM11, _MM_SHUFFLE(1, 0, 1, 0));
			XMM12 = _mm_shuffle_ps(XMM12, XMM12, _MM_SHUFFLE(1, 0, 1, 0));
			XMM13 = _mm_shuffle_ps(XMM13, XMM13, _MM_SHUFFLE(1, 0, 1, 0));
#else
			GIdouble v1[2], v2[2];
			GIdouble dA, dB, dC, d1Sum;
			p[0] = (GIdouble)x;
			GI_VEC2_SUB(v1, param[1], p);
			GI_VEC2_SUB(v2, param[2], p);
			dA = GI_VEC2_DET(v1, v2);
			GI_VEC2_SUB(v1, param[2], p);
			GI_VEC2_SUB(v2, param[0], p);
			dB = GI_VEC2_DET(v1, v2);
			GI_VEC2_SUB(v1, param[0], p);
			GI_VEC2_SUB(v2, param[1], p);
			dC = GI_VEC2_DET(v1, v2);
			d1Sum = 1.0 / (dA+dB+dC);
			dA *= d1Sum; dB *= d1Sum; dC *= d1Sum;
#if OPENGI_SSE >= 1
			XMM11 = _mm_set1_ps(dA);
			XMM12 = _mm_set1_ps(dB);
			XMM13 = _mm_set1_ps(dC);
#endif
#endif

			/* interpolate and render attributes */
			for(i=0; i<data->num_attribs; ++i)
			{
#if OPENGI_SSE >= 1
				/* interpolate attribute */
				__m128 XMM0 = _mm_load_ps(temp+(i*12));
				__m128 XMM1 = _mm_load_ps(temp+(i*12)+4);
				__m128 XMM2 = _mm_load_ps(temp+(i*12)+8);
				XMM0 = _mm_mul_ps(XMM0, XMM11);
				XMM1 = _mm_mul_ps(XMM1, XMM12);
				XMM2 = _mm_mul_ps(XMM2, XMM13);
				XMM0 = _mm_add_ps(XMM0, XMM1);
				XMM0 = _mm_add_ps(XMM0, XMM2);

				/* normalize if neccessary */
				a = data->attribs[i];
				if(pSampler->attrib_mode[a] == GI_SAMPLE_NORMALIZED)
				{
					static const GIuint ones = 0xFFFFFFFF;
					__m128 XMM3;
					XMM1 = _mm_set_ss(*(const GIfloat*)&ones);
					XMM1 = _mm_shuffle_ps(XMM1, XMM1, _MM_SHUFFLE(1, 0, 0, 0));
					XMM0 = _mm_and_ps(XMM0, XMM1);
					XMM1 = _mm_mul_ps(XMM0, XMM0);
#if OPENGI_SSE >= 3
					XMM1 = _mm_hadd_ps(XMM1, XMM1);
					XMM1 = _mm_hadd_ps(XMM1, XMM1);
#else
					XMM2 = _mm_shuffle_ps(XMM1, XMM1, _MM_SHUFFLE(0, 1, 2, 3));
					XMM1 = _mm_add_ps(XMM1, XMM2);
					XMM2 = _mm_shuffle_ps(XMM1, XMM1, _MM_SHUFFLE(2, 3, 0, 1));
					XMM1 = _mm_add_ps(XMM1, XMM2);
#endif
					XMM2 = _mm_rsqrt_ps(XMM1);
					XMM1 = _mm_mul_ps(XMM1, XMM2);
					XMM1 = _mm_mul_ps(XMM1, XMM2);
					XMM3 = _mm_set1_ps(3.0f);
					XMM3 = _mm_sub_ps(XMM3, XMM1);
					XMM1 = _mm_set1_ps(0.5f);
					XMM2 = _mm_mul_ps(XMM2, XMM1);
					XMM2 = _mm_mul_ps(XMM2, XMM3);
					XMM0 = _mm_mul_ps(XMM0, XMM2);
//					XMM1 = _mm_sqrt_ps(XMM1);
//					XMM0 = _mm_div_ps(XMM0, XMM1);
					if(pContext->attrib_image[a]->type == GI_UNSIGNED_BYTE)
					{
//						XMM1 = _mm_set1_ps(0.5f);
						XMM0 = _mm_mul_ps(XMM0, XMM1);
						XMM0 = _mm_add_ps(XMM0, XMM1);
					}
				}
				XOUT = XMM0;
				if(pSampler->attrib_mode[a] == GI_SAMPLE_TEXTURED)
					data->textures[a].filter(data->textures+a, 
						(const GIfloat*)&XOUT, (GIfloat*)&XOUT);
				else if(pSampler->attrib_mode[a] == GI_SAMPLE_NORMALIZED)
					((GIfloat*)&XOUT)[3] = 1.0f;
				data->setpixel_fn[a](pContext->attrib_image[a], x, y, (const GIfloat*)&XOUT);
#else
				/* interpolate attribute */
				GIfloat out[4];
				GIfloat *a0 = temp + (i*12), *a1 = a0 + 4, *a2 = a0 + 8;
				out[0] = dA*a0[0] + dB*a1[0] + dC*a2[0];
				out[1] = dA*a0[1] + dB*a1[1] + dC*a2[1];
				out[2] = dA*a0[2] + dB*a1[2] + dC*a2[2];
				out[3] = dA*a0[3] + dB*a1[3] + dC*a2[3];

				/* normalize if neccessary */
				if(pSampler->attrib_mode[a] == GI_SAMPLE_NORMALIZED)
				{
					GIvec3f_normalize(out);
					if(pContext->attrib_image[a]->type == GI_UNSIGNED_BYTE)
					{
						GI_VEC3_ADD_SCALED(out, vcHalf, out, 0.5f);
					}
					out[3] = 1.0f;
				}
				else if(pSampler->attrib_mode[a] == GI_SAMPLE_TEXTURED)
					data->textures[a].filter(data->textures+a, out, out);
				data->setpixel_fn[a](pContext->attrib_image[a], x, y, out);
#endif
			}
		}

		/* advance active edge data */
		++y;
		if(!(--lines[e1]) && y < yMax)
		{
			e1 = ((e1+e2)<<1) % 3;
			dX1 = minX[e1];
		}
		else
			dX1 += incX[e1];
		if(!(--lines[e2]) && y < yMax)
		{
			e2 = ((e1+e2)<<1) % 3;
			dX2 = minX[e2];
		}
		else
			dX2 += incX[e2];
	}
}

/** \internal
 *  \brief Thread execution function for rasterizing triangles.
 *  \param arg thread parameters
 *  \return number of rasterized faces
 *  \ingroup sampling
 */
GIthreadret GITHREADENTRY GISampler_rasterize_thread(void *arg)
{
	GIRasterizerData *pData = (GIRasterizerData*)arg;
	GIFace *pFace;
	GIfloat *pTemp = GI_MALLOC_ALIGNED(GI_SSE_SIZE(
		pData->num_attribs*12*sizeof(GIfloat)), GI_SSE_ALIGN_FLOAT);
	GIuint uiCount = 0;

#if OPENGI_NUM_THREADS > 1
	/* rasterize triangles */
	for(;;)
	{
		/* fetch next face */
		GIMutex_lock(&pData->mutex);
		pFace = pData->next_face;
		if(pFace == pData->end_face)
		{
			GIMutex_unlock(&pData->mutex);
			break;
		}
		pData->next_face = pFace->next;
		GIMutex_unlock(&pData->mutex);

		/* sample face */
		GISampler_rasterize_triangle(pFace, pData, pTemp);
		++uiCount;
	}
#endif

	GI_FREE_ALIGNED(pTemp);
	return (GIthreadret)uiCount;
}

/** \internal
 *  \brief Texture constructor.
 *  \param texture texture to construct
 *  \param gl_texture OpenGL texture object to construct from
 *  \param dim texture dimension in [1,4]
 *  \ingroup sampling
 */
void GITexture_construct(GITexture *texture, GIuint gl_texture, GIuint dim)
{
	static const GIenum target[4] = { GL_TEXTURE_1D, 
		GL_TEXTURE_2D, GL_TEXTURE_3D, GL_TEXTURE_CUBE_MAP };
	static const GIenum binding[4] = { GL_TEXTURE_BINDING_1D, 
		GL_TEXTURE_BINDING_2D, GL_TEXTURE_BINDING_3D, GL_TEXTURE_BINDING_CUBE_MAP };
	static const GIenum image[4] = { GL_TEXTURE_1D, GL_TEXTURE_2D, 
		GL_TEXTURE_3D, GL_TEXTURE_CUBE_MAP_POSITIVE_X };
	GLint iBoundTex, iSize = 0, iLuminance = 0, 
		iIntensity = 0, i, j, iNumImages = (dim<4) ? 1 : 6;
	GLint iImageSize[6];
	GIfloat *pData;

	/* retrieve image data */
	memset(texture, 0, sizeof(GITexture));
	texture->dimension = dim;
	glGetIntegerv(binding[dim-1], &iBoundTex);
	glBindTexture(target[dim-1], gl_texture);
	for(i=0; i<iNumImages; ++i)
	{
		glGetTexLevelParameteriv(image[dim-1]+i, 0, GL_TEXTURE_WIDTH, texture->width+i);
		glGetTexLevelParameteriv(image[dim-1]+i, 0, GL_TEXTURE_HEIGHT, texture->height+i);
		glGetTexLevelParameteriv(image[dim-1]+i, 0, GL_TEXTURE_DEPTH, texture->depth+i);
		texture->image[i] = (GIfloat*)iSize;
		iSize += iImageSize[i] = texture->width[i] * 
			texture->height[i] * texture->depth[i] * 4;
	}
	texture->data = (GIfloat*)GI_MALLOC_ALIGNED(GI_SSE_SIZE(
		iSize*4*sizeof(GIfloat)), GI_SSE_ALIGN_FLOAT);
	for(i=0; i<iNumImages; ++i)
	{
		texture->image[i] = texture->data + (GIint)texture->image[i];
		glGetTexImage(image[dim-1]+i, 0, GL_RGBA, GL_FLOAT, texture->image[i]);
		glGetTexLevelParameteriv(image[dim-1]+i, 0, GL_TEXTURE_LUMINANCE_SIZE, &iLuminance);
		glGetTexLevelParameteriv(image[dim-1]+i, 0, GL_TEXTURE_INTENSITY_SIZE, &iIntensity);
		if(iLuminance)
		{
			for(j=0,pData=texture->image[i]; j<iImageSize[i]; j+=4)
			{
#if OPENGI_SSE >= 1
				__m128 XMM0 = _mm_load_ps(pData+j);
				XMM0 = _mm_shuffle_ps(XMM0, XMM0, _MM_SHUFFLE(3, 0, 0, 0));
				_mm_store_ps(pData+j, XMM0);
#else
				pData[j+1] = pData[j+2] = pData[j];
#endif
			}
		}
		else if(iIntensity)
		{
			for(j=0,pData=texture->image[i]; j<iImageSize[i]; j+=4)
			{
#if OPENGI_SSE >= 1
				__m128 XMM0 = _mm_load1_ps(pData+j);
				_mm_store_ps(pData+j, XMM0);
#else
				pData[j+1] = pData[j+2] = 
					pData[j+3] = pData[j];
#endif
			}
		}
	}

	/* read filter modes */
	glGetTexParameteriv(target[dim-1], GL_TEXTURE_MAG_FILTER, &texture->mag_filter);
	glGetTexParameteriv(target[dim-1], GL_TEXTURE_MIN_FILTER, &texture->min_filter);
	if(texture->min_filter == GL_NEAREST || 
	   texture->min_filter == GL_NEAREST_MIPMAP_NEAREST || 
	   texture->min_filter == GL_NEAREST_MIPMAP_LINEAR)
		texture->min_filter = GL_NEAREST;
	else
		texture->min_filter = GL_LINEAR;
	texture->filter_mode = texture->min_filter==
		texture->mag_filter ? texture->min_filter : 0;

	/* read wrapping modes */
	glGetTexParameteriv(target[dim-1], GL_TEXTURE_WRAP_S, &texture->wrap_s);
	glGetTexParameteriv(target[dim-1], GL_TEXTURE_WRAP_T, &texture->wrap_t);
	glGetTexParameteriv(target[dim-1], GL_TEXTURE_WRAP_R, &texture->wrap_r);
	if(dim == 4 || texture->wrap_s != GL_REPEAT)
		texture->wrap_s = GL_CLAMP;
	if(dim == 4 || texture->wrap_t != GL_REPEAT)
		texture->wrap_t = GL_CLAMP;
	if(dim == 4 || texture->wrap_r != GL_REPEAT)
		texture->wrap_r = GL_CLAMP;
	glBindTexture(target[dim-1], iBoundTex);

	/* set filtering function */
	switch(dim)
	{
	case 1:
		texture->filter = &GITexture_filter1D;
		break;
	case 2:
		texture->filter = &GITexture_filter2D;
		break;
	case 3:
		texture->filter = &GITexture_filter3D;
		break;
	case 4:
		texture->filter = &GITexture_filterCube;
		break;
	default:
		texture->filter = NULL;
	}
}

/** \internal
 *  \brief Texture destructor.
 *  \param texture texture to destruct
 *  \ingroup sampling
 */
void GITexture_destruct(GITexture *texture)
{
	/* clean up */
	if(texture->data)
		GI_FREE_ALIGNED(texture->data);
	memset(texture, 0, sizeof(GITexture));
}

/** \internal
 *  \brief Read pixel from 1D texture.
 *  \param texture texture to read
 *  \param x x-coordinate of pixel
 *  \return address of pixel
 *  \ingroup sampling
 */
const GIfloat* GITexture_getpixel1D(GITexture *texture, GIint x)
{
	/* wrapping modes */
	if(texture->wrap_s == GL_REPEAT)
		x %= texture->width[0];
	else
		x = GI_CLAMP(x, 0, texture->width[0]-1);
	return texture->data + (x<<2);
}

/** \internal
 *  \brief Read pixel from 2D texture.
 *  \param texture texture to read
 *  \param x x-coordinate of pixel
 *  \param y y-coordinate of pixel
 *  \return address of pixel
 *  \ingroup sampling
 */
const GIfloat* GITexture_getpixel2D(GITexture *texture, GIint x, GIint y)
{
	/* wrapping modes */
	if(texture->wrap_s == GL_REPEAT)
		x %= texture->width[0];
	else
		x = GI_CLAMP(x, 0, texture->width[0]-1);
	if(texture->wrap_t == GL_REPEAT)
		y %= texture->height[0];
	else
		y = GI_CLAMP(y, 0, texture->height[0]-1);
	return texture->data + ((y*texture->width[0]+x)<<2);
}

/** \internal
 *  \brief Read pixel from 3D texture.
 *  \param texture texture to read
 *  \param x x-coordinate of pixel
 *  \param y y-coordinate of pixel
 *  \param z z-coordinate of pixel
 *  \return address of pixel
 *  \ingroup sampling
 */
const GIfloat* GITexture_getpixel3D(GITexture *texture, 
									GIint x, GIint y, GIint z)
{
	/* wrapping modes */
	if(texture->wrap_s == GL_REPEAT)
		x %= texture->width[0];
	else
		x = GI_CLAMP(x, 0, texture->width[0]-1);
	if(texture->wrap_t == GL_REPEAT)
		y %= texture->height[0];
	else
		y = GI_CLAMP(y, 0, texture->height[0]-1);
	if(texture->wrap_r == GL_REPEAT)
		z %= texture->depth[0];
	else
		z = GI_CLAMP(z, 0, texture->depth[0]-1);
	return texture->data + (((z*texture->height[0]+y)*texture->width[0]+x)<<2);
}

/** \internal
 *  \brief Read pixel from cube map side.
 *  \param texture texture to read
 *  \param face face to read from
 *  \param x x-coordinate of pixel
 *  \param y y-coordinate of pixel
 *  \return address of pixel
 *  \ingroup sampling
 */
const GIfloat* GITexture_getpixelCube(GITexture *texture, 
									  GIuint face, GIint x, GIint y)
{
	/* wrapping modes */
	x = GI_CLAMP(x, 0, texture->width[face]-1);
	y = GI_CLAMP(y, 0, texture->height[face]-1);
	return texture->image[face] + ((y*texture->width[face]+x)<<2);
}

/** \internal
 *  \brief Sample 1D texture at given texture coordinate.
 *  \param texture texture to sample
 *  \param coords texture coordinates
 *  \param out vector to place results in
 *  \ingroup sampling
 */
void GITexture_filter1D(GITexture *texture, const GIfloat *coords, GIfloat *out)
{
	GIfloat s = coords[0]*(GIfloat)texture->width[0] - 0.5f;

	/* make positive if repeating */
	if(texture->wrap_s == GL_REPEAT)
		while(s < 0.0f)
			s += (GIfloat)texture->width[0];

	/* filter */
	if(texture->filter_mode == GL_NEAREST)
	{
#if OPENGI_SSE >= 1
		__m128 XMM0 = _mm_load_ps(GITexture_getpixel1D(texture, GI_ROUND(s)));
		_mm_store_ps(out, XMM0);
#else
		const GIfloat *pixel = GITexture_getpixel1D(texture, GI_ROUND(s));
		out[0] = pixel[0];
		out[1] = pixel[1];
		out[2] = pixel[2];
		out[3] = pixel[3];
#endif
	}
	else
	{
		GIint x = floor(s);
#if OPENGI_SSE >= 1
		__m128 XMM0 = _mm_set1_ps(1.0f);
		__m128 XMM1 = _mm_set1_ps(s-(GIfloat)x);
		__m128 XMM2 = _mm_load_ps(GITexture_getpixel1D(texture, x));
		__m128 XMM3 = _mm_load_ps(GITexture_getpixel1D(texture, x+1));
		XMM0 = _mm_sub_ps(XMM0, XMM1);
		XMM2 = _mm_mul_ps(XMM2, XMM0);
		XMM3 = _mm_mul_ps(XMM3, XMM1);
		XMM2 = _mm_add_ps(XMM2, XMM3);
		_mm_store_ps(out, XMM2);
#else
		GIfloat f1 = s - (GIfloat)x;
		GIfloat f0 = 1.0f - f1;
		const GIfloat *t0 = GITexture_getpixel1D(texture, x);
		const GIfloat *t1 = GITexture_getpixel1D(texture, x+1);
		out[0] = f0*t0[0] + f1*t1[0];
		out[1] = f0*t0[1] + f1*t1[1];
		out[2] = f0*t0[2] + f1*t1[2];
		out[3] = f0*t0[3] + f1*t1[3];
#endif
	}
}

/** \internal
 *  \brief Sample 2D texture at given texture coordinate.
 *  \param texture texture to sample
 *  \param coords texture coordinates
 *  \param out vector to place results in
 *  \ingroup sampling
 */
void GITexture_filter2D(GITexture *texture, const GIfloat *coords, GIfloat *out)
{
	GIfloat s = coords[0]*(GIfloat)texture->width[0] - 0.5f;
	GIfloat t = coords[1]*(GIfloat)texture->height[0] - 0.5f;

	/* make positive if repeating */
	if(texture->wrap_s == GL_REPEAT)
		while(s < 0.0f)
			s += (GIfloat)texture->width[0];
	if(texture->wrap_t == GL_REPEAT)
		while(t < 0.0f)
			t += (GIfloat)texture->height[0];

	/* filter */
	if(texture->filter_mode == GL_NEAREST)
	{
#if OPENGI_SSE >= 1
		__m128 XMM0 = _mm_load_ps(GITexture_getpixel2D(
			texture, GI_ROUND(s), GI_ROUND(t)));
		_mm_store_ps(out, XMM0);
#else
		const GIfloat *pixel = GITexture_getpixel2D(
			texture, GI_ROUND(s), GI_ROUND(t));
		out[0] = pixel[0];
		out[1] = pixel[1];
		out[2] = pixel[2];
		out[3] = pixel[3];
#endif
	}
	else
	{
		GIint x = floor(s), y = floor(t);
		GIfloat fRatioS = s - (GIfloat)x, fRatioT = t - (GIfloat)y;
		GIfloat f1RatioS = 1.0f - fRatioS, f1RatioT = 1.0f - fRatioT;
#if OPENGI_SSE >= 1
		__m128 XMM0 = _mm_set1_ps(f1RatioS*f1RatioT);
		__m128 XMM1 = _mm_set1_ps(fRatioS*f1RatioT);
		__m128 XMM2 = _mm_set1_ps(f1RatioS*fRatioT);
		__m128 XMM3 = _mm_set1_ps(fRatioS*fRatioT);
		__m128 XMM4 = _mm_load_ps(GITexture_getpixel2D(texture, x, y));
		__m128 XMM5 = _mm_load_ps(GITexture_getpixel2D(texture, x+1, y));
		__m128 XMM6 = _mm_load_ps(GITexture_getpixel2D(texture, x, y+1));
		__m128 XMM7 = _mm_load_ps(GITexture_getpixel2D(texture, x+1, y+1));
		XMM4 = _mm_mul_ps(XMM4, XMM0);
		XMM5 = _mm_mul_ps(XMM5, XMM1);
		XMM6 = _mm_mul_ps(XMM6, XMM2);
		XMM7 = _mm_mul_ps(XMM7, XMM3);
		XMM4 = _mm_add_ps(XMM4, XMM5);
		XMM4 = _mm_add_ps(XMM4, XMM6);
		XMM4 = _mm_add_ps(XMM4, XMM7);
		_mm_store_ps(out, XMM4);
#else
		GIfloat f00 = f1RatioT * f1RatioS, f10 = f1RatioT * fRatioS;
		GIfloat f01 = fRatioT * f1RatioS, f11 = fRatioT * fRatioS;
		const GIfloat *t00 = GITexture_getpixel2D(texture, x, y);
		const GIfloat *t10 = GITexture_getpixel2D(texture, x+1, y);
		const GIfloat *t01 = GITexture_getpixel2D(texture, x, y+1);
		const GIfloat *t11 = GITexture_getpixel2D(texture, x+1, y+1);
		out[0] = f00*t00[0] + f10*t10[0] + f01*t01[0] + f11*t11[0];
		out[1] = f00*t00[1] + f10*t10[1] + f01*t01[1] + f11*t11[1];
		out[2] = f00*t00[2] + f10*t10[2] + f01*t01[2] + f11*t11[2];
		out[3] = f00*t00[3] + f10*t10[3] + f01*t01[3] + f11*t11[3];
#endif
	}
}

/** \internal
 *  \brief Sample 3D texture at given texture coordinate.
 *  \param texture texture to sample
 *  \param coords texture coordinates
 *  \param out vector to place results in
 *  \ingroup sampling
 */
void GITexture_filter3D(GITexture *texture, const GIfloat *coords, GIfloat *out)
{
	GIfloat s = coords[0]*(GIfloat)texture->width[0] - 0.5f;
	GIfloat t = coords[1]*(GIfloat)texture->height[0] - 0.5f;
	GIfloat r = coords[2]*(GIfloat)texture->depth[0] - 0.5f;

	/* make positive if repeating */
	if(texture->wrap_s == GL_REPEAT)
		while(s < 0.0f)
			s += (GIfloat)texture->width[0];
	if(texture->wrap_t == GL_REPEAT)
		while(t < 0.0f)
			t += (GIfloat)texture->height[0];
	if(texture->wrap_r == GL_REPEAT)
		while(r < 0.0f)
			r += (GIfloat)texture->depth[0];

	/* filter */
	if(texture->filter_mode == GL_NEAREST)
	{
#if OPENGI_SSE >= 1
		__m128 XMM0 = _mm_load_ps(GITexture_getpixel3D(
			texture, GI_ROUND(s), GI_ROUND(t), GI_ROUND(r)));
		_mm_store_ps(out, XMM0);
#else
		const GIfloat *pixel = GITexture_getpixel3D(
			texture, GI_ROUND(s), GI_ROUND(t), GI_ROUND(r));
		out[0] = pixel[0];
		out[1] = pixel[1];
		out[2] = pixel[2];
		out[3] = pixel[3];
#endif
	}
	else
	{
		GIint x = floor(s), y = floor(t), z = floor(r);
		GIfloat fRatioS = s-(GIfloat)x, fRatioT = t-(GIfloat)y, fRatioR = r-(GIfloat)z;
		GIfloat f1RatioS = 1.0f-fRatioS, f1RatioT = 1.0f-fRatioT, f1RatioR = 1.0f-fRatioR;
#if OPENGI_SSE >= 1
		__m128 XOUT, XMM0 = _mm_set1_ps(f1RatioT*f1RatioS);
		__m128 XMM1 = _mm_set1_ps(f1RatioT*fRatioS);
		__m128 XMM2 = _mm_set1_ps(fRatioT*f1RatioS);
		__m128 XMM3 = _mm_set1_ps(fRatioT*fRatioS);
		__m128 XMM4 = _mm_load_ps(GITexture_getpixel3D(texture, x, y, z));
		__m128 XMM5 = _mm_load_ps(GITexture_getpixel3D(texture, x+1, y, z));
		__m128 XMM6 = _mm_load_ps(GITexture_getpixel3D(texture, x, y+1, z));
		__m128 XMM7 = _mm_load_ps(GITexture_getpixel3D(texture, x+1, y+1, z));
		XMM4 = _mm_mul_ps(XMM4, XMM0);
		XMM5 = _mm_mul_ps(XMM5, XMM1);
		XMM6 = _mm_mul_ps(XMM6, XMM2);
		XMM7 = _mm_mul_ps(XMM7, XMM3);
		XMM4 = _mm_add_ps(XMM4, XMM5);
		XMM4 = _mm_add_ps(XMM4, XMM6);
		XMM4 = _mm_add_ps(XMM4, XMM7);
		XMM5 = _mm_set1_ps(f1RatioR);
		XOUT = _mm_mul_ps(XMM4, XMM5);
		XMM4 = _mm_load_ps(GITexture_getpixel3D(texture, x, y, z+1));
		XMM5 = _mm_load_ps(GITexture_getpixel3D(texture, x+1, y, z+1));
		XMM6 = _mm_load_ps(GITexture_getpixel3D(texture, x, y+1, z+1));
		XMM7 = _mm_load_ps(GITexture_getpixel3D(texture, x+1, y+1, z+1));
		XMM4 = _mm_mul_ps(XMM4, XMM0);
		XMM5 = _mm_mul_ps(XMM5, XMM1);
		XMM6 = _mm_mul_ps(XMM6, XMM2);
		XMM7 = _mm_mul_ps(XMM7, XMM3);
		XMM4 = _mm_add_ps(XMM4, XMM5);
		XMM4 = _mm_add_ps(XMM4, XMM6);
		XMM4 = _mm_add_ps(XMM4, XMM7);
		XMM5 = _mm_set1_ps(fRatioR);
		XMM4 = _mm_mul_ps(XMM4, XMM5);
		XOUT = _mm_add_ps(XOUT, XMM4);
		_mm_store_ps(out, XOUT);
#else
		GIfloat f00 = f1RatioT*f1RatioS, f10 = f1RatioT*fRatioS;
		GIfloat f01 = fRatioT*f1RatioS, f11 = fRatioT*fRatioS;
		const GIfloat *t000 = GITexture_getpixel3D(texture, x, y, z);
		const GIfloat *t100 = GITexture_getpixel3D(texture, x+1, y, z);
		const GIfloat *t010 = GITexture_getpixel3D(texture, x, y+1, z);
		const GIfloat *t110 = GITexture_getpixel3D(texture, x+1, y+1, z);
		const GIfloat *t001 = GITexture_getpixel3D(texture, x, y, z+1);
		const GIfloat *t101 = GITexture_getpixel3D(texture, x+1, y, z+1);
		const GIfloat *t011 = GITexture_getpixel3D(texture, x, y+1, z+1);
		const GIfloat *t111 = GITexture_getpixel3D(texture, x+1, y+1, z+1);
		out[0] = f1RatioR*(f00*t000[0]+f10*t100[0]+f01*t010[0]+f11*t110[0]) + 
			fRatioR*(f00*t001[0]+f10*t101[0]+f01*t011[0]+f11*t111[0]);
		out[1] = f1RatioR*(f00*t000[1]+f10*t100[1]+f01*t010[1]+f11*t110[1]) + 
			fRatioR*(f00*t001[1]+f10*t101[1]+f01*t011[1]+f11*t111[1]);
		out[2] = f1RatioR*(f00*t000[2]+f10*t100[2]+f01*t010[2]+f11*t110[2]) + 
			fRatioR*(f00*t001[2]+f10*t101[2]+f01*t011[2]+f11*t111[2]);
		out[3] = f1RatioR*(f00*t000[3]+f10*t100[3]+f01*t010[3]+f11*t110[3]) + 
			fRatioR*(f00*t001[3]+f10*t101[3]+f01*t011[3]+f11*t111[3]);
#endif
	}
}

/** \internal
 *  \brief Sample cube map texture at given texture coordinate.
 *  \param texture texture to sample
 *  \param coords texture coordinates
 *  \param out vector to place results in
 *  \ingroup sampling
 */
void GITexture_filterCube(GITexture *texture, const GIfloat *coords, GIfloat *out)
{
	GIfloat ax = fabs(coords[0]), ay = fabs(coords[1]), az = fabs(coords[2]);
	GIfloat s, t, m;
	GIint i;

	/* compute face and coordinates in face */
	if(ay > ax)
		i = (az>ay) ? 4 : 2;
	else
		i = (az>ax) ? 4 : 0;
	if(coords[i>>1] < 0.0)
		i += 1;
	switch(i)
	{
	case 0: s = -coords[2]; t = -coords[1]; m = 1.0f / ax; break;
	case 1: s = coords[2]; t = -coords[1]; m = 1.0f / ax; break;
	case 2: s = coords[0]; t = coords[2]; m = 1.0f / ay; break;
	case 3: s = coords[0]; t = -coords[2]; m = 1.0f / ay; break;
	case 4: s = coords[0]; t = -coords[1]; m = 1.0f / az; break;
	case 5: s = -coords[0]; t = -coords[1]; m = 1.0f / az; break;
	}
	s = 0.5f*(s*m+1.0f)*(GIfloat)texture->width[i] - 0.5f;
	t = 0.5f*(t*m+1.0f)*(GIfloat)texture->height[i] - 0.5f;

	/* filter */
	if(texture->filter_mode == GL_NEAREST)
	{
#if OPENGI_SSE >= 1
		__m128 XMM0 = _mm_load_ps(GITexture_getpixelCube(
			texture, i, GI_ROUND(s), GI_ROUND(t)));
		_mm_store_ps(out, XMM0);
#else
		const GIfloat *pixel = GITexture_getpixelCube(
			texture, i, GI_ROUND(s), GI_ROUND(t));
		out[0] = pixel[0];
		out[1] = pixel[1];
		out[2] = pixel[2];
		out[3] = pixel[3];
#endif
	}
	else
	{
		GIint x = floor(s), y = floor(t);
		GIfloat fRatioS = s - (GIfloat)x, fRatioT = t - (GIfloat)y;
		GIfloat f1RatioS = 1.0f - fRatioS, f1RatioT = 1.0f - fRatioT;
#if OPENGI_SSE >= 1
		__m128 XMM0 = _mm_set1_ps(f1RatioS*f1RatioT);
		__m128 XMM1 = _mm_set1_ps(fRatioS*f1RatioT);
		__m128 XMM2 = _mm_set1_ps(f1RatioS*fRatioT);
		__m128 XMM3 = _mm_set1_ps(fRatioS*fRatioT);
		__m128 XMM4 = _mm_load_ps(GITexture_getpixelCube(texture, i, x, y));
		__m128 XMM5 = _mm_load_ps(GITexture_getpixelCube(texture, i, x+1, y));
		__m128 XMM6 = _mm_load_ps(GITexture_getpixelCube(texture, i, x, y+1));
		__m128 XMM7 = _mm_load_ps(GITexture_getpixelCube(texture, i, x+1, y+1));
		XMM4 = _mm_mul_ps(XMM4, XMM0);
		XMM5 = _mm_mul_ps(XMM5, XMM1);
		XMM6 = _mm_mul_ps(XMM6, XMM2);
		XMM7 = _mm_mul_ps(XMM7, XMM3);
		XMM4 = _mm_add_ps(XMM4, XMM5);
		XMM4 = _mm_add_ps(XMM4, XMM6);
		XMM4 = _mm_add_ps(XMM4, XMM7);
		_mm_store_ps(out, XMM4);
#else
		GIfloat f00 = f1RatioT * f1RatioS, f10 = f1RatioT * fRatioS;
		GIfloat f01 = fRatioT * f1RatioS, f11 = fRatioT * fRatioS;
		const GIfloat *t00 = GITexture_getpixelCube(texture, i, x, y);
		const GIfloat *t10 = GITexture_getpixelCube(texture, i, x+1, y);
		const GIfloat *t01 = GITexture_getpixelCube(texture, i, x, y+1);
		const GIfloat *t11 = GITexture_getpixelCube(texture, i, x+1, y+1);
		out[0] = f00*t00[0] + f10*t10[0] + f01*t01[0] + f11*t11[0];
		out[1] = f00*t00[1] + f10*t10[1] + f01*t01[1] + f11*t11[1];
		out[2] = f00*t00[2] + f10*t10[2] + f01*t01[2] + f11*t11[2];
		out[3] = f00*t00[3] + f10*t10[3] + f01*t01[3] + f11*t11[3];
#endif
	}
}
#endif