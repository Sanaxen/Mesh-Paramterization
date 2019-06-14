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
 *  \brief Implementation of structures and functions for context handling.
 */

#include "gi_context.h"
#include "gi_memory.h"

#include <stdlib.h>
#include <stdio.h>

#include "gi_blas.h"


/** \internal
 *  \brief Current OpenGI context.
 */
//#ifdef _DEBUG
	static GIContext *gi_CurrentContext = NULL;
//#else
//	GIContext *gi_CurrentContext = NULL;
//#endif

/** Create new context.
 *  \return created GI context
 *  \ingroup context
 */
GIcontext GIAPIENTRY giCreateContext()
{
	GIContext *pContext;
	GIuint i;

	/* create allocator if neccessary */
	if(!g_SmallObjAlloc.pool)
		GISmallObjectAllocator_construct(&g_SmallObjAlloc);
	if(!g_PersistentAlloc.pool)
		GISmallObjectAllocator_construct(&g_PersistentAlloc);

	/* create new context and initialize states */
	pContext = (GIContext*)GI_CALLOC_SINGLE(sizeof(GIContext));
	GIHash_construct(&pContext->mesh_hash, 32, 0.0f, sizeof(GIuint), 
		hash_uint, compare_uint, copy_uint);
	GIHash_construct(&pContext->image_hash, 64, 0.0f, sizeof(GIuint), 
		hash_uint, compare_uint, copy_uint);
	pContext->next_mid = 1;
	pContext->next_iid = 1;
	pContext->use_threads = GI_TRUE;
	pContext->error = GI_NO_ERROR;
	pContext->semantic[0] = 0;
	pContext->attrib_semantic[0] = GI_POSITION_ATTRIB;
	for(i=1; i<GI_SEMANTIC_COUNT; ++i)
	{
		pContext->semantic[i] = GI_ATTRIB_COUNT - GI_SEMANTIC_COUNT + i;
		pContext->attrib_semantic[pContext->semantic[i]] = i + GI_SEMANTIC_BASE;
	}
	GICutter_construct(&pContext->cutter, pContext);
	GIParameterizer_construct(&pContext->parameterizer, pContext);
//	GISampler_construct(&pContext->sampler, pContext);
	GIRenderer_construct(&pContext->renderer, pContext);
	pContext->gl_manager = (GIGLManager*)GI_CALLOC_SINGLE(sizeof(GIGLManager));
	return pContext;
}

/** Set active Context.
 *  \param context GI context to make current
 *  \ingroup context
 */
void GIAPIENTRY giMakeCurrent(GIcontext context)
{
	/* set context */
	gi_CurrentContext = (GIContext*)context;
}

/** Get active context.
 *  \return active GI context
 *  \ingroup context
 */
GIcontext GIAPIENTRY giGetCurrent()
{
	/* return context */
	return gi_CurrentContext;
}

/** Delete context and all associated data.
 *  \param context GI context to destroy
 *  \ingroup context
 */
void GIAPIENTRY giDestroyContext(GIcontext context)
{
	GIContext *pContext = (GIContext*)context;
	GIMesh *pMesh;
	GIuint i;

	/* clean up */
	for(i=1; i<pContext->next_mid; ++i)
	{
		pMesh = (GIMesh*)GIHash_remove(&pContext->mesh_hash, &i);
		if(pMesh)
		{
			GIMesh_destruct(pMesh);
			GI_FREE_SINGLE(pMesh, sizeof(GIMesh));
		}
	}
	GIHash_destruct(&pContext->mesh_hash, 0);
//	GIHash_destruct(&pContext->image_hash, sizeof(GIImage));
	GIGLManager_destruct(pContext->gl_manager);
	GI_FREE_SINGLE(pContext->gl_manager, sizeof(GIGLManager));

	/* unbind if current and delete */
	if(gi_CurrentContext == pContext)
		gi_CurrentContext = NULL;
	GI_FREE_SINGLE(pContext, sizeof(GIContext));
}

/** Retrieve boolean state value.
 *  \param pname state to query
 *  \param params address of variable to store value
 *  \ingroup state
 */
void GIAPIENTRY giGetBooleanv(GIenum pname, GIboolean *params)
{
	GIContext *pContext = GIContext_current();

	/* select state and get value */
	switch(pname)
	{
	case GI_EXACT_MAPPING_SUBSET_SORTED:
	case GI_PARAM_CORNER_SUBSET_SORTED:
		*params = pContext->subset_sorted[pname-GI_EXACT_MAPPING_SUBSET_SORTED];
		break;
/*	case GI_STRAIGHTEN_CUT:
		*params = pContext->cutter.straighten;
		break;
*/	case GI_SAMPLER_USE_SHADER:
//		*params = pContext->sampler.use_shader;
		break;
	case GI_SAMPLER_USE_RENDER_TO_TEXTURE:
//		*params = pContext->sampler.use_fbo;
		break;
	case GI_GL_USE_VERTEX_TEXTURE:
	case GI_GL_USE_GEOMETRY_SHADER:
		*params = ((pContext->renderer.gim_flags&pname) != 0);
		break;
	default:
		*params = giIsEnabled(pname);
	}
}

/** Retrieve integer state value.
 *  \param pname state to query
 *  \param params address of variable to store value
 *  \ingroup state
 */
void GIAPIENTRY giGetIntegerv(GIenum pname, GIint *params)
{
	GIContext *pContext = GIContext_current();

	/* select state and get value */
	switch(pname)
	{
	case GI_VERSION:
		params[0] = GI_VERSION_MAJOR;
		params[1] = GI_VERSION_MINOR;
		params[2] = GI_VERSION_REVISION;
		break;
	case GI_MAX_ATTRIBS:
		*params = GI_ATTRIB_COUNT;
		break;
	case GI_MESH_BINDING:
		*params = pContext->mesh ? pContext->mesh->id : 0;
		break;
	case GI_IMAGE_BINDING:
//		*params = pContext->image ? pContext->image->id : 0;
		break;
	case GI_EXACT_MAPPING_SUBSET_COUNT:
	case GI_PARAM_CORNER_SUBSET_COUNT:
		*params = pContext->subset_count[pname-GI_EXACT_MAPPING_SUBSET_COUNT];
		break;
	case GI_POSITION_ATTRIB:
	case GI_PARAM_ATTRIB:
	case GI_PARAM_STRETCH_ATTRIB:
		*params = pContext->semantic[pname-GI_SEMANTIC_BASE];
		break;
	case GI_CUTTER:
		*params = pContext->cutter.cutter;
		break;
	case GI_SUBDIVISION_ITERATIONS:
		*params = pContext->cutter.iterations;
		break;
	case GI_PARAMETERIZER:
		*params = pContext->parameterizer.parameterizer;
		break;
	case GI_INITIAL_PARAMETERIZATION:
		*params = pContext->parameterizer.initial_param;
		break;
	case GI_STRETCH_METRIC:
		*params = pContext->parameterizer.stretch_metric;
		break;
	case GI_PARAM_RESOLUTION:
		*params = pContext->parameterizer.sampling_res;
		break;
	case GI_UNSYMMETRIC_SOLVER:
		*params = pContext->parameterizer.solver;
		break;
	case GI_PARAM_SOURCE_ATTRIB:
		*params = pContext->parameterizer.source_attrib;
		break;
	case GI_SAMPLER:
//		*params = pContext->sampler.sampler;
		break;
	case GI_SAMPLED_ATTRIB_COUNT:
		{
			//GIint i, c = 0;
			//for(i=0; i<GI_ATTRIB_COUNT; ++i)
			//	if(pContext->sampler.sampled_attribs & (1<<i))
			//		++c;
			//*params = c;
		}
		break;
	case GI_SAMPLED_ATTRIBS:
		//*params = pContext->sampler.sampled_attribs;
		break;
	case GI_RENDER_RESOLUTION_U:
	case GI_RENDER_RESOLUTION_V:
		*params = pContext->renderer.render_res[pname-
			GI_RENDER_RESOLUTION_U];
		break;
	case GI_RENDER_CACHE_SIZE:
		*params = pContext->renderer.gim_cache_size;
		break;
	default:
		{
			GIboolean bBool = GI_FALSE;
			giGetBooleanv(pname, &bBool);
			*params = bBool;
		}
	}
}

/** Retrieve floating point state value.
 *  \param pname state to query
 *  \param params address of variable to store value
 *  \ingroup state
 */
void GIAPIENTRY giGetFloatv(GIenum pname, GIfloat *params)
{
	GIContext *pContext = GIContext_current();

	/* select state and get value */
	switch(pname)
	{
	case GI_VERSION:
		params[0] = GI_VERSION_MAJOR;
		params[1] = GI_VERSION_MINOR;
		params[2] = GI_VERSION_REVISION;
		break;
	case GI_CONFORMAL_WEIGHT:
		*params = pContext->parameterizer.conformal_weight;
		break;
	case GI_AUTHALIC_WEIGHT:
		*params = pContext->parameterizer.authalic_weight;
		break;
	case GI_STRETCH_WEIGHT:
		*params = pContext->parameterizer.stretch_weight;
		break;
	case GI_AREA_WEIGHT:
		*params = pContext->parameterizer.area_weight;
		break;
/*	case GI_ORIENTATION_WEIGHT:
		*params = pContext->cutter.orientation_weight;
		break;
	case GI_SHAPE_WEIGHT:
		*params = pContext->cutter.shape_weight;
		break;
*/	default:
		{
			GIint iInt = 0;
			giGetIntegerv(pname, &iInt);
			*params = iInt;
		}
	}
}

/** \internal
 *  \brief Enable/disable feature
 *  \param pname feature to enable/disable
 *  \param enable enable flag or -1 for query only
 *  \retval GI_TRUE feature is enabled
 *  \retval GI_FALSE feature is disbaled
 *  \ingroup state
 */
GIboolean giEnableDisable(GIenum pname, GIbyte enable)
{
	GIContext *pContext = GIContext_current();

	/* select state and set value */
	switch(pname)
	{
	case GI_MULTITHREADING:
		if(enable >= 0)
			pContext->use_threads = enable;
		return pContext->use_threads;
	case GI_EXACT_MAPPING_SUBSET:
	case GI_PARAM_CORNER_SUBSET:
		if(enable >= 0)
			pContext->subset_enabled[pname-GI_SUBSET_BASE] = enable;
		return pContext->subset_enabled[pname-GI_SUBSET_BASE];
	default:
		GIContext_error(pContext, GI_INVALID_ENUM);
		return GI_FALSE;
	}
}

/** Enable feature.
 *  \param pname feature to enable
 *  \ingroup state
 */
void GIAPIENTRY giEnable(GIenum pname)
{
	/* set state */
	giEnableDisable(pname, GI_TRUE);
}

/** Disable feture.
 *  \param pname feature to disable
 *  \ingroup state
 */
void GIAPIENTRY giDisable(GIenum pname)
{
	/* set state */
	giEnableDisable(pname, GI_FALSE);
}

/** Check if feature is enabled.
 *  \param pname feature to check
 *  \retval GI_TRUE feature is enabled
 *  \retval GI_FALSE feature is disbaled
 *  \ingroup state
 */
GIboolean GIAPIENTRY giIsEnabled(GIenum pname)
{
	/* query state */
	return giEnableDisable(pname, -1);
}

/** Set vertex subset.
 *  \param subset subset to set
 *  \param count number of vertices ins subset
 *  \param sorted GI_TRUE if indices sorted, GI_FALSE else
 *  \param indices indices of vertices in subset
 *  \ingroup state
 */
void GIAPIENTRY giVertexSubset(GIenum subset, GIsizei count, GIboolean sorted, const GIuint *indices)
{
	GIContext *pContext = GIContext_current();

	/* error checking */
	if(subset < GI_SUBSET_BASE || subset > GI_SUBSET_END)
	{
		GIContext_error(pContext, GI_INVALID_ENUM);
		return;
	}
	if(count <= 0)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}
	if(!indices)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}

	/* set subset data */
	subset -= GI_SUBSET_BASE;
	pContext->subset[subset] = indices;
	pContext->subset_count[subset] = count;
	pContext->subset_sorted[subset] = sorted;
}

/** Get pointer to attribute array or vertex subset.
 *  \param pname array to query
 *  \param params address to store pointer at
 *  \ingroup state
 */
void GIAPIENTRY giGetPointerv(GIenum pname, GIvoid **params)
{
	GIContext *pContext = GIContext_current();

	/* get pointer */
	if(pname >= GI_SUBSET_BASE && pname <= GI_SUBSET_END)
		/**params = pContext->subset[pname-GI_SUBSET_BASE]*/;
	else if(pname == GI_IMAGE_DATA)
	{
		//if(!pContext->image)
			/*GIContext_error(pContext, GI_INVALID_OPERATION)*/;
		//else
			/**params = pContext->image->data*/;
	}
	else
		GIContext_error(pContext, GI_INVALID_ENUM);
}

/** Bind image to attribute channel.
 *  \param attrib attribute to bind image to
 *  \param image image to bind
 */
void GIAPIENTRY giAttribImage(GIuint attrib, GIuint image)
{
#if 0
	GIContext *pContext = GIContext_current();
	GIImage *pImage;

	/* image allready bound */
	if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}
	if(pContext->attrib_image[attrib] && 
		pContext->attrib_image[attrib]->id == image)
		return;
	if(!image)
	{
		pContext->attrib_image[attrib] = NULL;
		return;
	}

	/* search and bind image */
	pImage = (GIImage*)GIHash_find(&pContext->image_hash, &image);
	if(!pImage)
		GIContext_error(pContext, GI_INVALID_ID);
	else
		pContext->attrib_image[attrib] = pImage;
#endif
}

/** Bind attribute to special semantic.
 *  \param semantic semantic to bind attribute to
 *  \param attrib attrib to use channel for semantic
 */
void GIAPIENTRY giBindAttrib(GIenum semantic, GIuint attrib)
{
	GIContext *pContext = GIContext_current();

	if(semantic < GI_SEMANTIC_BASE || semantic > GI_SEMANTIC_END)
		GIContext_error(pContext, GI_INVALID_ENUM);
	else if(attrib >= GI_ATTRIB_COUNT)
		GIContext_error(pContext, GI_INVALID_VALUE);
	else
	{
		pContext->attrib_semantic[pContext->semantic[
			semantic-GI_SEMANTIC_BASE]] = GI_NONE;
		pContext->semantic[semantic-GI_SEMANTIC_BASE] = attrib;
		pContext->attrib_semantic[attrib] = semantic;
	}
}

/** Enable attribute array.
 *  \param attrib attribute to enable array for
 */
void GIAPIENTRY giEnableAttribArray(GIuint attrib)
{
	if(attrib < GI_ATTRIB_COUNT)
		GIContext_current()->attrib_enabled[attrib] = GI_TRUE;
	else
		GIContext_error(GIContext_current(), GI_INVALID_VALUE);
}

/** Disable attribute array.
 *  \param attrib attribute to disable array for
 */
void GIAPIENTRY giDisableAttribArray(GIuint attrib)
{
	if(attrib < GI_ATTRIB_COUNT)
		GIContext_current()->attrib_enabled[attrib] = GI_FALSE;
	else
		GIContext_error(GIContext_current(), GI_INVALID_VALUE);
}

/** Set attribute array to use during mesh creation and retrieval.
 *  This function sets the pointer to the attribute data that is used to 
 *  create and retrieve the mesh. This state is set for the entire context 
 *  and not per mesh as the actual mesh creation and therefore use of the 
 *  array is not done till either giIndexedMesh() or giNonIndexedMesh() is called.
 *  \param attrib attribute index to set pointer for
 *  \param size number of components per attribute
 *  \param normalized normalization flag
 *  \param stride distance between two consecutive vertices or 0 if tightly packed
 *  \param pointer pointer to attribute data
 *  \ingroup state
 */
void GIAPIENTRY giAttribPointer(GIuint attrib, GIint size, 
								GIboolean normalized, GIsizei stride, 
								GIfloat *pointer)
{
	GIContext *pContext = GIContext_current();

	/* set data, if attribute valid */
	if(attrib >= GI_ATTRIB_COUNT || size < 1 || size > 4)
		GIContext_error(pContext, GI_INVALID_VALUE);
	pContext->attrib_pointer[attrib] = pointer;
	pContext->attrib_size[attrib] = size;
	pContext->attrib_stride[attrib] = stride ? stride : size;
	pContext->attrib_normalized[attrib] = normalized;
}

/** Retrieve boolean attribute state.
 *  \param attrib attribute to query
 *  \param pname state to query
 *  \param params address of variable to store value
 *  \ingroup state
 */
void GIAPIENTRY giGetAttribbv(GIuint attrib, GIenum pname, GIboolean *params)
{
	GIContext *pContext = GIContext_current();
	if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}

	/* select state and get value */
	switch(pname)
	{
	case GI_ATTRIB_ARRAY_ENABLED:
		*params = pContext->attrib_enabled[attrib];
		break;
	case GI_ATTRIB_ARRAY_NORMALIZED:
		*params = pContext->attrib_normalized[attrib];
		break;
	default:
		GIContext_error(pContext, GI_INVALID_ENUM);
	}
}

/** Retrieve integer attribute state.
 *  \param attrib attribute to query
 *  \param pname state to query
 *  \param params address of variable to store value
 *  \ingroup state
 */
void GIAPIENTRY giGetAttribiv(GIuint attrib, GIenum pname, GIint *params)
{
	GIContext *pContext = GIContext_current();
	if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}

	/* select state and get value */
	switch(pname)
	{
	case GI_ATTRIB_ARRAY_SIZE:
		*params = pContext->attrib_size[attrib];
		break;
	case GI_ATTRIB_ARRAY_STRIDE:
		*params = pContext->attrib_stride[attrib];
		break;
	case GI_ATTRIB_ARRAY_SEMANTIC:
		*params = pContext->attrib_semantic[attrib];
		break;
	case GI_ATTRIB_IMAGE:
		//*params = pContext->attrib_image[attrib] ? 
		//	pContext->attrib_image[attrib]->id : 0;
		break;
	case GI_GL_RENDER_SEMANTIC:
		*params = pContext->renderer.attrib_semantic[attrib];
		break;
	case GI_GL_RENDER_CHANNEL:
		*params = pContext->renderer.attrib_channel[attrib];
		break;
	case GI_TEXTURE_COORD_DOMAIN:
		*params = pContext->renderer.image_domain[attrib];
		break;
	case GI_SAMPLING_MODE:
//		*params = pContext->sampler.attrib_mode[attrib];
		break;
	case GI_TEXTURE_DIMENSION:
//		*params = pContext->sampler.texture_dim[attrib];
		break;
	case GI_GL_SAMPLE_TEXTURE:
//		*params = pContext->sampler.attrib_texture[attrib];
		break;
	default:
		{
			GIboolean bBool = GI_FALSE;
			giGetBooleanv(pname, &bBool);
			*params = bBool;
		}
	}
}

/** Retrieve floating point attribute state.
 *  \param attrib attribute to query
 *  \param pname state to query
 *  \param params address of variable to store value
 *  \ingroup state
 */
void GIAPIENTRY giGetAttribfv(GIuint attrib, GIenum pname, GIfloat *params)
{
	GIContext *pContext = GIContext_current();
	if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}

	/* select state and get value */
	switch(pname)
	{
	case GI_SAMPLING_TRANSFORM:
		//memcpy(params, pContext->sampler
		//	.attrib_matrix[attrib], 16*sizeof(GIfloat));
		break;
	default:
		{
			GIint iInt = 0;
			giGetAttribiv(attrib, pname, &iInt);
			*params = iInt;
		}
	}
}

/** Get pointer to attribute array.
 *  \param attrib attribute to query
 *  \param params address to store pointer at
 *  \ingroup state
 */
void GIAPIENTRY giGetAttribPointerv(GIuint attrib, GIvoid **params)
{
	/* get pointer */
	if(attrib < GI_ATTRIB_COUNT)
		*params = GIContext_current()->attrib_pointer[attrib];
	else
		GIContext_error(GIContext_current(), GI_INVALID_VALUE);
}

/** Retrieve and clear last error.
 *  \return error code of last encountered error
 *  \ingroup error
 */
GIenum GIAPIENTRY giGetError()
{
	/* return and clear error code */
	GIuint uiError = GIContext_current()->error;
	GIContext_current()->error = GI_NO_ERROR;
	return uiError;
}

/** Get descriptive error string.
 *  \param errorCode error code to get description for
 *  \return descriptive string for given error code
 *  \ingroup error
 */
const GIchar* GIAPIENTRY giErrorString(GIenum errorCode)
{
	static const GIchar *szNoError = "正常終了";
#if 0
	static const GIchar *szInvalidEnum = "enumeration value out of range";
	static const GIchar *szInvalidOp = "operation illegal in current state";
	static const GIchar *szInvalidValue = "numeric value out of range";
	static const GIchar *szInvalidID = "object with specified id does not exist";
#else
	static const GIchar *szInvalidEnum = "Invalid";
	static const GIchar *szInvalidOp = "数値計算エラー";
	static const GIchar *szInvalidValue = "数値計算エラー";
	static const GIchar *szInvalidID = "Invalid";
#endif
	static const GIchar *szInvalidMesh = "非多様体";
	static const GIchar *szInvalidCut = "カット出来ない";
	static const GIchar *szInvalidParameterization = "パラメータ化失敗";
	static const GIchar *szNumericalError = "数値計算エラー";
	static const GIchar *szUnsupportedOperation = "内部エラー";

	/* select appropriate description */
	switch(errorCode)
	{
	case GI_INVALID_ENUM:
		return szInvalidEnum;
		break;
	case GI_INVALID_OPERATION:
		return szInvalidOp;
		break;
	case GI_INVALID_VALUE:
		return szInvalidValue;
		break;
	case GI_INVALID_ID:
		return szInvalidID;
		break;
	case GI_INVALID_MESH:
		return szInvalidMesh;
		break;
	case GI_INVALID_CUT:
		return szInvalidCut;
		break;
	case GI_INVALID_PARAMETERIZATION:
		return szInvalidParameterization;
		break;
	case GI_NUMERICAL_ERROR:
		return szNumericalError;
		break;
	case GI_UNSUPPORTED_OPERATION:
		return szUnsupportedOperation;
		break;
	default:
		return szNoError;
	}
}

/** Set error callback function.
 *  \param fn function to use
 *  \param data user data
 *  \ingroup error
 */
void GIAPIENTRY giErrorCallback(GIerrorcb fn, GIvoid *data)
{
	/* set callback */
	GIContext_current()->error_cb = fn;
	GIContext_current()->edata = data;
}

/** Get enum value by name.
 *  \param name enumeration name as string
 *  \return enumeration value
 *  \ingroup state
 */
GIenum GIAPIENTRY giGetEnumValue(const GIchar *name)
{
	static GIHash hEnumMap = { NULL, 0, 0, 0, 0, 0.0f, NULL, NULL, NULL };
#if 0
	if(!hEnumMap.size)
	{
		/* create and fill map */
		GIHash_construct(&hEnumMap, 256, 1.0f, 48*sizeof(GIchar*), 
			hash_string, compare_string, copy_string);
		GIHash_insert(&hEnumMap, "GI_TRUE", (GIvoid*)GI_TRUE);
		GIHash_insert(&hEnumMap, "GI_FALSE", (GIvoid*)GI_FALSE);
		GIHash_insert(&hEnumMap, "GI_NULL", (GIvoid*)GI_NULL);
		GIHash_insert(&hEnumMap, "GI_NONE", (GIvoid*)GI_NONE);
		GIHash_insert(&hEnumMap, "GI_ALL_PATCHES", (GIvoid*)GI_ALL_PATCHES);
		GIHash_insert(&hEnumMap, "GI_BYTE", (GIvoid*)GI_BYTE);
		GIHash_insert(&hEnumMap, "GI_UNSIGNED_BYTE", (GIvoid*)GI_UNSIGNED_BYTE);
		GIHash_insert(&hEnumMap, "GI_SHORT", (GIvoid*)GI_SHORT);
		GIHash_insert(&hEnumMap, "GI_UNSIGNED_SHORT", (GIvoid*)GI_UNSIGNED_SHORT);
		GIHash_insert(&hEnumMap, "GI_INT", (GIvoid*)GI_INT);
		GIHash_insert(&hEnumMap, "GI_UNSIGNED_INT", (GIvoid*)GI_UNSIGNED_INT);
		GIHash_insert(&hEnumMap, "GI_FLOAT", (GIvoid*)GI_FLOAT);
		GIHash_insert(&hEnumMap, "GI_DOUBLE", (GIvoid*)GI_DOUBLE);
		GIHash_insert(&hEnumMap, "GI_HALF_FLOAT", (GIvoid*)GI_HALF_FLOAT);
		GIHash_insert(&hEnumMap, "GI_MULTITHREADING", (GIvoid*)GI_MULTITHREADING);
		GIHash_insert(&hEnumMap, "GI_EXACT_MAPPING_SUBSET", (GIvoid*)GI_EXACT_MAPPING_SUBSET);
		GIHash_insert(&hEnumMap, "GI_PARAM_CORNER_SUBSET", (GIvoid*)GI_PARAM_CORNER_SUBSET);
		GIHash_insert(&hEnumMap, "GI_VERSION", (GIvoid*)GI_VERSION);
		GIHash_insert(&hEnumMap, "GI_MAX_ATTRIBS", (GIvoid*)GI_MAX_ATTRIBS);
		GIHash_insert(&hEnumMap, "GI_MESH_BINDING", (GIvoid*)GI_MESH_BINDING);
		GIHash_insert(&hEnumMap, "GI_IMAGE_BINDING", (GIvoid*)GI_IMAGE_BINDING);
		GIHash_insert(&hEnumMap, "GI_SAMPLED_ATTRIB_COUNT", (GIvoid*)GI_SAMPLED_ATTRIB_COUNT);
		GIHash_insert(&hEnumMap, "GI_SAMPLED_ATTRIBS", (GIvoid*)GI_SAMPLED_ATTRIBS);
		GIHash_insert(&hEnumMap, "GI_EXACT_MAPPING_SUBSET_COUNT", (GIvoid*)GI_EXACT_MAPPING_SUBSET_COUNT);
		GIHash_insert(&hEnumMap, "GI_PARAM_CORNER_SUBSET_COUNT", (GIvoid*)GI_PARAM_CORNER_SUBSET_COUNT);
		GIHash_insert(&hEnumMap, "GI_EXACT_MAPPING_SUBSET_SORTED", (GIvoid*)GI_EXACT_MAPPING_SUBSET_SORTED);
		GIHash_insert(&hEnumMap, "GI_PARAM_CORNER_SUBSET_SORTED", (GIvoid*)GI_PARAM_CORNER_SUBSET_SORTED);
		GIHash_insert(&hEnumMap, "GI_POSITION_ATTRIB", (GIvoid*)GI_POSITION_ATTRIB);
		GIHash_insert(&hEnumMap, "GI_PARAM_ATTRIB", (GIvoid*)GI_PARAM_ATTRIB);
		GIHash_insert(&hEnumMap, "GI_PARAM_STRETCH_ATTRIB", (GIvoid*)GI_PARAM_STRETCH_ATTRIB);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_ARRAY_ENABLED", (GIvoid*)GI_ATTRIB_ARRAY_ENABLED);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_ARRAY_SIZE", (GIvoid*)GI_ATTRIB_ARRAY_SIZE);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_ARRAY_STRIDE", (GIvoid*)GI_ATTRIB_ARRAY_STRIDE);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_ARRAY_NORMALIZED", (GIvoid*)GI_ATTRIB_ARRAY_NORMALIZED);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_ARRAY_SEMANTIC", (GIvoid*)GI_ATTRIB_ARRAY_SEMANTIC);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_IMAGE", (GIvoid*)GI_ATTRIB_IMAGE);
		GIHash_insert(&hEnumMap, "GI_HAS_CUT", (GIvoid*)GI_HAS_CUT);
		GIHash_insert(&hEnumMap, "GI_HAS_PARAMS", (GIvoid*)GI_HAS_PARAMS);
		GIHash_insert(&hEnumMap, "GI_PATCH_COUNT", (GIvoid*)GI_PATCH_COUNT);
		GIHash_insert(&hEnumMap, "GI_FACE_COUNT", (GIvoid*)GI_FACE_COUNT);
		GIHash_insert(&hEnumMap, "GI_EDGE_COUNT", (GIvoid*)GI_EDGE_COUNT);
		GIHash_insert(&hEnumMap, "GI_VERTEX_COUNT", (GIvoid*)GI_VERTEX_COUNT);
		GIHash_insert(&hEnumMap, "GI_PARAM_COUNT", (GIvoid*)GI_PARAM_COUNT);
		GIHash_insert(&hEnumMap, "GI_AABB_MIN", (GIvoid*)GI_AABB_MIN);
		GIHash_insert(&hEnumMap, "GI_AABB_MAX", (GIvoid*)GI_AABB_MAX);
		GIHash_insert(&hEnumMap, "GI_RADIUS", (GIvoid*)GI_RADIUS);
		GIHash_insert(&hEnumMap, "GI_ACTIVE_PATCH", (GIvoid*)GI_ACTIVE_PATCH);
		GIHash_insert(&hEnumMap, "GI_MIN_PARAM_STRETCH", (GIvoid*)GI_MIN_PARAM_STRETCH);
		GIHash_insert(&hEnumMap, "GI_MAX_PARAM_STRETCH", (GIvoid*)GI_MAX_PARAM_STRETCH);
		GIHash_insert(&hEnumMap, "GI_PARAM_STRETCH_METRIC", (GIvoid*)GI_PARAM_STRETCH_METRIC);
		GIHash_insert(&hEnumMap, "GI_TOPOLOGICAL_SIDEBAND_LENGTH", (GIvoid*)GI_TOPOLOGICAL_SIDEBAND_LENGTH);
		GIHash_insert(&hEnumMap, "GI_TOPOLOGICAL_SIDEBAND", (GIvoid*)GI_TOPOLOGICAL_SIDEBAND);
		GIHash_insert(&hEnumMap, "GI_HAS_ATTRIB", (GIvoid*)GI_HAS_ATTRIB);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_SIZE", (GIvoid*)GI_ATTRIB_SIZE);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_NORMALIZED", (GIvoid*)GI_ATTRIB_NORMALIZED);
		GIHash_insert(&hEnumMap, "GI_ATTRIB_SEMANTIC", (GIvoid*)GI_ATTRIB_SEMANTIC);
		GIHash_insert(&hEnumMap, "GI_CUTTER", (GIvoid*)GI_CUTTER);
		GIHash_insert(&hEnumMap, "GI_SUBDIVISION_ITERATIONS", (GIvoid*)GI_SUBDIVISION_ITERATIONS);
		GIHash_insert(&hEnumMap, "GI_INITIAL_GIM", (GIvoid*)GI_INITIAL_GIM);
		GIHash_insert(&hEnumMap, "GI_CATMULL_CLARK_SUBDIVISION", (GIvoid*)GI_CATMULL_CLARK_SUBDIVISION);
		GIHash_insert(&hEnumMap, "GI_PARAMETERIZER", (GIvoid*)GI_PARAMETERIZER);
		GIHash_insert(&hEnumMap, "GI_INITIAL_PARAMETERIZATION", (GIvoid*)GI_INITIAL_PARAMETERIZATION);
		GIHash_insert(&hEnumMap, "GI_STRETCH_METRIC", (GIvoid*)GI_STRETCH_METRIC);
		GIHash_insert(&hEnumMap, "GI_CONFORMAL_WEIGHT", (GIvoid*)GI_CONFORMAL_WEIGHT);
		GIHash_insert(&hEnumMap, "GI_AUTHALIC_WEIGHT", (GIvoid*)GI_AUTHALIC_WEIGHT);
		GIHash_insert(&hEnumMap, "GI_STRETCH_WEIGHT", (GIvoid*)GI_STRETCH_WEIGHT);
		GIHash_insert(&hEnumMap, "GI_PARAM_RESOLUTION", (GIvoid*)GI_PARAM_RESOLUTION);
		GIHash_insert(&hEnumMap, "GI_UNSYMMETRIC_SOLVER", (GIvoid*)GI_UNSYMMETRIC_SOLVER);
		GIHash_insert(&hEnumMap, "GI_AREA_WEIGHT", (GIvoid*)GI_AREA_WEIGHT);
		GIHash_insert(&hEnumMap, "GI_PARAM_SOURCE_ATTRIB", (GIvoid*)GI_PARAM_SOURCE_ATTRIB);
		GIHash_insert(&hEnumMap, "GI_FROM_ATTRIB", (GIvoid*)GI_FROM_ATTRIB);
		GIHash_insert(&hEnumMap, "GI_TUTTE_BARYCENTRIC", (GIvoid*)GI_TUTTE_BARYCENTRIC);
		GIHash_insert(&hEnumMap, "GI_SHAPE_PRESERVING", (GIvoid*)GI_SHAPE_PRESERVING);
		GIHash_insert(&hEnumMap, "GI_DISCRETE_HARMONIC", (GIvoid*)GI_DISCRETE_HARMONIC);
		GIHash_insert(&hEnumMap, "GI_MEAN_VALUE", (GIvoid*)GI_MEAN_VALUE);
		GIHash_insert(&hEnumMap, "GI_DISCRETE_AUTHALIC", (GIvoid*)GI_DISCRETE_AUTHALIC);
		GIHash_insert(&hEnumMap, "GI_INTRINSIC", (GIvoid*)GI_INTRINSIC);
		GIHash_insert(&hEnumMap, "GI_STRETCH_MINIMIZING", (GIvoid*)GI_STRETCH_MINIMIZING);
		GIHash_insert(&hEnumMap, "GI_GIM", (GIvoid*)GI_GIM);
		GIHash_insert(&hEnumMap, "GI_SOLVER_BICGSTAB", (GIvoid*)GI_SOLVER_BICGSTAB);
		GIHash_insert(&hEnumMap, "GI_SOLVER_GMRES", (GIvoid*)GI_SOLVER_GMRES);
		GIHash_insert(&hEnumMap, "GI_PARAM_STARTED", (GIvoid*)GI_PARAM_STARTED);
		GIHash_insert(&hEnumMap, "GI_PARAM_CHANGED", (GIvoid*)GI_PARAM_CHANGED);
		GIHash_insert(&hEnumMap, "GI_PARAM_FINISHED", (GIvoid*)GI_PARAM_FINISHED);
		GIHash_insert(&hEnumMap, "GI_MAX_GEOMETRIC_STRETCH", (GIvoid*)GI_MAX_GEOMETRIC_STRETCH);
		GIHash_insert(&hEnumMap, "GI_RMS_GEOMETRIC_STRETCH", (GIvoid*)GI_RMS_GEOMETRIC_STRETCH);
		GIHash_insert(&hEnumMap, "GI_COMBINED_STRETCH", (GIvoid*)GI_COMBINED_STRETCH);
		GIHash_insert(&hEnumMap, "GI_IMAGE_WIDTH", (GIvoid*)GI_IMAGE_WIDTH);
		GIHash_insert(&hEnumMap, "GI_IMAGE_HEIGHT", (GIvoid*)GI_IMAGE_HEIGHT);
		GIHash_insert(&hEnumMap, "GI_IMAGE_COMPONENTS", (GIvoid*)GI_IMAGE_COMPONENTS);
		GIHash_insert(&hEnumMap, "GI_IMAGE_TYPE", (GIvoid*)GI_IMAGE_TYPE);
		GIHash_insert(&hEnumMap, "GI_GL_IMAGE_TEXTURE", (GIvoid*)GI_GL_IMAGE_TEXTURE);
		GIHash_insert(&hEnumMap, "GI_GL_IMAGE_BUFFER", (GIvoid*)GI_GL_IMAGE_BUFFER);
		GIHash_insert(&hEnumMap, "GI_IMAGE_DATA", (GIvoid*)GI_IMAGE_DATA);
		GIHash_insert(&hEnumMap, "GI_IMAGE_STORAGE", (GIvoid*)GI_IMAGE_STORAGE);
		GIHash_insert(&hEnumMap, "GI_SUBIMAGE_X", (GIvoid*)GI_SUBIMAGE_X);
		GIHash_insert(&hEnumMap, "GI_SUBIMAGE_Y", (GIvoid*)GI_SUBIMAGE_Y);
		GIHash_insert(&hEnumMap, "GI_SUBIMAGE_WIDTH", (GIvoid*)GI_SUBIMAGE_WIDTH);
		GIHash_insert(&hEnumMap, "GI_SUBIMAGE_HEIGHT", (GIvoid*)GI_SUBIMAGE_HEIGHT);
		GIHash_insert(&hEnumMap, "GI_SUBIMAGE", (GIvoid*)GI_SUBIMAGE);
		GIHash_insert(&hEnumMap, "GI_EXTERNAL_DATA", (GIvoid*)GI_EXTERNAL_DATA);
		GIHash_insert(&hEnumMap, "GI_GL_TEXTURE_DATA", (GIvoid*)GI_GL_TEXTURE_DATA);
		GIHash_insert(&hEnumMap, "GI_GL_BUFFER_DATA", (GIvoid*)GI_GL_BUFFER_DATA);
		GIHash_insert(&hEnumMap, "GI_NO_IMAGE_DATA", (GIvoid*)GI_NO_IMAGE_DATA);
		GIHash_insert(&hEnumMap, "GI_SAMPLER", (GIvoid*)GI_SAMPLER);
		GIHash_insert(&hEnumMap, "GI_SAMPLER_USE_SHADER", (GIvoid*)GI_SAMPLER_USE_SHADER);
		GIHash_insert(&hEnumMap, "GI_SAMPLER_USE_RENDER_TO_TEXTURE", (GIvoid*)GI_SAMPLER_USE_RENDER_TO_TEXTURE);
		GIHash_insert(&hEnumMap, "GI_SAMPLER_SOFTWARE", (GIvoid*)GI_SAMPLER_SOFTWARE);
		GIHash_insert(&hEnumMap, "GI_SAMPLER_OPENGL", (GIvoid*)GI_SAMPLER_OPENGL);
		GIHash_insert(&hEnumMap, "GI_SAMPLING_MODE", (GIvoid*)GI_SAMPLING_MODE);
		GIHash_insert(&hEnumMap, "GI_SAMPLING_TRANSFORM", (GIvoid*)GI_SAMPLING_TRANSFORM);
		GIHash_insert(&hEnumMap, "GI_TEXTURE_DIMENSION", (GIvoid*)GI_TEXTURE_DIMENSION);
		GIHash_insert(&hEnumMap, "GI_GL_SAMPLE_TEXTURE", (GIvoid*)GI_GL_SAMPLE_TEXTURE);
		GIHash_insert(&hEnumMap, "GI_SAMPLE_DEFAULT", (GIvoid*)GI_SAMPLE_DEFAULT);
		GIHash_insert(&hEnumMap, "GI_SAMPLE_NORMALIZED", (GIvoid*)GI_SAMPLE_NORMALIZED);
		GIHash_insert(&hEnumMap, "GI_SAMPLE_TEXTURED", (GIvoid*)GI_SAMPLE_TEXTURED);
		GIHash_insert(&hEnumMap, "GI_GL_USE_VERTEX_TEXTURE", (GIvoid*)GI_GL_USE_VERTEX_TEXTURE);
		GIHash_insert(&hEnumMap, "GI_GL_USE_GEOMETRY_SHADER", (GIvoid*)GI_GL_USE_GEOMETRY_SHADER);
		GIHash_insert(&hEnumMap, "GI_RENDER_RESOLUTION_U", (GIvoid*)GI_RENDER_RESOLUTION_U);
		GIHash_insert(&hEnumMap, "GI_RENDER_RESOLUTION_V", (GIvoid*)GI_RENDER_RESOLUTION_V);
		GIHash_insert(&hEnumMap, "GI_RENDER_CACHE_SIZE", (GIvoid*)GI_RENDER_CACHE_SIZE);
		GIHash_insert(&hEnumMap, "GI_GL_RENDER_SEMANTIC", (GIvoid*)GI_GL_RENDER_SEMANTIC);
		GIHash_insert(&hEnumMap, "GI_GL_RENDER_CHANNEL", (GIvoid*)GI_GL_RENDER_CHANNEL);
		GIHash_insert(&hEnumMap, "GI_TEXTURE_COORD_DOMAIN", (GIvoid*)GI_TEXTURE_COORD_DOMAIN);
		GIHash_insert(&hEnumMap, "GI_GL_VERTEX", (GIvoid*)GI_GL_VERTEX);
		GIHash_insert(&hEnumMap, "GI_GL_NORMAL", (GIvoid*)GI_GL_NORMAL);
		GIHash_insert(&hEnumMap, "GI_GL_COLOR", (GIvoid*)GI_GL_COLOR);
		GIHash_insert(&hEnumMap, "GI_GL_SECONDARY_COLOR", (GIvoid*)GI_GL_SECONDARY_COLOR);
		GIHash_insert(&hEnumMap, "GI_GL_FOG_COORD", (GIvoid*)GI_GL_FOG_COORD);
		GIHash_insert(&hEnumMap, "GI_GL_EVAL_COORD", (GIvoid*)GI_GL_EVAL_COORD);
		GIHash_insert(&hEnumMap, "GI_GL_TEXTURE_COORD", (GIvoid*)GI_GL_TEXTURE_COORD);
		GIHash_insert(&hEnumMap, "GI_GL_VERTEX_ATTRIB", (GIvoid*)GI_GL_VERTEX_ATTRIB);
		GIHash_insert(&hEnumMap, "GI_UNIT_SQUARE", (GIvoid*)GI_UNIT_SQUARE);
		GIHash_insert(&hEnumMap, "GI_HALF_TEXEL_INDENT", (GIvoid*)GI_HALF_TEXEL_INDENT);
		GIHash_insert(&hEnumMap, "GI_NO_ERROR", (GIvoid*)GI_NO_ERROR);
		GIHash_insert(&hEnumMap, "GI_INVALID_ENUM", (GIvoid*)GI_INVALID_ENUM);
		GIHash_insert(&hEnumMap, "GI_INVALID_OPERATION", (GIvoid*)GI_INVALID_OPERATION);
		GIHash_insert(&hEnumMap, "GI_INVALID_VALUE", (GIvoid*)GI_INVALID_VALUE);
		GIHash_insert(&hEnumMap, "GI_INVALID_ID", (GIvoid*)GI_INVALID_ID);
		GIHash_insert(&hEnumMap, "GI_INVALID_MESH", (GIvoid*)GI_INVALID_MESH);
		GIHash_insert(&hEnumMap, "GI_INVALID_CUT", (GIvoid*)GI_INVALID_CUT);
		GIHash_insert(&hEnumMap, "GI_NUMERICAL_ERROR", (GIvoid*)GI_NUMERICAL_ERROR);
		GIHash_insert(&hEnumMap, "GI_UNSUPPORTED_OPERATION", (GIvoid*)GI_UNSUPPORTED_OPERATION);
		GIHash_insert(&hEnumMap, "GI_INVALID_PARAMETERIZATION", (GIvoid*)GI_INVALID_PARAMETERIZATION);
	}
#else
printf("aaaaaaaaaaaaaaaa\n");
#endif

	/* return enum value */
	return (GIuint)GIHash_find(&hEnumMap, name);
}

/** \internal
 *  \brief Get current context.
 *  \return pointer to current context
 *  \ingroup context
 */
GIContext* GIContext_current()
{
	/* return current context */
	return gi_CurrentContext;
}

/** \internal
 *  \brief Set last error and put out if in debug mode.
 *  \param context context to work on
 *  \param error error code
 *  \ingroup error
 */
void GIContext_error(GIContext *context, GIenum error)
{
	/* set and print error */
	context->error = error;
	if(context->error_cb)
		context->error_cb(error, context->edata);
	else
	{
		GIDebug(fprintf(stderr, "エラー[%s]\n", giErrorString(error)));
		fprintf(stderr, "エラー[%s]\n", giErrorString(error));
		printf("エラー[%s]\n", giErrorString(error));
	}
}
