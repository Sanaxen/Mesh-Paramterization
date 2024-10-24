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
 *  \brief Implementation of structures and functions for parameterizing triangular meshes.
 */

#include "gi_context.h"
#include "gi_parameterizer.h"
#include "gi_math.h"
#include "gi_memory.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>

#define GI_HALF_SQRT_3				0.8660254037844386


/** \internal
 *  \brief Param rescue structure.
 */
typedef struct _GIParamSave
{
	GIuint		id;									/**< Id of param. */
	GIdouble	params[2];							/**< Parameter coordinates. */
	GIdouble	stretch;							/**< Stretch value. */
} GIParamSave;

/** \internal
 *  \brief Information about path during A* or Dijkstra.
*/
typedef struct _GIPathInfo
{
	GIHalfEdge	*source;							/**< Source half edge. */
	GIdouble	distance;							/**< Distance to root. */
} GIPathInfo;

/** \internal
 *  \brief Information about local parameterization.
*/
typedef struct _GILocalInfo
{
	GIdouble	param[2];							/**< Parameter coordinates. */
	GIdouble	angle;								/**< Surface angle. */
} GILocalInfo;


/** \internal
 *  \brief Compare path pointers for qsort.
 *  \param a first value
 *  \param b second value
 *  \return negative value if a > b, positive value if b > a and 0 if equal
 */
static int compare(const void *a, const void *b)
{
	GICutPath** p = (GICutPath**)a;
	GICutPath** q = (GICutPath**)b;
	return ((*q)->glength - (*p)->glength);
}

/** Set boolean configuration parameter of parameterization.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup parameterization
 */
void GIAPIENTRY giParameterizerParameterb(GIenum pname, GIboolean param)
{
	GIParameterizer *pPar = &(GIContext_current()->parameterizer);

	/* select state and set value */
	switch(pname)
	{
	default:
		GIContext_error(pPar->context, GI_INVALID_ENUM);
	}
}

/** Set integer configuration parameter of parameterization.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup parameterization
 */
void GIAPIENTRY giParameterizerParameteri(GIenum pname, GIint param)
{
	GIParameterizer *pPar = &(GIContext_current()->parameterizer);

	/* select state and set value */
	switch(pname)
	{
	case GI_PARAMETERIZER:
		switch(param)
		{
		case GI_FROM_ATTRIB:
		case GI_TUTTE_BARYCENTRIC:
		case GI_SHAPE_PRESERVING:
		case GI_DISCRETE_HARMONIC:
		case GI_MEAN_VALUE:
		case GI_DISCRETE_AUTHALIC:
		case GI_INTRINSIC:
		case GI_STRETCH_MINIMIZING:
		case GI_GIM:
			pPar->parameterizer = param;
			break;
		default:
			GIContext_error(pPar->context, GI_INVALID_ENUM);
		}
		break;
	case GI_INITIAL_PARAMETERIZATION:
		switch(param)
		{
		case GI_SHAPE_PRESERVING:
		case GI_DISCRETE_HARMONIC:
		case GI_MEAN_VALUE:
		case GI_DISCRETE_AUTHALIC:
			pPar->initial_param = param;
			break;
		default:
			GIContext_error(pPar->context, GI_INVALID_ENUM);
		}
		break;
	case GI_STRETCH_METRIC:
		if(param < GI_STRETCH_BASE || param > GI_STRETCH_END)
			GIContext_error(pPar->context, GI_INVALID_ENUM);
		else
			pPar->stretch_metric = param;
		break;
	case GI_PARAM_RESOLUTION:
		if((GIuint)param < 2)
		{
			GIContext_error(pPar->context, GI_INVALID_VALUE);
			return;
		}
		pPar->sampling_res = param;
		break;
	case GI_UNSYMMETRIC_SOLVER:
		if(param == GI_SOLVER_BICGSTAB || param == GI_SOLVER_GMRES)
			pPar->solver = param;
		else
			GIContext_error(pPar->context, GI_INVALID_ENUM);
		break;
	case GI_PARAM_SOURCE_ATTRIB:
		if(param < GI_ATTRIB_COUNT)
			pPar->source_attrib = param;
		else
			GIContext_error(pPar->context, GI_INVALID_VALUE);
		break;
	default:
		GIContext_error(pPar->context, GI_INVALID_ENUM);
	}
}

/** Set floating point configuration parameter of parameterization.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup parameterization
 */
void GIAPIENTRY giParameterizerParameterf(GIenum pname, GIfloat param)
{
	GIParameterizer *pPar = &(GIContext_current()->parameterizer);

	/* select state and set value */
	switch(pname)
	{
	case GI_CONFORMAL_WEIGHT:
		if(param >= 0.0f)
			pPar->conformal_weight = param;
		else
			GIContext_error(pPar->context, GI_INVALID_VALUE);
		break;
	case GI_AUTHALIC_WEIGHT:
		if(param >= 0.0f)
			pPar->authalic_weight = param;
		else
			GIContext_error(pPar->context, GI_INVALID_VALUE);
		break;
	case GI_STRETCH_WEIGHT:
		if(param >= 0.0f && param <= 1.0f)
			pPar->stretch_weight = param;
		else
			GIContext_error(pPar->context, GI_INVALID_VALUE);
		break;
	case GI_AREA_WEIGHT:
		if(param >= 0.0f)
			pPar->area_weight = param;
		else
			GIContext_error(pPar->context, GI_INVALID_VALUE);
		break;
	default:
		GIContext_error(pPar->context, GI_INVALID_ENUM);
	}
}

/** Set callback function for parameterization.
 *  \param which callback to set
 *  \param fn function to use
 *  \param data custom user data
 *  \ingroup parameterization
 */
void GIAPIENTRY giParameterizerCallback(GIenum which, GIparamcb fn, GIvoid *data)
{
	GIParameterizer *pPar = &(GIContext_current()->parameterizer);

	/* error checking */
	if(which < GI_CALLBACK_BASE || which > GI_CALLBACK_END)
	{
		GIContext_error(pPar->context, GI_INVALID_ENUM);
		return;
	}

	/* set callback data */
	pPar->callback[which-GI_CALLBACK_BASE] = fn;
	pPar->cdata[which-GI_CALLBACK_BASE] = data;
}

/** Parameterize current mesh.
 *  This function computes parameter coordinates for the current bound mesh.
 *  \ingroup parameterization
 */
void GIAPIENTRY giParameterize()
{
	GIParameterizer *pPar = &(GIContext_current()->parameterizer);
	GIMesh *pMesh = pPar->context->mesh;
	GIPatch *pPatch, *pPStart, *pPEnd;
	GIboolean bSuccess, bActive = GI_FALSE;

	/* error checking */
	if(!pMesh || (!pMesh->patches && pPar->parameterizer!=GI_FROM_ATTRIB))
	{
		GIContext_error(pPar->context, GI_INVALID_OPERATION);
		return;
	}

	/* call back */
	if(pPar->callback[GI_PARAM_STARTED-GI_CALLBACK_BASE] && 
		!pPar->callback[GI_PARAM_STARTED-GI_CALLBACK_BASE](
		pPar->cdata[GI_PARAM_STARTED-GI_CALLBACK_BASE]))
		return;

	/* 1 or more patches? */
	if(pMesh->active_patch)
	{
		pPatch = pPStart = pMesh->active_patch;
		pPEnd = pPatch->next;
		bActive = GI_TRUE;
	}
	else
		pPatch = pPStart = pPEnd = pMesh->patches;

	/* texCoords as parameter coordinates */
	if(pPar->parameterizer == GI_FROM_ATTRIB)
	{
		GIFace *pFace;
		GIHalfEdge *pHalfEdge;
		GIfloat *pParams, *p;
		GIuint h, uiAttrib = pPar->source_attrib;
		GIint iOffset = pMesh->aoffset[uiAttrib];
		GIboolean bPos = pMesh->asemantic[uiAttrib] == GI_POSITION_ATTRIB;
		if(bActive || pMesh->asemantic[uiAttrib] == GI_PARAM_ATTRIB || 
			pMesh->asemantic[uiAttrib] == GI_PARAM_STRETCH_ATTRIB || 
			(!bPos && iOffset < 0) || pMesh->asize[uiAttrib] < 2)
		{
			GIContext_error(pPar->context, GI_INVALID_OPERATION);
			if(pPar->callback[GI_PARAM_FINISHED-GI_CALLBACK_BASE])
				pPar->callback[GI_PARAM_FINISHED-GI_CALLBACK_BASE](
					pPar->cdata[GI_PARAM_FINISHED-GI_CALLBACK_BASE]);
			return;
		}

		/* prepare for automatic patch generation and do it */
		GIMesh_destroy_cut(pMesh);
		if(bPos)
			pParams = (GIfloat*)GI_MALLOC_ARRAY(6*pMesh->fcount, sizeof(GIfloat));
		GI_LIST_FOREACH(pMesh->faces, pFace)
			h = 0;
			GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
				if(bPos)
				{
					p = pParams + (6*pFace->id+((h++)<<1));
					GI_VEC2_COPY(p, pHalfEdge->vstart->coords);
					pHalfEdge->pstart = (GIParam*)p;
				}
				else
					pHalfEdge->pstart = (GIParam*)((GIbyte*)pHalfEdge->astart+iOffset);
			GI_LIST_NEXT(pFace->hedges, pHalfEdge)
		GI_LIST_NEXT(pMesh->faces, pFace)
		if(GICutter_from_params(&pPar->context->cutter, pMesh))
		{
			pMesh->resolution = UINT_MAX;
			if(pPar->callback[GI_PARAM_CHANGED-GI_CALLBACK_BASE] && 
			   !pPar->callback[GI_PARAM_CHANGED-GI_CALLBACK_BASE](
			   pPar->cdata[GI_PARAM_CHANGED-GI_CALLBACK_BASE]))
				GIMesh_destroy_cut(pMesh);
		}
		if(bPos)
			GI_FREE_ARRAY(pParams);
		if(pPar->callback[GI_PARAM_FINISHED-GI_CALLBACK_BASE])
			pPar->callback[GI_PARAM_FINISHED-GI_CALLBACK_BASE](
				pPar->cdata[GI_PARAM_FINISHED-GI_CALLBACK_BASE]);
		return;
	}


	/* reverse boundary splits if reparameterizing all patches */
	if(!bActive && pMesh->resolution != pPar->sampling_res)
		GIMesh_revert_splits(pMesh, pMesh->split_hedges.size-pMesh->cut_splits);

	/* process patches */
	memset(pMesh->stretch, 0, GI_STRETCH_COUNT*sizeof(GIdouble));
	pMesh->param_metric = 0;
	do
	{
		GIDebug(printf("parameterizing patch %d\n", pPatch->id));
		pMesh->active_patch = pPatch;
		memset(pPatch->stretch, 0, GI_STRETCH_COUNT*sizeof(GIdouble));
		pPatch->param_metric = 0;

		/* former parameterization external -> search corners */
		if(pPatch->resolution == UINT_MAX)
		{
			memset(pPatch->corners, 0, 4*sizeof(GIParam*));
			GIPatch_find_corners(pPatch);
		}
		if(!bActive && pMesh->resolution != pPar->sampling_res)
			pPatch->resolution = 0;

		/* parameterize boundary */
		if(!GIParameterizer_arc_length_square(pPar, pPatch))
		{
			/* will be decremented again soon */
			if(!pPatch->parameterized)
				++pMesh->param_patches;
			bSuccess = GI_FALSE;
		}
		else if(pPatch->pcount > pPatch->hcount)
		{
			/* parameterize interior */
			switch(pPar->parameterizer)
			{
				case GI_TUTTE_BARYCENTRIC:
				case GI_SHAPE_PRESERVING:
				case GI_DISCRETE_HARMONIC:
				case GI_MEAN_VALUE:
				case GI_DISCRETE_AUTHALIC:
				case GI_INTRINSIC:
					{
						GILinearSystem system;
						GILinearSystem_construct(&system, pPar, pPatch, 
							pPar->parameterizer, GI_FALSE, GI_FALSE);
						bSuccess = GILinearSystem_solve(&system);
						GILinearSystem_unknowns_to_params(&system);
						GILinearSystem_destruct(&system);
					}
					break;
				case GI_STRETCH_MINIMIZING:
					bSuccess = GIParameterizer_stretch_minimizing(pPar, pPatch);
					break;
				case GI_GIM:
					if(pMesh->patch_count > 1)
					{
						GIContext_error(pPar->context, GI_INVALID_OPERATION);
						bSuccess = GI_FALSE;
					}
					else
						bSuccess = GIParameterizer_gim(pPar, pPatch);
			}
		}
		else
		{
			if(!pPatch->parameterized)
			{
				pPatch->parameterized = GI_TRUE;
				++pMesh->param_patches;
			}
			bSuccess = GI_TRUE;
		}

		/* parameterization successful or cancelled? */
		if(!bSuccess || (pPar->callback[GI_PARAM_CHANGED-GI_CALLBACK_BASE] && 
			!pPar->callback[GI_PARAM_CHANGED-GI_CALLBACK_BASE](
			pPar->cdata[GI_PARAM_CHANGED-GI_CALLBACK_BASE])))
		{
			--pMesh->param_patches;
			pPatch->parameterized = GI_FALSE;
			memset(pPatch->stretch, 0, GI_STRETCH_COUNT*sizeof(GIdouble));
			pPatch->param_metric = 0;
		}
		pPatch = pPatch->next;
		printf("--next.\n");
	}while(pPatch != pPEnd);

	/* finish parameterization */
	if(bActive)
		pMesh->resolution = 0;
	else
	{
		pMesh->active_patch = NULL;
		pMesh->resolution = pPar->sampling_res;
		if(pMesh->param_patches == pMesh->patch_count)
			pMesh->param_metric = pMesh->patches->param_metric;
	}
	if(pPar->callback[GI_PARAM_FINISHED-GI_CALLBACK_BASE])
		pPar->callback[GI_PARAM_FINISHED-GI_CALLBACK_BASE](
			pPar->cdata[GI_PARAM_FINISHED-GI_CALLBACK_BASE]);
	if(0)
	{
		GIParam *pParam;
		GI_LIST_FOREACH(pPatch->params, pParam)
			if(!pParam->cut_hedge)
			{
				GI_VEC2_SET(pParam->params, 0.5, 0.5);
			}
		GI_LIST_NEXT(pPatch->params, pParam)
	}
}

/** \internal
 *  \brief Parameterizer constructor.
 *  \param par parameterizer to construct
 *  \param context context to construct in
 *  \ingroup parameterization
 */
void GIParameterizer_construct(GIParameterizer *par, GIContext *context)
{
	/* initialize state */
	par->context = context;
	par->parameterizer = GI_STRETCH_MINIMIZING;
	par->initial_param = GI_MEAN_VALUE;
	par->stretch_metric = GI_RMS_GEOMETRIC_STRETCH;
	par->conformal_weight = 1.0f;
	par->authalic_weight = 1.0f;
	par->stretch_weight = 1.0f;
	par->area_weight = 1.0f;
	par->source_attrib = 0;
	par->sampling_res = 33;
	par->solver = GI_SOLVER_BICGSTAB;
	memset(par->callback, 0, GI_CALLBACK_COUNT*sizeof(GIparamcb));
	memset(par->cdata, 0, GI_CALLBACK_COUNT*sizeof(GIvoid*));
}

/** \internal
 *  \brief Parameterize border on unit circle with proportionally sampled points.
 *  \param par parameterizer to use
 *  \param patch patch to parameterize
 *  \retval GI_TRUE if border parameterized successfully
 *  \retval GI_FALSE if error occurred
 *  \ingroup parameterization
 */
GIboolean GIParameterizer_arc_length_circle(GIParameterizer *par, GIPatch *patch)
{
	GIParam *pParam = patch->params;
	GIdouble dAngle = 0.0, dTo2Pi = GI_TWO_PI / patch->hlength;
	if(!pParam)
		return GI_FALSE;

	/* parameterize cut */
	do
	{
		GI_VEC2_SET(pParam->params, cos(dAngle), sin(dAngle));
		dAngle += pParam->cut_hedge->edge->length * dTo2Pi;
		pParam = pParam->cut_hedge->next->pstart;
	}while(pParam != patch->params);
	patch->param_area = GI_PI;
	patch->resolution = UINT_MAX;
	return GI_TRUE;
}

/** \internal
 *  \brief Parameterize border on unit square with proportionally sampled points.
 *  \param par parameterizer to use
 *  \param patch patch to parameterize
 *  \retval GI_TRUE if border parameterized successfully
 *  \retval GI_FALSE if error occurred
 *  \ingroup parameterization
 */
GIboolean GIParameterizer_arc_length_square(GIParameterizer *par, GIPatch *patch)
{
	GIHalfEdge *pHalfEdge, *pCutHedge;
	GIParam *pParam = patch->params, *pPEnd;
	GICutPath *pPath = patch->paths, *pPath2;
	GICutPath **pPathArray, **pFreePaths;
	GIdouble *pLengthArray = NULL;
	GIint iSide = 0, iSidePos, iRes = par->sampling_res, iSideLen = iRes - 1, iLength;
	GIdouble dSidePos, dSideLen = (GIdouble)iSideLen, dInvRes = 1.0 / (GIdouble)iRes;
	GIdouble dLength, dRelLength, dLenToGLen;
	GIint iCorners[8] = { 0, 0, iSideLen, 0, iSideLen, iSideLen, 0, iSideLen };
	GIint iDirs[8] = { 1, 0, 0, 1, -1, 0, 0, -1 };
	GIint iVec[2];
	GIdouble dCorners[8] = { 0.0, 0.0, dSideLen, 0.0, dSideLen, dSideLen, 0.0, dSideLen };
	GIdouble dDirs[8] = { 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0 };
	GIdouble dEPS = 1e-9 * dSideLen;
	GIuint i, uiSides, uiGroups, uiOldStack = patch->split_paths.size;

	/* allready parameterized */
	if(patch->parameterized && patch->resolution == iRes)
		return GI_TRUE;

	/* corners given? */
	pPathArray = (GICutPath**)GI_CALLOC_ARRAY(patch->groups, sizeof(GICutPath*));
	if(patch->fixed_corners)
	{
		pPEnd = patch->corners[1];
		iLength = iSideLen;
		dLength = patch->side_lengths[0];
		uiSides = 4;
		pFreePaths = (GICutPath**)GI_MALLOC_ARRAY(patch->groups, sizeof(GICutPath*));
	}
	else
	{
		memset(patch->corners, 0, 4*sizeof(GIParam*));
		pPEnd = patch->params;
		iLength = iSideLen << 2;
		dLength = patch->hlength;
		uiSides = 1;
		pFreePaths = pPathArray;
		pLengthArray = (GIdouble*)GI_MALLOC_ARRAY(patch->groups, sizeof(GIdouble));
	}
	dLenToGLen = (GIdouble)iLength / dLength;

	/* process sides (or whole boundary) */
	for(i=0; i<uiSides; ++i)
	{
		/* allocate length */
		uiGroups = patch->fixed_corners ? 0 : patch->groups;
		do
		{
			if(pPathArray[pPath->group])
			{
				/* mate allready allocated */
				pPath->glength = pPath->twin->glength;
				iLength -= pPath->glength;
				dLength -= pPath->elength;
			}
			else
			{
				/* allocate length */
				dRelLength = pPath->elength * (iLength>0 ? ((GIdouble)iLength/dLength) : dLenToGLen);
				pPath->glength = (GIint)GI_ROUND(dRelLength);
				if(pPath->glength < 1)
					pPath->glength = 1;
				iLength -= pPath->glength;
				dLength -= pPath->elength;
				pPathArray[pPath->group] = pPath;
				if(patch->fixed_corners)
					pFreePaths[uiGroups++] = pPath;
				else
					pLengthArray[pPath->group] = pPath->elength;
			}
			pPath = pPath->next;
		}while(pPath->pstart != pPEnd);

		/* cut over- or underallocated? */
		if(iLength)
		{
			int loop = 0;
			GICutPath *pLastPath = NULL;
			GIint i = 0, iInc = (iLength > 0 ? 1 : -1);
			if(!uiGroups)
			{
				/* no free cut paths? */
				GI_FREE_ARRAY(pPathArray);
				GI_FREE_ARRAY(pFreePaths);
				return GI_FALSE;
			}
			qsort(pFreePaths, uiGroups, sizeof(GICutPath*), compare);
			while(iLength)
			{
				loop++;
				pPath2 = pFreePaths[i];
				if((pPath2->glength + iInc) != 0 && (!pPath2->twin || 
					pPath2->twin->patch != patch || abs(iLength) > 1))
				{
					/* shorten/lengthen cut path */
					pPath2->glength += iInc;
					iLength -= iInc;
					if(pPath2->twin && pPath2->twin->patch == patch)
					{
						pPath2->twin->glength += iInc;
						iLength -= iInc;
					}
					pLastPath = pPath2;
				}
				else if(pPath2 == pLastPath )
				{
					/* too much cut paths? */
					GI_FREE_ARRAY(pPathArray);
					if(patch->fixed_corners)
						GI_FREE_ARRAY(pFreePaths);
					else
						GI_FREE_ARRAY(pLengthArray);
					return GI_FALSE;
				}
				i = (i+1) % uiGroups;
				
				if ( loop > 10000 )
				{
					printf("位相不正=>無限ループ\n");
					break;
				}
			}
		}

		/* next side */
		if(patch->fixed_corners)
		{
			iLength = iSideLen;
			dLength = patch->side_lengths[i+1];
			pPEnd = patch->corners[(i+2)%4];
		}
	}
	GI_FREE_ARRAY(pPathArray);

	/* search fo(u)r corners if neccessary */
	if(!patch->fixed_corners)
	{
		GIint iNextCorner = iSideLen;
		iLength = 0;
		patch->corners[0] = patch->params;
		GI_LIST_FOREACH(patch->paths, pPath)
			if(iLength+pPath->glength > iNextCorner)
			{
				/* find param near corner or split half edge */
				dLength = 0.0;
				dLenToGLen = (GIdouble)pPath->glength / pPath->elength;
				pParam = pPath->pstart;
				do
				{
					pHalfEdge = pParam->cut_hedge;
					pParam = pHalfEdge->next->pstart;
					dLength += pHalfEdge->edge->length;
					dRelLength = (GIdouble)iLength + dLength*dLenToGLen;
					dSidePos = dRelLength - (GIdouble)iNextCorner;
					if(dSidePos > dEPS)
					{
						GIHalfEdge_split(pHalfEdge, patch, 
							pPath->twin ? pPath->twin->patch : NULL,  
							1.0-dSidePos/(pHalfEdge->edge->length*dLenToGLen), NULL);
						pParam = pHalfEdge->next->pstart;
						dLength -= pParam->cut_hedge->edge->length;
					}
				}while(dSidePos < -dEPS);

				/* split path */
				if(GIPatch_split_path(patch, pPath, pParam))
				{
					pPath2 = pPath->next;
					pPath2->glength = pPath->glength + iLength - iNextCorner;
					pPath2->elength = pPath->elength - dLength;
					pPath->elength = dLength;
					pPath->glength = iNextCorner - iLength;
					if(pPath->twin && pPath->twin->patch == patch)
					{
						pPath->twin->elength = pPath->elength;
						pPath->twin->glength = pPath->glength;
						pPath2->twin->elength = pPath2->elength;
						pPath2->twin->glength = pPath2->glength;
					}
				}
				patch->corners[iNextCorner/iSideLen] = pParam;
				iNextCorner += iSideLen;
			}
			else if(iLength+pPath->glength == iNextCorner)
			{
				iSide = iNextCorner / iSideLen;
				if(iSide < 4)
					patch->corners[iSide] = pPath->next->pstart;
				iNextCorner += iSideLen;
			}
			iLength += pPath->glength;
		GI_LIST_NEXT(patch->paths, pPath)
	}

	/* walk along cut paths */
	iLength = iSide = 0;
	GI_LIST_FOREACH(patch->paths, pPath)
		/* parameterize cut node */
		pParam = pPath->pstart;
		iSide = (iLength/iSideLen) << 1;
		dSidePos = iSidePos = iLength % iSideLen;
		GI_VEC2_ADD_SCALED(iVec, iCorners+iSide, iDirs+iSide, iSidePos);
		GI_VEC2_SCALE(pParam->params, iVec, 1.0/dSideLen);

		/* parameterize non-cut nodes */
		dLength = 0.0;
		dLenToGLen = (GIdouble)pPath->glength / pPath->elength;
		pHalfEdge = pParam->cut_hedge;
		pParam = pHalfEdge->next->pstart;
		while(pParam != pPath->next->pstart)
		{
			/* compute length and parameterize */
			pCutHedge = pParam->cut_hedge;
			dLength += pHalfEdge->edge->length;
			dRelLength = (GIdouble)iLength + dLength*dLenToGLen;
			iSide = (GIint)(dRelLength/dSideLen) << 1;
			dSidePos = fmod(dRelLength, dSideLen);
			GI_VEC2_ADD_SCALED(pParam->params, dCorners+iSide, dDirs+iSide, dSidePos);
			GI_VEC2_SCALE(pParam->params, pParam->params, 1.0/dSideLen);
			pHalfEdge = pCutHedge;
			pParam = pCutHedge->next->pstart;
		}
		iLength += pPath->glength;
	GI_LIST_NEXT(patch->paths, pPath)

	/* reverse temporary splits and reset path lengths (stability!) */
	if(!patch->fixed_corners)
	{
		GIPatch_revert_splits(patch, patch->split_paths.size-uiOldStack);
		GI_LIST_FOREACH(patch->paths, pPath)
			pPath->elength = pLengthArray[pPath->group];
		GI_LIST_NEXT(patch->paths, pPath)
		GI_FREE_ARRAY(pLengthArray);
	}
	else
		GI_FREE_ARRAY(pFreePaths);
	patch->resolution = iRes;
	patch->param_area = 1.0;
	return GI_TRUE;
}

/** \internal
 *  \brief Yoshizawa's iterative stretch minimizing parameterization.
 *  \param par parameterizer to use
 *  \param patch patch to parameterize
 *  \retval GI_TRUE if parameterized successfully
 *  \retval GI_FALSE on error or abort
 *  \ingroup parameterization
 */
GIboolean GIParameterizer_stretch_minimizing(GIParameterizer *par, 
											 GIPatch *patch)
{
	GILinearSystem system;
	GISparseMatrixCSR *A;
	GISparseMatrixLIL *B;
	GIVectorElement *pElement;
	GIParam **pBorderParams;
	GIdouble dSum;
	GIuint i, j, ii, ij, N, uiMetric = par->stretch_metric;
	GIParam *pParam;
	GIdouble *pOldParams;
	GIdouble *pPowStretches, *pOldStretches;
	GIdouble dOldStretch, dOldMin, dOldMax, dEta = par->stretch_weight;
	GIuint uiPCount = patch->pcount - patch->hcount;
	GIboolean bSuccess = GI_TRUE, bEta = (fabs(dEta-1.0) > 1e-4);

	/* create temporary arrays */
	pBorderParams = (GIParam**)GI_MALLOC_ARRAY(patch->hcount, sizeof(GIParam*));
	pOldParams = (GIdouble*)GI_MALLOC_ARRAY(uiPCount, 2*sizeof(GIdouble));
	pOldStretches = (GIdouble*)GI_MALLOC_ARRAY(patch->pcount, sizeof(GIdouble));
	if(bEta)
		pPowStretches = (GIdouble*)GI_MALLOC_ARRAY(patch->pcount, sizeof(GIdouble));
	else
		pPowStretches = pOldStretches;

	/* compute initial parameterization and stretch */
	GI_LIST_FOREACH(patch->params, pParam)
		pParam->stretch = 0.0;
		if(pParam->cut_hedge)
			pBorderParams[pParam->id-uiPCount] = pParam;
	GI_LIST_NEXT(patch->params, pParam)
	GILinearSystem_construct(&system, par, patch, par->initial_param, GI_TRUE, GI_TRUE);
	bSuccess = GILinearSystem_solve(&system);
	GILinearSystem_unknowns_to_params(&system);
	GIPatch_compute_stretch(patch, uiMetric, GI_TRUE, GI_FALSE);
	A = system.A;
	B = system.B;
	N = A->n;

	do
	{
		/* notify of changes */
		if (!bSuccess || (par->callback[GI_PARAM_CHANGED - GI_CALLBACK_BASE] &&
			!par->callback[GI_PARAM_CHANGED - GI_CALLBACK_BASE](
				par->cdata[GI_PARAM_CHANGED - GI_CALLBACK_BASE])))
		{
			/* clean up */
			GILinearSystem_destruct(&system);
			GI_FREE_ARRAY(pBorderParams);
			GI_FREE_ARRAY(pOldParams);
			GI_FREE_ARRAY(pOldStretches);
			if (bEta)
				GI_FREE_ARRAY(pPowStretches);

			//printf("notify of changes\n");
			return GI_FALSE;
		}

		/* save old parameterization and stretch */
		GI_LIST_FOREACH(patch->params, pParam)
			if (!pParam->cut_hedge)
			{
				GI_VEC2_COPY(pOldParams + (pParam->id << 1), pParam->params);
			}
		pOldStretches[pParam->id] = pParam->stretch;
		if (bEta)
			pPowStretches[pParam->id] = pow(pParam->stretch, dEta);
		pParam->stretch = 0.0;
		GI_LIST_NEXT(patch->params, pParam)
			dOldStretch = patch->stretch[uiMetric - GI_STRETCH_BASE];
		dOldMin = patch->min_param_stretch;
		dOldMax = patch->max_param_stretch;

		/* adjust weights and reparameterize */
		for(i=0,ij=0; i<N; ++i)
		{
			system.bU[i] = system.bV[i] = dSum = 0.0;
			for(; A->idx[ij]<i; ++ij)
			{
				A->values[ij] /= pPowStretches[A->idx[ij]];
				dSum -= A->values[ij];
			}
			ii = ij;
			for(++ij; ij<A->ptr[i+1]; ++ij)
			{
				A->values[ij] /= pPowStretches[A->idx[ij]];
				dSum -= A->values[ij];
			}
			for(j=0,pElement=B->rows[i].elements; j<B->rows[i].size; ++j,++pElement)
			{
				pElement->value /= pPowStretches[pElement->index];
				dSum += pElement->value;
				pParam = pBorderParams[pElement->index-uiPCount];
				system.bU[i] += pElement->value * pParam->params[0];
				system.bV[i] += pElement->value * pParam->params[1];
			}
			A->values[ii] = dSum;
		}
		bSuccess = GILinearSystem_solve(&system);
		GILinearSystem_unknowns_to_params(&system);
		if(bSuccess)
			GIPatch_compute_stretch(patch, uiMetric, GI_TRUE, GI_FALSE);
		printf("stretch %f => %f\n", (double)dOldStretch, (double)patch->stretch[uiMetric-GI_STRETCH_BASE]);

	}while(patch->stretch[uiMetric-GI_STRETCH_BASE] < dOldStretch );

	printf("\n");

	/* restore last parameterization and stretch */
	GI_LIST_FOREACH(patch->params, pParam)
		if(!pParam->cut_hedge)
		{
			GI_VEC2_COPY(pParam->params, pOldParams+(pParam->id<<1));
		}
		pParam->stretch = pOldStretches[pParam->id];
	GI_LIST_NEXT(patch->params, pParam)
	patch->stretch[uiMetric-GI_STRETCH_BASE] = dOldStretch;
	patch->min_param_stretch = dOldMin;
	patch->max_param_stretch = dOldMax;

	/* clean up */
	GILinearSystem_destruct(&system);
	GI_FREE_ARRAY(pBorderParams);
	GI_FREE_ARRAY(pOldParams);
	GI_FREE_ARRAY(pOldStretches);
	if(bEta)
		GI_FREE_ARRAY(pPowStretches);
	return GI_TRUE;
}

/** \internal
 *  \brief Dong's iterative stretch minimizing parameterization.
 *  \param par parameterizer to use
 *  \param patch patch to parameterize
 *  \retval GI_TRUE if parameterized successfully
 *  \retval GI_FALSE on error or abort
 *  \ingroup parameterization
 */
GIboolean GIParameterizer_stretch_minimizing2(GIParameterizer *par, 
											  GIPatch *patch)
{
	GILinearSystem system;
	GISparseVector row;
	GIFace *pFace = patch->faces, *pFEnd = patch->next->faces;
	GIHalfEdge *pHalfEdge, *pHEnd, *pPrev, *pNext;
	GIuint i, j, ij, uiMetric = par->stretch_metric;
	GIParam *pParam;
	GIHash hFaceAreas;
	GIdouble *pOldParams, *pArea;
	GIdouble *pOldStretches;
	GIdouble dOldStretch, dOldMin, dOldMax;
	GIuint uiPCount = patch->pcount - patch->hcount;
	GIboolean bSuccess = GI_TRUE;
	GIdouble dLength, dDiameter, dEPS;
	GIdouble dSimplex[6], dMid[2], dVec[2], dVNew[2], *dOld;
	GIdouble dValues[3], dNew, dExp;
	GIuint i1, i2, i3, j1, j2, j3, uiEval, uiEvalSum = 0;
	GIdouble v0[2], v1[2], v2[2];
	GIdouble dCos, dTan1, dTan2, dInvR0, dInvR1, dInvR2, dCoord, dSum;
	GIvoid *pArgs = &par->area_weight;

	/* create temporary datastructures */
	GISparseVector_construct(&row);
	pOldParams = (GIdouble*)GI_MALLOC_ARRAY(uiPCount, 2*sizeof(GIdouble));
	pOldStretches = (GIdouble*)GI_MALLOC_ARRAY(patch->pcount, sizeof(GIdouble));
	GIHash_construct(&hFaceAreas, patch->fcount, 0.0f, 
		sizeof(GIuint), hash_uint, compare_uint, copy_uint);
	do
	{
		pArea = (GIdouble*)GI_MALLOC_SINGLE(sizeof(GIdouble));
		*pArea = GIFace_area(pFace);
		GIHash_insert(&hFaceAreas, pFace, pArea);
		pFace = pFace->next;
	}while(pFace != pFEnd);

	/* compute initial parameterization and stretch */
	GILinearSystem_construct(&system, par, patch, par->initial_param, GI_TRUE, GI_FALSE);
	bSuccess = GILinearSystem_solve(&system);
	GILinearSystem_unknowns_to_params(&system);
	GIPatch_compute_stretch(patch, uiMetric, GI_TRUE, GI_TRUE);

	int count = 0;
	do
	{
		/* notify of changes */
		if(!bSuccess || (par->callback[GI_PARAM_CHANGED-GI_CALLBACK_BASE] && 
			!par->callback[GI_PARAM_CHANGED-GI_CALLBACK_BASE](
			par->cdata[GI_PARAM_CHANGED-GI_CALLBACK_BASE])))
		{
			/* clean up */
			GILinearSystem_destruct(&system);
			GI_FREE_ARRAY(pOldParams);
			GI_FREE_ARRAY(pOldStretches);
			GIHash_destruct(&hFaceAreas, sizeof(GIdouble));
			return GI_FALSE;
		}

		/* save old parameterization and compute new matrix */
		ij = 0;
		GI_LIST_FOREACH(patch->params, pParam)
			if(!pParam->cut_hedge)
			{
				dOld = pOldParams + (pParam->id<<1);
				GI_VEC2_COPY(dOld, pParam->params);

				/* compute diameter */
				dDiameter = 0.0f;
				pHalfEdge = pHEnd = pParam->vertex->hedge->twin;
				do
				{
					dLength = GIvec2d_dist_sqr(pParam->params, 
						pHalfEdge->pstart->params);
					if(dLength > dDiameter)
						dDiameter = dLength;
					pHalfEdge = pHalfEdge->next->twin;
				}while(pHalfEdge != pHEnd);
				dDiameter = 2.0 * sqrt(dDiameter);
				dEPS = 0.0025 * dDiameter;

				/* initial simplex (equilateral, edge length = 0.01*diameter) */
				dLength = 0.005 * dDiameter / GI_HALF_SQRT_3;
				pParam->params[1] += dLength;
				GI_VEC2_COPY(dSimplex, pParam->params);
				dValues[0] = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);
				GI_VEC2_SET(dVec, GI_HALF_SQRT_3, -0.5);
				GI_VEC2_ADD_SCALED(pParam->params, dOld, dVec, dLength);
				GI_VEC2_COPY(dSimplex+2, pParam->params);
				dValues[1] = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);
				dVec[0] = -GI_HALF_SQRT_3;
				GI_VEC2_ADD_SCALED(pParam->params, dOld, dVec, dLength);
				GI_VEC2_COPY(dSimplex+4, pParam->params);
				dValues[2] = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);

				/* Nelder-Mead-Optimization */
				uiEval = 3;
				do
				{
					/* sort simplex corners */
					if(dValues[0] < dValues[1])
					{
						if(dValues[1] < dValues[2]) { i1 = 0; i2 = 1; i3 = 2; }
						else if(dValues[0] < dValues[2]) { i1 = 0; i2 = 2; i3 = 1; }
						else { i1 = 2; i2 = 0; i3 = 1; }
					}
					else
					{
						if(dValues[0] < dValues[2]) { i1 = 1; i2 = 0; i3 = 2; }
						else if(dValues[1] < dValues[2]) { i1 = 1; i2 = 2; i3 = 0; }
						else { i1 = 2; i2 = 1; i3 = 0; }
					}
					j1 = i1 << 1; j2 = i2 << 1; j3 = i3 << 1;
					GI_VEC2_ADD(dMid, dSimplex+j1, dSimplex+j2);
					GI_VEC2_SCALE(dMid, dMid, 0.5);
					GI_VEC2_SUB(dVec, dMid, dSimplex+j3);

					/* reflection */
					GI_VEC2_ADD(pParam->params, dMid, dVec);
					dNew = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);
					if(dNew < dValues[i2] && dNew >= dValues[i1])
					{
						dLength = GIvec2d_dist(dSimplex+j3, pParam->params);
						GI_VEC2_COPY(dSimplex+j3, pParam->params);
						dValues[i3] = dNew;
					}
					else if(dNew < dValues[i1])
					{
						/* expansion */
						GI_VEC2_COPY(dVNew, pParam->params);
						GI_VEC2_ADD_SCALED(pParam->params, dMid, dVec, 2.0);
						dExp = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);
						if(dExp < dNew)
						{
							dLength = GIvec2d_dist(dSimplex+j3, pParam->params);
							GI_VEC2_COPY(dSimplex+j3, pParam->params);
							dValues[i3] = dExp;
						}
						else
						{
							dLength = GIvec2d_dist(dSimplex+j3, dVNew);
							GI_VEC2_COPY(dSimplex+j3, dVNew);
							dValues[i3] = dNew;
						}
						++uiEval;
					}
					else
					{
						/* contraction */
						GI_VEC2_ADD_SCALED(pParam->params, dSimplex+j3, dVec, 0.5);
						dNew = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);
						if(dNew < dValues[i3])
						{
							dLength = GIvec2d_dist(dSimplex+j3, pParam->params);
							GI_VEC2_COPY(dSimplex+j3, pParam->params);
							dValues[i3] = dNew;
						}
						else
						{
							/* reduction */
							GI_VEC2_COPY(pParam->params, dMid);
							dValues[i2] = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);
							GI_VEC2_ADD(pParam->params, dSimplex+j1, dSimplex+j3);
							GI_VEC2_SCALE(pParam->params, pParam->params, 0.5);
							dValues[i3] = GIParam_stretch(pParam, uiMetric, pArgs, &hFaceAreas);
							dLength = GIvec2d_dist(dSimplex+j2, dMid) + 
								GIvec2d_dist(dSimplex+j3, pParam->params);
							GI_VEC2_COPY(dSimplex+j2, dMid);
							GI_VEC2_COPY(dSimplex+j3, pParam->params);
							uiEval += 2;
						}
						++uiEval;
					}
					++uiEval;
				}while(uiEval < 100 && dLength > dEPS);

				/* take minimum */
				if(dValues[0] < dValues[1])
				{
					if(dValues[0] < dValues[2])
					{
						GI_VEC2_COPY(dVNew, dSimplex);
					}
					else
					{
						GI_VEC2_COPY(dVNew, dSimplex+4);
					}
				}
				else
				{
					if(dValues[1] < dValues[2])
					{
						GI_VEC2_COPY(dVNew, dSimplex+2);
					}
					else
					{
						GI_VEC2_COPY(dVNew, dSimplex+4);
					}
				}

				/* compute new matrix row */
				i = pParam->id;
				GISparseVector_to_zero(&row);
				system.bU[i] = system.bV[i] = 0.0;
				dSum = 0.0;
				pHalfEdge = pHEnd = pParam->vertex->hedge->twin;
				pPrev = pHalfEdge->twin->prev;
				dInvR0 = 1.0 / GIvec2d_dist(dVNew, pHalfEdge->pstart->params);
				dInvR1 = 1.0 / GIvec2d_dist(dVNew, pPrev->pstart->params);
				do
				{
					/* compute weights Wij = (tan(Aij/2) + tan(Bji/2)) / Rij */
					pNext = pHalfEdge->next->twin;
					GI_VEC2_SUB(v1, pPrev->pstart->params, dVNew);
					GI_VEC2_SUB(v0, pHalfEdge->pstart->params, dVNew);
					GI_VEC2_SUB(v2, pNext->pstart->params, dVNew);
					dInvR2 = 1.0 / GI_VEC2_LENGTH(v2);
					dCos = GI_VEC2_DOT(v0, v1) * dInvR0 * dInvR1;
					dTan1 = sqrt((1.0-dCos)/(1.0+dCos));
					dCos = GI_VEC2_DOT(v0, v2) * dInvR0 * dInvR2;
					dTan2 = sqrt((1.0-dCos)/(1.0+dCos));
					dCoord = (dTan1+dTan2) * dInvR0;
					dSum += dCoord;

					/* set coefficients */
					j = pHalfEdge->pstart->id;
					if(j < system.A->n)
						GISparseVector_append(&row, j, -dCoord);
					else
					{
						system.bU[i] += dCoord * pHalfEdge->pstart->params[0];
						system.bV[i] += dCoord * pHalfEdge->pstart->params[1];
					}
					pPrev = pHalfEdge;
					pHalfEdge = pNext;
					dInvR1 = dInvR0;
					dInvR0 = dInvR2;
				}while(pHalfEdge != pHEnd);
				GISparseVector_append(&row, i, dSum);
				for(j=0; j<row.size; ++j)
					system.A->values[ij++] = row.elements[j].value;
				uiEvalSum += uiEval;
			}
			pOldStretches[pParam->id] = pParam->stretch;
			pParam->stretch = 0.0;
		GI_LIST_NEXT(patch->params, pParam)
		dOldStretch = patch->stretch[uiMetric-GI_STRETCH_BASE];
		dOldMin = patch->min_param_stretch;
		dOldMax = patch->max_param_stretch;
		GIDebug(printf("evaluations per vertex: %5.2f\n", 
			(GIfloat)uiEvalSum/(GIfloat)uiPCount));
		uiEvalSum = 0;

		/* reparameterize */
		bSuccess = GILinearSystem_solve(&system);
		GILinearSystem_unknowns_to_params(&system);
		if(bSuccess)
			GIPatch_compute_stretch(patch, uiMetric, GI_TRUE, GI_FALSE);
		printf("stretch2 %f -> %f\n", (double)dOldStretch, (double)patch->stretch[uiMetric - GI_STRETCH_BASE]);
		count++;
	}while (dOldStretch - patch->stretch[uiMetric - GI_STRETCH_BASE] > 1e-4);

	/* restore last parameterization and stretch */
	GI_LIST_FOREACH(patch->params, pParam)
		if(!pParam->cut_hedge)
		{
			GI_VEC2_COPY(pParam->params, pOldParams+(pParam->id<<1));
		}
		pParam->stretch = pOldStretches[pParam->id];
	GI_LIST_NEXT(patch->params, pParam)
	patch->stretch[uiMetric-GI_STRETCH_BASE] = dOldStretch;
	patch->min_param_stretch = dOldMin;
	patch->max_param_stretch = dOldMax;

	/* clean up */
	GILinearSystem_destruct(&system);
	GI_FREE_ARRAY(pOldParams);
	GI_FREE_ARRAY(pOldStretches);
	GIHash_destruct(&hFaceAreas, sizeof(GIdouble));
	return GI_TRUE;
}

/** \internal
 *  \brief Gu's original iterative cutting and parameterization algorithm.
 *  \param par parameterizer to use
 *  \param patch patch to parameterize
 *  \retval GI_TRUE if parameterized successfully
 *  \retval GI_FALSE on error or abort
 *  \ingroup parameterization
 */
GIboolean GIParameterizer_gim(GIParameterizer *par, GIPatch *patch)
{
	GILinearSystem system;
	GIMesh *pMesh = patch->mesh;
	GICutPath *pPath, *pNewPath, *pNewPath2;
	GIFace *pFace;
	GIHalfEdge *pHalfEdge, *pHSource;
	GIVertex *pVertex, *pVStart, *pVEnd, *pVOther;
	GIParam *pParam, *pPNew, *pPSplit;
	GIParamSave *pSave;
	GIPathInfo *pInfo;
	GIParam *pOldCorners[4];
	GIint *pOldGLengths = NULL;
	GIdouble dOldStretch, dOldMin, dOldMax;
	GIdouble dDist, dOldDist, dMaxDist, dLength, dOldHLength, dOldSideLength;
	GIuint uiOldHStack, uiOldPStack, uiOldPCount, uiOldHCount, uiOldCutSplits;
	GIuint uiMetric = par->stretch_metric;
	GIboolean bSuccess;
	GIHeap qFringe;
	GIHash hPathInfos;
	GIHash hParamSaves;
	GIuint uiSide;

	/* create datastructures */
	GIHeap_construct(&qFringe, pMesh->ecount, lessd, -DBL_MAX, GI_TRUE);
	GIHash_construct(&hParamSaves, patch->pcount, 0.0f, sizeof(GIParam*), 
		hash_pointer, compare_pointer, copy_pointer);
	GIHash_construct(&hPathInfos, patch->pcount, 0.0f, sizeof(GIuint), 
		hash_uint, compare_uint, copy_uint);

	/* initial parameterization */
	bSuccess = GIParameterizer_stretch_minimizing(par, patch);


	do
	{
		/* save old parameterization and cut */
		GIDebug(printf("save old parameterization\n"));
		uiOldHStack = pMesh->split_hedges.size;
		uiOldPStack = patch->split_paths.size;
		uiOldCutSplits = pMesh->cut_splits;
		uiOldPCount = patch->pcount;
		uiOldHCount = patch->hcount;
		dOldHLength = patch->hlength;
		pOldGLengths = (GIint*)GI_REALLOC_ARRAY(pOldGLengths, patch->groups, sizeof(GIint));
		GI_LIST_FOREACH(patch->paths, pPath)
			pOldGLengths[pPath->group] = pPath->glength;
		GI_LIST_NEXT(patch->paths, pPath)
		GI_LIST_FOREACH(patch->params, pParam)
			pSave = (GIParamSave*)GIHash_find(&hParamSaves, &pParam);
			if(!pSave)
			{
				pSave = (GIParamSave*)GI_MALLOC_SINGLE(sizeof(GIParamSave));
				GIHash_insert(&hParamSaves, &pParam, pSave);
			}
			pSave->id = pParam->id;
			GI_VEC2_COPY(pSave->params, pParam->params);
			pSave->stretch = pParam->stretch;
		GI_LIST_NEXT(patch->params, pParam)
		dOldStretch = patch->stretch[uiMetric-GI_STRETCH_BASE];
		dOldMin = patch->min_param_stretch;
		dOldMax = patch->max_param_stretch;
		if(!patch->fixed_corners)
			memcpy(pOldCorners, patch->corners, 4*sizeof(GIParam*));

		/* conformal parameterization on circle */
		GIDebug(printf("conformal parameterization\n"));
		GIParameterizer_arc_length_circle(par, patch);
		GILinearSystem_construct(&system, par, patch, GI_MEAN_VALUE, GI_FALSE, GI_FALSE);
		GILinearSystem_solve(&system);
		GILinearSystem_unknowns_to_params(&system);
		GILinearSystem_destruct(&system);

		/* find new cut node to connect */
		GIDebug(printf("find max stretch vertex\n"));
		pFace = GIPatch_compute_stretch(patch, uiMetric, GI_FALSE, GI_FALSE);
		pVStart = pFace->hedges->vstart;
		dMaxDist = 0.0;
		GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
			pParam = pHalfEdge->pstart;
			if(!pParam->cut_hedge)
			{
				dDist = GI_VEC2_LENGTH_SQR(pParam->params);
				if(dDist > dMaxDist)
				{
					dMaxDist = dDist;
					pVStart = pHalfEdge->vstart;
				}
			}
		GI_LIST_NEXT(pFace->hedges, pHalfEdge)

		/* find shortest path (Dijkstra) */
		GIDebug(printf("dijkstra\n"));
		GIHeap_clear(&qFringe);
		GIHash_clear(&hPathInfos, sizeof(GIPathInfo));
		GIHash_insert(&hPathInfos, &pVStart->id, 
			GI_CALLOC_SINGLE(sizeof(GIPathInfo)));
		pVertex = pVStart;
		dOldDist = 0.0;
		while(!pVertex->cut_degree)
		{
			pHalfEdge = pVertex->hedge;
			do
			{
				pVOther = pHalfEdge->next->vstart;
				dDist = dOldDist + pHalfEdge->edge->length;
				pInfo = (GIPathInfo*)GIHash_find(&hPathInfos, &pVOther->id);
				if(!pInfo || dDist < pInfo->distance)
				{
					if(!pInfo)
						GIHash_insert(&hPathInfos, &pVOther->id, pInfo=
							(GIPathInfo*)GI_MALLOC_SINGLE(sizeof(GIPathInfo)));
					pInfo->source = pHalfEdge;
					pInfo->distance = dDist;
					GIHeap_enqueue(&qFringe, pVOther, dDist);
				}
				pHalfEdge = pHalfEdge->twin->next;
			}while(pHalfEdge != pVertex->hedge);
			pVertex = (GIVertex*)GIHeap_dequeue(&qFringe, &dOldDist);
		}
		pVEnd = pVertex;
		++pVEnd->cut_degree;

		/* find path to split */
		pHSource = ((GIPathInfo*)GIHash_find(&hPathInfos, pVEnd))->source;
		pPSplit = pHSource->next->pstart;
		pParam = patch->params;
		uiSide = 0;
		GI_LIST_FOREACH(patch->paths, pPath)
			dLength = 0.0;
			do
			{
				dLength += pParam->cut_hedge->edge->length;
				pParam = pParam->cut_hedge->next->pstart;
				if(pParam == pPSplit)
				{
					if(GIPatch_split_path(patch, pPath, pParam))
					{
						pPath->next->elength = pPath->elength - dLength;
						pPath->elength = dLength;
						if(pPath->twin)
						{
							pPath->twin->elength = pPath->elength;
							pPath->next->twin->elength = pPath->next->elength;
						}
					}
					break;
				}
				if(pParam == patch->corners[uiSide+1])
					++uiSide;
			}while(pParam != pPath->next->pstart);
			if(pParam == pPSplit)
				break;
		GI_LIST_NEXT(patch->paths, pPath)
		dOldSideLength = patch->side_lengths[uiSide];

		/* create new paths and connect to existing cut */
		pNewPath = (GICutPath*)GI_MALLOC_PERSISTENT(sizeof(GICutPath));
		pNewPath2 = (GICutPath*)GI_MALLOC_PERSISTENT(sizeof(GICutPath));
		GI_LIST_INSERT(patch->paths, pPath->next, pNewPath2);
		GI_LIST_INSERT(patch->paths, pNewPath2, pNewPath);
		pNewPath->id = patch->path_count++;
		pNewPath2->id = patch->path_count++;
		pNewPath->patch = pNewPath2->patch = patch;
		pNewPath->group = pNewPath2->group = patch->groups++;
		pNewPath->glength = pNewPath2->glength = 0;
		pNewPath->elength = 0.0;
		pNewPath->twin = pNewPath2;
		pNewPath2->twin = pNewPath;
		pNewPath2->pstart = pVStart->hedge->pstart;
		pNewPath->pstart = pPNew = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
		GI_LIST_ADD(patch->params, pPNew);
		pPNew->id = patch->pcount++;
		GI_VEC2_COPY(pPNew->params, pPSplit->params);
		pPNew->vertex = pVEnd;
		pPNew->cut_hedge = pHSource->twin;
		pHalfEdge = pHSource;
		while(pHalfEdge->twin->face && pHalfEdge != pHalfEdge->pstart->cut_hedge)
		{
			pHalfEdge->twin->pstart = pPNew;
			pHalfEdge = pHalfEdge->twin->prev;
		}

		/* assemble path */
		GIDebug(printf("assemble cut path\n"));
		pVertex = pHSource->vstart;
		while(pVertex != pVStart)
		{
			patch->hcount += 2;
			pNewPath->elength += pHSource->edge->length;
			pVertex->cut_degree = 2;
			pParam = pHSource->pstart;
			pParam->cut_hedge = pHSource;
			pHalfEdge = pHSource;
			pHSource = ((GIPathInfo*)GIHash_find(&hPathInfos, &pVertex->id))->source;
			pPNew = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
			GI_LIST_ADD(patch->params, pPNew);
			pPNew->id = patch->pcount++;
			GI_VEC2_COPY(pPNew->params, pParam->params);
			pPNew->vertex = pVertex;
			while(pHalfEdge != pHSource->twin)
			{
				pHalfEdge = pHalfEdge->twin->next;
				pHalfEdge->pstart = pPNew;
			}
			pPNew->cut_hedge = pHalfEdge;
			pVertex = pHSource->vstart;
		}
		pVStart->cut_degree = 1;
		pHSource->pstart->cut_hedge = pHSource;
		patch->hcount += 2;
		pNewPath->elength += pHSource->edge->length;
		pNewPath2->elength = pNewPath->elength;
		patch->hlength += 2.0 * pNewPath->elength;
		if(patch->fixed_corners)
			patch->side_lengths[uiSide] += 2.0 * pNewPath->elength;
		pMesh->varray_attribs = 0;

		/* reparameterize */
		GIDebug(printf("stretch minimization\n"));
		GIPatch_prevent_singularities(patch);
		GIPatch_renumerate_params(patch);
		pMesh->cut_splits = pMesh->split_hedges.size;
		patch->resolution = 0;
		GIboolean bSuccess1 = GIParameterizer_arc_length_square(par, patch);
		GIboolean bSuccess2 = GIParameterizer_stretch_minimizing(par, patch);

		bSuccess1 = (bSuccess1 && bSuccess2);
		printf("stretch %f => %f\n", (double)dOldStretch, (double)patch->stretch[uiMetric - GI_STRETCH_BASE]);
		if(!bSuccess)
			break;
	}while(patch->stretch[uiMetric-GI_STRETCH_BASE] < dOldStretch);

	/* reverse edge splits */
	GIMesh_revert_splits(pMesh, pMesh->split_hedges.size-uiOldHStack);
	pMesh->cut_splits = uiOldCutSplits;

	/* destroy new path */
	pVertex = pVStart;
	pParam = pVStart->hedge->pstart;
	pHalfEdge = pParam->cut_hedge;
	while(pVertex != pVEnd)
	{
		pVertex->cut_degree = 0;
		pParam->cut_hedge = NULL;
		pParam = pHalfEdge->next->pstart;
		pVertex = pParam->vertex;
		pHalfEdge = pHalfEdge->twin;
		pPNew = pHalfEdge->pstart;
		while(pHalfEdge->pstart == pPNew)
		{
			pHalfEdge->pstart = pParam;
			pHalfEdge = pHalfEdge->prev->twin;
		}
		GI_LIST_DELETE_PERSISTENT(patch->params, pPNew, sizeof(GIParam));
		pHalfEdge = pParam->cut_hedge;
	}
	--pVEnd->cut_degree;
	patch->pcount = uiOldPCount;
	patch->hcount = uiOldHCount;
	patch->hlength = dOldHLength;
	if(patch->fixed_corners)
		patch->side_lengths[uiSide] = dOldSideLength;
	else
		memcpy(patch->corners, pOldCorners, 4*sizeof(GIParam*));

	/* delete paths and reverse path split */
	GI_LIST_DELETE_PERSISTENT(patch->paths, pNewPath, sizeof(GICutPath));
	GI_LIST_DELETE_PERSISTENT(patch->paths, pNewPath2, sizeof(GICutPath));
	patch->path_count -= 2;
	--patch->groups;
	GIPatch_revert_splits(patch, patch->split_paths.size-uiOldPStack);

	/* restore grid lengths */
	GI_LIST_FOREACH(patch->paths, pPath)
		pPath->glength = pOldGLengths[pPath->group];
	GI_LIST_NEXT(patch->paths, pPath)
	GI_FREE_ARRAY(pOldGLengths);

	/* restore parameterization */
	GI_LIST_FOREACH(patch->params, pParam)
		pSave = (GIParamSave*)GIHash_remove(&hParamSaves, &pParam);
		pParam->id = pSave->id;
		GI_VEC2_COPY(pParam->params, pSave->params);
		pParam->stretch = pSave->stretch;
		GI_FREE_SINGLE(pSave, sizeof(GIParamSave));
	GI_LIST_NEXT(patch->params, pParam)
	patch->stretch[uiMetric-GI_STRETCH_BASE] = dOldStretch;
	patch->min_param_stretch = dOldMin;
	patch->max_param_stretch = dOldMax;

	/* clean up */
	GIHeap_destruct(&qFringe);
	GIHash_destruct(&hParamSaves, sizeof(GIParamSave));
	GIHash_destruct(&hPathInfos, sizeof(GIPathInfo));
	return bSuccess;
}

/** \internal
 *  \brief Linear system constructor.
 *  \param system system to construct
 *  \param par parameterizer to use
 *  \param patch patch to construct system for
 *  \param type mapping type to construct system for
 *  \param force_non_symmetric GI_TRUE to fully store symmetric matrices
 *  \param store_rhs GI_TRUE to store boundary coefficients in extra sparse matrix
 *  \ingroup parameterization
 */
void GILinearSystem_construct(GILinearSystem *system, GIParameterizer *par, 
							  GIPatch *patch, GIenum type, 
							  GIboolean force_non_symmetric, GIboolean store_rhs)
{
	GISparseMatrixLIL A;
	GISparseMatrixLIL *B = NULL;
	GIdouble *bU, *bV;
	GIuint N = patch->pcount - patch->hcount;
	GIHalfEdge *pHalfEdge, *pEnd, *pPrev, *pNext, *pWork;
	GIVertex *pVertex;
	GIParam *pParam;
	GIdouble v0[3], v1[3], v2[3], v3[3], v4[3];
	GIdouble *vec;
	GIdouble dCos, dTan1, dTan2, dCot1, dCot2, dInvR, dCoord, dSum;
	GIuint i, j, k;

	/* create system */
	system->parameterizer = par;
	system->patch = patch;
	system->A = (GISparseMatrixCSR*)GI_MALLOC_SINGLE(sizeof(GISparseMatrixCSR));
	bU = system->bU = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	bV = system->bV = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	system->u = (GIdouble*)GI_CALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	system->v = (GIdouble*)GI_CALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	if(store_rhs)
	{
		B = system->B = (GISparseMatrixLIL*)GI_MALLOC_SINGLE(sizeof(GISparseMatrixLIL));
		GISparseMatrixLIL_construct(B, N, GI_FALSE);
	}
	else
		system->B = NULL;

	/* compute coefficients */
	switch(type)
	{
	case GI_TUTTE_BARYCENTRIC:
		GISparseMatrixLIL_construct(&A, N, !force_non_symmetric);
		GI_LIST_FOREACH(patch->params, pParam)
			if(!pParam->cut_hedge)
			{
				i = pParam->id;
				bU[i] = bV[i] = 0.0;
				dSum = 0.0;
				pHalfEdge = pEnd = pParam->vertex->hedge->twin;
				do
				{
					/* compute weight Wij = 1 and set coefficients */
					dSum += 1.0;
					j = pHalfEdge->pstart->id;
					if(j < N)
						GISparseMatrixLIL_append(&A, i, j, -1.0);
					else
					{
						vec = pHalfEdge->pstart->params;
						if(B)
							GISparseMatrixLIL_append(B, i, j, 1.0);
						bU[i] += vec[0];
						bV[i] += vec[1];
					}
					pHalfEdge = pHalfEdge->next->twin;
				}while(pHalfEdge != pEnd);
				GISparseMatrixLIL_append(&A, i, i, dSum);
			}
		GI_LIST_NEXT(patch->params, pParam)
		break;
	case GI_SHAPE_PRESERVING:
		{
		GIHash hLocal;
		GILocalInfo *pLocal;
		GIParam *pPk[3];
		GIdouble *p0, *p1, *p2;
		GIdouble dTheta, dDet, dMu[3], d1Sum;
		GISparseMatrixLIL_construct(&A, N, GI_FALSE);
		GIHash_construct(&hLocal, 16, 0.0f, sizeof(GIuint), 
			hash_uint, compare_uint, copy_uint);
		GI_LIST_FOREACH(patch->params, pParam)
			if(!pParam->cut_hedge)
			{
				/* compute sum of surface angles */
				pVertex = pParam->vertex;
				pPrev = pEnd = pVertex->hedge->twin;
				pHalfEdge = pPrev->twin->prev;
				dTheta = 0.0;
				while(pHalfEdge != pEnd)
				{
					GI_VEC3_SUB(v0, pPrev->vstart->coords, pVertex->coords);
					GI_VEC3_SUB(v1, pHalfEdge->vstart->coords, pVertex->coords);
					dTheta += acos(GI_VEC3_DOT(v0, v1)/
						(pPrev->edge->length*pHalfEdge->edge->length));
					pLocal = (GILocalInfo*)GI_MALLOC_SINGLE(sizeof(GILocalInfo));
					pLocal->angle = dTheta;
					GIHash_insert(&hLocal, pHalfEdge->pstart, pLocal);
					pPrev = pHalfEdge;
					pHalfEdge = pHalfEdge->twin->prev;
				}
				GI_VEC3_SUB(v0, pHalfEdge->vstart->coords, pVertex->coords);
				dTheta = GI_TWO_PI / (dTheta+acos(GI_VEC3_DOT(v0, v1)/
					(pPrev->edge->length*pHalfEdge->edge->length)));

				/* compute local parameterization */
				pLocal = (GILocalInfo*)GI_CALLOC_SINGLE(sizeof(GILocalInfo));
				pLocal->param[0] = pHalfEdge->edge->length;
				GIHash_insert(&hLocal, pHalfEdge->pstart, pLocal);
				pHalfEdge = pHalfEdge->twin->prev;
				while(pHalfEdge != pEnd)
				{
					pLocal = (GILocalInfo*)GIHash_find(&hLocal, pHalfEdge->pstart);
					pLocal->angle *= dTheta;
					GI_VEC2_SET(pLocal->param, pHalfEdge->edge->length*cos(pLocal->angle), 
						pHalfEdge->edge->length*sin(pLocal->angle));
					pHalfEdge = pHalfEdge->twin->prev;
				}

				/* compute weights */
				i = pParam->id;
				bU[i] = bV[i] = 0.0;
				dSum = 0.0;
				do
				{
					/* compute barycentric coordinates for enclosing triangles */
					pPk[0] = pHalfEdge->pstart;
					p0 = ((GILocalInfo*)GIHash_find(&hLocal, pPk[0]))->param;
					pWork = pHalfEdge->twin->prev;
					pPk[2] = pWork->pstart;
					p2 = ((GILocalInfo*)GIHash_find(&hLocal, pPk[2]))->param;
					unsigned int count = 0;
					do
					{
						pPk[1] = pPk[2];
						p1 = p2;
						pWork = pWork->twin->prev;
						pPk[2] = pWork->pstart;
						p2 = ((GILocalInfo*)GIHash_find(&hLocal, pPk[2]))->param;
						dDet = GI_VEC2_DET(p0, p2);
						count++;
						if (count == 1000000) break;
					}while(dDet >= 0.0);
					if (count == 1000000)
					{
						printf("Waring:三角形を囲む重心座標計算\n");
					}
					dMu[0] = GI_VEC2_DET(p1, p2);
					dMu[1] = GI_VEC2_DET(p2, p0);
					dMu[2] = GI_VEC2_DET(p0, p1);
					d1Sum = 1.0 / (dMu[0]+dMu[1]+dMu[2]);
					GI_VEC3_SCALE(dMu, dMu, d1Sum);

					/* add coefficients */
					for(k=0; k<3; ++k)
					{
						dCoord = dMu[k];
						dSum += dCoord;
						j = pPk[k]->id;
						if(j < N)
							GISparseMatrixLIL_add(&A, i, j, -dCoord);
						else
						{
							vec = pPk[k]->params;
							if(B)
								GISparseMatrixLIL_add(B, i, j, dCoord);
							bU[i] += dCoord * vec[0];
							bV[i] += dCoord * vec[1];
						}
					}
					pHalfEdge = pHalfEdge->twin->prev;
				}while(pHalfEdge != pEnd);
				GISparseMatrixLIL_append(&A, i, i, dSum);
				GIHash_clear(&hLocal, sizeof(GILocalInfo));
			}
		GI_LIST_NEXT(patch->params, pParam)
		GIHash_destruct(&hLocal, sizeof(GILocalInfo));
		}
		break;
	case GI_DISCRETE_HARMONIC:
		GISparseMatrixLIL_construct(&A, N, !force_non_symmetric);
		GI_LIST_FOREACH(patch->params, pParam)
			if(!pParam->cut_hedge)
			{
				i = pParam->id;
				bU[i] = bV[i] = 0.0;
				dSum = 0.0;
				pVertex = pParam->vertex;
				pHalfEdge = pEnd = pVertex->hedge->twin;
				pPrev = pHalfEdge->twin->prev;
				do
				{
					/* compute weights Wij = cot(Yij) + cot(Yji) */
					pNext = pHalfEdge->next->twin;
					GI_VEC3_SUB(v0, pVertex->coords, pPrev->vstart->coords);
					GI_VEC3_SUB(v1, pHalfEdge->vstart->coords, pPrev->vstart->coords);
					GI_VEC3_SUB(v2, pVertex->coords, pNext->vstart->coords);
					GI_VEC3_SUB(v3, pHalfEdge->vstart->coords, pNext->vstart->coords);
					dCos = GI_VEC3_DOT(v0, v1) / (pPrev->edge->length*pPrev->prev->edge->length);
					dCot1 = dCos / sqrt(1.0-dCos*dCos);
					dCos = GI_VEC3_DOT(v2, v3) / (pNext->edge->length*pHalfEdge->prev->edge->length);
					dCot2 = dCos / sqrt(1.0-dCos*dCos);
					dCoord = dCot1 + dCot2;
					dSum += dCoord;

					/* set coefficients */
					j = pHalfEdge->pstart->id;
					if(j < N)
						GISparseMatrixLIL_append(&A, i, j, -dCoord);
					else
					{
						vec = pHalfEdge->pstart->params;
						if(B)
							GISparseMatrixLIL_append(B, i, j, dCoord);
						bU[i] += dCoord * vec[0];
						bV[i] += dCoord * vec[1];
					}
					pPrev = pHalfEdge;
					pHalfEdge = pNext;
				}while(pHalfEdge != pEnd);
				GISparseMatrixLIL_append(&A, i, i, dSum);
			}
		GI_LIST_NEXT(patch->params, pParam)
		break;
	case GI_MEAN_VALUE:
		GISparseMatrixLIL_construct(&A, N, GI_FALSE);
		GI_LIST_FOREACH(patch->params, pParam)
			if(!pParam->cut_hedge)
			{
				i = pParam->id;
				bU[i] = bV[i] = 0.0;
				dSum = 0.0;
				pVertex = pParam->vertex;
				pHalfEdge = pEnd = pVertex->hedge->twin;
				pPrev = pHalfEdge->twin->prev;
				do
				{
					/* compute weights Wij = (tan(Aij/2) + tan(Bji/2)) / Rij */
					pNext = pHalfEdge->next->twin;
					GI_VEC3_SUB(v1, pPrev->vstart->coords, pVertex->coords);
					GI_VEC3_SUB(v0, pHalfEdge->vstart->coords, pVertex->coords);
					GI_VEC3_SUB(v2, pNext->vstart->coords, pVertex->coords);
					dCos = GI_VEC3_DOT(v0, v1) / 
						(pHalfEdge->edge->length*pPrev->edge->length);
					dTan1 = sqrt((1.0-dCos)/(1.0+dCos));
					dCos = GI_VEC3_DOT(v0, v2) / 
						(pHalfEdge->edge->length*pNext->edge->length);
					dTan2 = sqrt((1.0-dCos)/(1.0+dCos));
					dCoord = (dTan1+dTan2) * 
						(patch->mesh->mean_edge/pHalfEdge->edge->length);
					dSum += dCoord;

					/* set coefficients */
					j = pHalfEdge->pstart->id;
					if(j < N)
						GISparseMatrixLIL_append(&A, i, j, -dCoord);
					else
					{
						vec = pHalfEdge->pstart->params;
						if(B)
							GISparseMatrixLIL_append(B, i, j, dCoord);
						bU[i] += dCoord * vec[0];
						bV[i] += dCoord * vec[1];
					}
					pPrev = pHalfEdge;
					pHalfEdge = pNext;
				}while(pHalfEdge != pEnd);
				GISparseMatrixLIL_append(&A, i, i, dSum);
			}
		GI_LIST_NEXT(patch->params, pParam)
		break;
	case GI_DISCRETE_AUTHALIC:
		GISparseMatrixLIL_construct(&A, N, GI_FALSE);
		GI_LIST_FOREACH(patch->params, pParam)
			if(!pParam->cut_hedge)
			{
				i = pParam->id;
				bU[i] = bV[i] = 0.0;
				dSum = 0.0;
				pVertex = pParam->vertex;
				pHalfEdge = pEnd = pVertex->hedge->twin;
				pPrev = pHalfEdge->twin->prev;
				do
				{
					/* compute weight Wij = (cot(Aji) + cot(Bij)) / Rij^2 */
					pNext = pHalfEdge->next->twin;
					GI_VEC3_SUB(v1, pPrev->vstart->coords, pHalfEdge->vstart->coords);
					GI_VEC3_SUB(v0, pVertex->coords, pHalfEdge->vstart->coords);
					GI_VEC3_SUB(v2, pNext->vstart->coords, pHalfEdge->vstart->coords);
					dCos = GI_VEC3_DOT(v0, v1) / 
						(pHalfEdge->edge->length*pPrev->prev->edge->length);
					dCot1 = dCos / sqrt(1.0-dCos*dCos);
					dCos = GI_VEC3_DOT(v0, v2) / 
						(pHalfEdge->edge->length*pHalfEdge->prev->edge->length);
					dCot2 = dCos / sqrt(1.0-dCos*dCos);
					dInvR = patch->mesh->mean_edge / pHalfEdge->edge->length;
					dCoord = (dCot1+dCot2) * dInvR * dInvR;
					dSum += dCoord;

					/* set coefficients */
					j = pHalfEdge->pstart->id;
					if(j < N)
						GISparseMatrixLIL_append(&A, i, j, -dCoord);
					else
					{
						vec = pHalfEdge->pstart->params;
						if(B)
							GISparseMatrixLIL_append(B, i, j, dCoord);
						bU[i] += dCoord * vec[0];
						bV[i] += dCoord * vec[1];
					}
					pPrev = pHalfEdge;
					pHalfEdge = pNext;
				}while(pHalfEdge != pEnd);
				GISparseMatrixLIL_append(&A, i, i, dSum);
			}
		GI_LIST_NEXT(patch->params, pParam)
		break;
	case GI_INTRINSIC:
		{
		GIdouble dLambda = par->conformal_weight;
		GIdouble dMu = par->authalic_weight;
		GISparseMatrixLIL_construct(&A, N, GI_FALSE);
		GI_LIST_FOREACH(patch->params, pParam)
			if(!pParam->cut_hedge)
			{
				i = pParam->id;
				bU[i] = bV[i] = 0.0;
				dSum = 0.0;
				pVertex = pParam->vertex;
				pHalfEdge = pEnd = pVertex->hedge->twin;
				pPrev = pHalfEdge->twin->prev;
				do
				{
					/* compute weights Wij = lambda*(cot(Yij) + cot(Yji)) + mu*((cot(Aji) + cot(Bij)) / Rij^2) */
					pNext = pHalfEdge->next->twin;
					GI_VEC3_SUB(v0, pVertex->coords, pPrev->vstart->coords);
					GI_VEC3_SUB(v1, pHalfEdge->vstart->coords, pPrev->vstart->coords);
					GI_VEC3_SUB(v2, pVertex->coords, pNext->vstart->coords);
					GI_VEC3_SUB(v3, pHalfEdge->vstart->coords, pNext->vstart->coords);
					GI_VEC3_SUB(v4, pVertex->coords, pHalfEdge->vstart->coords);
					dCos = GI_VEC3_DOT(v0, v1) / (pPrev->edge->length*pPrev->prev->edge->length);
					dCot1 = dCos / sqrt(1.0-dCos*dCos);
					dCos = GI_VEC3_DOT(v2, v3) / (pNext->edge->length*pHalfEdge->prev->edge->length);
					dCot2 = dCos / sqrt(1.0-dCos*dCos);
					dCoord = dLambda * (dCot1+dCot2);
					GI_VEC3_NEGATE(v1, v1);
					GI_VEC3_NEGATE(v3, v3);
					dCos = GI_VEC3_DOT(v4, v1) / 
						(pHalfEdge->edge->length*pPrev->prev->edge->length);
					dCot1 = dCos / sqrt(1.0-dCos*dCos);
					dCos = GI_VEC3_DOT(v4, v3) / 
						(pHalfEdge->edge->length*pHalfEdge->prev->edge->length);
					dCot2 = dCos / sqrt(1.0-dCos*dCos);
					dInvR = patch->mesh->mean_edge / pHalfEdge->edge->length;
					dCoord += dMu * (dCot1+dCot2) * dInvR * dInvR;
					dSum += dCoord;

					/* set coefficients */
					j = pHalfEdge->pstart->id;
					if(j < N)
						GISparseMatrixLIL_append(&A, i, j, -dCoord);
					else
					{
						vec = pHalfEdge->pstart->params;
						if(B)
							GISparseMatrixLIL_append(B, i, j, dCoord);
						bU[i] += dCoord * vec[0];
						bV[i] += dCoord * vec[1];
					}
					pPrev = pHalfEdge;
					pHalfEdge = pNext;
				}while(pHalfEdge != pEnd);
				GISparseMatrixLIL_append(&A, i, i, dSum);
			}
		GI_LIST_NEXT(patch->params, pParam)
		}
		break;
	default:
		GISparseMatrixLIL_construct(&A, N, GI_FALSE);
		memset(bU, 0, N*sizeof(GIdouble));
		memset(bV, 0, N*sizeof(GIdouble));
	}

	/* convert matrix */
	GISparseMatrixCSR_construct(system->A, &A);
	GISparseMatrixLIL_destruct(&A);
}

/** \internal
 *  \brief Linear system destructor.
 *  \param system system to destruct
 *  \ingroup parameterization
 */
void GILinearSystem_destruct(GILinearSystem *system)
{
	/* clean up */
	if(system->A)
	{
		GISparseMatrixCSR_destruct(system->A);
		GI_FREE_SINGLE(system->A, sizeof(GISparseMatrixCSR));
	}
	if(system->B)
	{
		GISparseMatrixLIL_destruct(system->B);
		GI_FREE_SINGLE(system->B, sizeof(GISparseMatrixLIL));
	}
	if(system->bU)
		GI_FREE_ALIGNED(system->bU);
	if(system->bV)
		GI_FREE_ALIGNED(system->bV);
	if(system->u)
		GI_FREE_ALIGNED(system->u);
	if(system->v)
		GI_FREE_ALIGNED(system->v);
	memset(system, 0, sizeof(GILinearSystem));
}

/** \internal
 *  \brief Copy unknown vectors to params.
 *  \param system system to work on
 *  \ingroup parameterization
 */
void GILinearSystem_unknowns_to_params(GILinearSystem *system)
{
	GIParam *pParam, *pPHead = system->patch->params;
	GIdouble *u = system->u;
	GIdouble *v = system->v;

	/* copy solution to params */
	GI_LIST_FOREACH(pPHead, pParam)
		if(!pParam->cut_hedge)
		{
			pParam->params[0] = u[pParam->id];
			pParam->params[1] = v[pParam->id];
		}
	GI_LIST_NEXT(pPHead, pParam)
	if(!system->patch->parameterized)
	{
		++system->patch->mesh->param_patches;
		system->patch->parameterized = GI_TRUE;
	}
	memset(system->patch->stretch, 0, GI_STRETCH_COUNT*sizeof(GIdouble));
	system->patch->param_metric = 0;
}

/** \internal
 *  \brief Solve linear system.
 *  \param system system to solve
 *  \retval GI_TRUE if solved successfully
 *  \retval GI_FALSE if system could not be solved
 *  \ingroup parameterization
 */
GIboolean GILinearSystem_solve(GILinearSystem *system)
{
	GISolverData dataU;
	GIboolean bBICGSTAB = (system->parameterizer->solver == GI_SOLVER_BICGSTAB);
	GIsolverfunc pfnUnsymmetric = (bBICGSTAB ? GISolver_bicgstab : GISolver_gmres);
	GIuint uiMaxIter = ((system->A->symmetric || bBICGSTAB) ? 13 : 3) * sqrt((GIdouble)system->A->n);
	GIuint uiIterU, uiIterV;

	/* assemble configuration */
	GISparseMatrixCSR_prepare_ilu(system->A);
	dataU.solver_func = (system->A->symmetric ? GISolver_cg : pfnUnsymmetric);
	dataU.A = (GISparseMatrix*)system->A;
	dataU.b = system->bU;
	dataU.x = system->u;
	dataU.ax = GISparseMatrixCSR_ax;
	dataU.pc = GISparseMatrixCSR_pc_ilu;
	dataU.eps = 1e-6;
	dataU.max_iter = (dataU.solver_func==GISolver_gmres ? ((uiMaxIter<<8)|25) : uiMaxIter);

#if OPENGI_NUM_THREADS > 1
	if(system->parameterizer->context->use_threads)
	{
		/* start solver threads and wait for completion */
		GIthread threadU = GIthread_create(GILinearSystem_solve_thread, &dataU);
		uiIterV = dataU.solver_func((GISparseMatrix*)system->A, system->bV, system->v, 
			dataU.ax, dataU.pc, dataU.eps, dataU.max_iter);
		uiIterU = (GIuint)GIthread_join(threadU);
	}
	else
#endif
	{
		/* solve systems sequentially */
#pragma omp parallel
#pragma omp sections
		{
#pragma omp section
			{
				uiIterU = dataU.solver_func((GISparseMatrix*)system->A, system->bU, system->u,
					dataU.ax, dataU.pc, dataU.eps, dataU.max_iter);
			}
#pragma omp section
			{
				uiIterV = dataU.solver_func((GISparseMatrix*)system->A, system->bV, system->v,
					dataU.ax, dataU.pc, dataU.eps, dataU.max_iter);
			}
		}
	}

	/* check results */
	GIDebug(printf("iterations: %d , %d (%d)\n", uiIterU, uiIterV, uiMaxIter));
	if(uiIterU > dataU.max_iter || uiIterV > dataU.max_iter)
	{
		GIContext_error(system->parameterizer->context, GI_NUMERICAL_ERROR);
		return GI_FALSE;
	}
	return GI_TRUE;
}

/** \internal
 *  \brief Thread execution function for solving linear system.
 *  \param arg thread parameters
 *  \return number of iterations for solving linear system
 *  \ingroup parameterization
 */
GIthreadret GITHREADENTRY GILinearSystem_solve_thread(GIvoid *arg)
{
	/* solve system */
	GISolverData *pData = (GISolverData*)arg;
	return (GIthreadret)pData->solver_func(pData->A, pData->b, 
		pData->x, pData->ax, pData->pc, pData->eps, pData->max_iter);
}
