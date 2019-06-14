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
 *  \brief Implementation of structures and functions for mesh cutting.
 */

#include "gi_cutter.h"
#include "gi_context.h"
#include "gi_memory.h"
#include "gi_multiresolution.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

#define GI_EDGE_WHITE			0x00		/**< Edge not processed (on mesh boundary or adjacent to 2 triangles). */
#define GI_EDGE_GREY			0x01		/**< Edge in processing (adjacent to exactly 1 triangle). */
#define GI_EDGE_BLACK			0x02		/**< Edge processed (removed, not on cut). */
#define GI_EDGE_FACE0			0x04		/**< Remove triangle of hedge[0] (rather than hedge[1]). */
#define GI_EDGE_USEABLE			0x10		/**< Edge may be on cut (or removed in further step). */
#define GI_EDGE_USED			0x30		/**< Edge surely on cut (on mesh boundary). */

/** \internal
 *  \brief Information about the faces belonging to a patch.
 */
typedef struct _GIPatchFaces
{
	GIuint					fcount;					/**< Number of faces in patch. */
	GIuint					hcount;					/**< Number of cut edges. */
	GIFace					*faces;					/**< Face list. */
	struct _GIPatchFaces	*next;					/**< Next info in list. */
} GIPatchFaces;

/** \internal
 *  \brief Information about path during A* or Dijkstra.
*/
typedef struct _GIPathInfo
{
	GIHalfEdge	*source;					/**< Source half edge. */
	GIdouble	distance;					/**< Distance to root. */
} GIPathInfo;

typedef struct _GIPatchBoundary
{
	GIHalfEdge	*hedge;
	GIVertex	*median;
} GIPatchBoundary;

typedef struct _GIClusterPatch
{
	GIuint			bcount;
	GIPatchBoundary	*boundaries;
	GIVertex		*root;
} GIClusterPatch;

typedef struct _GIEdgeAngle
{
	GIEdge		*edge;
	GIdouble	angle;
} GIEdgeAngle;

/** \internal
 *  \brief Compare edges for qsort.
 *  \param a first value
 *  \param b second value
 *  \return negative value if a < b, positive value if a > b and 0 if equal
 */
static int compare(const void *a, const void *b)
{
	if(((const GIEdgeAngle*)a)->angle < ((const GIEdgeAngle*)b)->angle)
		return -1;
	if(((const GIEdgeAngle*)b)->angle < ((const GIEdgeAngle*)a)->angle)
		return 1;
	return 0;
}
#if 0
/** Set boolean configuration parameter of cutting.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup cutting
 */
void GIAPIENTRY giCutterParameterb(GIenum pname, GIboolean param)
{
	GICutter *pCutter = &(GIContext_current()->cutter);

	/* select state and set value */
	switch(pname)
	{
	case GI_STRAIGHTEN_CUT:
		pCutter->straighten = param;
		break;
	default:
		GIContext_error(pCutter->context, GI_INVALID_ENUM);
	}
}
#endif
/** Set integer configuration parameter of cutting.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup cutting
 */
void GIAPIENTRY giCutterParameteri(GIenum pname, GIint param)
{
	GICutter *pCutter = &(GIContext_current()->cutter);

	/* select state and set value */
	switch(pname)
	{
	case GI_CUTTER:
		switch(param)
		{
		case GI_INITIAL_GIM:
		case GI_CATMULL_CLARK_SUBDIVISION:
//		case GI_FACE_CLUSTERING:
			pCutter->cutter = param;
			break;
		default:
			GIContext_error(pCutter->context, GI_INVALID_ENUM);
		}
		break;
	case GI_SUBDIVISION_ITERATIONS:
		if(param > 0)
			pCutter->iterations = param;
		else
			GIContext_error(pCutter->context, GI_INVALID_ENUM);
		break;
	default:
		GIContext_error(pCutter->context, GI_INVALID_ENUM);
	}
}
#if 0
/** Set floating point configuration parameter of cutting.
 *  \param pname state to set
 *  \param param value to set
 *  \ingroup cutting
 */
void GIAPIENTRY giCutterParameterf(GIenum pname, GIfloat param)
{
	GICutter *pCutter = &(GIContext_current()->cutter);

	/* select state and set value */
	switch(pname)
	{
	case GI_ORIENTATION_WEIGHT:
		pCutter->orientation_weight = param;
		break;
	case GI_SHAPE_WEIGHT:
		pCutter->shape_weight = param;
		break;
	default:
		GIContext_error(pCutter->context, GI_INVALID_ENUM);
	}
}
#endif
/** Cut current mesh.
 *  This function cuts the current bound mesh into one or more topological disks.
 *  \ingroup cutting
 */
void GIAPIENTRY giCut()
{
	GICutter *pCutter = &(GIContext_current()->cutter);
	GIMesh *pMesh = pCutter->context->mesh;

	/* error checking */
	if(!pMesh)
	{
		GIContext_error(pCutter->context, GI_INVALID_OPERATION);
		return;
	}
	GIMesh_destroy_cut(pMesh);

	/* select algorithm and cut */
	switch(pCutter->cutter)
	{
	case GI_INITIAL_GIM:
		GICutter_initial_gim(pCutter, pMesh);
		//GICutter_face_clustering(pCutter, pMesh);
		//GICutter_seamless_atlas(pCutter, pMesh);
		break;
	case GI_CATMULL_CLARK_SUBDIVISION:
		GICutter_catmull_clark(pCutter, pMesh);
/*		break;
	case GI_FACE_CLUSTERING:
		GICutter_face_clustering(pCutter, pMesh);
*/	}
}

/** \internal
 *  \brief Cutter constructor.
 *  \param cutter cutter to construct
 *  \param context context to construct in
 *  \ingroup cutting
 */
void GICutter_construct(GICutter *cutter, GIContext *context)
{
	/* initialize state */
	cutter->context = context;
	cutter->cutter = GI_INITIAL_GIM;
	cutter->straighten = GI_TRUE;
	cutter->iterations = 1;
	cutter->orientation_weight = 1.0f;
	cutter->shape_weight = 1.0f;
}

/** \internal
 *  \brief Compute cut from existing parameterization discontinuities.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \return number of created patches (0 on error)
 *  \ingroup cutting
 */
GIint GICutter_from_params(GICutter *cutter, GIMesh *mesh)
{
	GIPatch *pPatch;
	GIFace *pFace, *pFTwin, *pFaces = NULL;
	GIHalfEdge *pHalfEdge, *pHTwin;
	GIPatchFaces *pPFacesList = NULL, *pPFaces;
	GIint *pFacePatchMap;
	GIubyte *pCutEdgeFlags;
	GIuint *pHCounts;
	GIDynamicQueue qFaces;
	GIboolean bSuccess;
	GIint i;

	/* auxiliary datastructures */
	pFacePatchMap = (GIint*)GI_MALLOC_ARRAY(mesh->fcount, sizeof(GIint));
	pCutEdgeFlags = (GIubyte*)GI_CALLOC_ARRAY(mesh->ecount, sizeof(GIubyte));
	GIDynamicQueue_construct(&qFaces);
	for(i=0; i<mesh->fcount; ++i)
		pFacePatchMap[i] = -1;

	/* sort faces into patches by region growing */
	while(mesh->faces)
	{
		pPFaces = (GIPatchFaces*)GI_CALLOC_SINGLE(sizeof(GIPatchFaces));
		pPFaces->next = pPFacesList;
		pPFacesList = pPFaces;
		pPFaces->faces = mesh->faces;
		GIDynamicQueue_enqueue(&qFaces, mesh->faces);
		pFacePatchMap[mesh->faces->id] = mesh->patch_count;
		while(qFaces.size)
		{
			pFace = (GIFace*)GIDynamicQueue_dequeue(&qFaces);
			GI_LIST_REMOVE(mesh->faces, pFace);
			GI_LIST_ADD(pFaces, pFace);
			++pPFaces->fcount;

			/* look at neighbouring faces */
			GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
				pHTwin = pHalfEdge->twin;
				pFTwin = pHTwin->face;
				if(pFTwin && GI_VEC2_EQUAL((GIfloat*)pHalfEdge->pstart, 
					(GIfloat*)pHTwin->next->pstart) && 
					GI_VEC2_EQUAL((GIfloat*)pHalfEdge->next->pstart, 
					(GIfloat*)pHTwin->pstart))
				{
					if(pFacePatchMap[pFTwin->id] == -1)
					{
						GIDynamicQueue_enqueue(&qFaces, pFTwin);
						pFacePatchMap[pHTwin->face->id] = mesh->patch_count;
					}
				}
				else
				{
					if(!pCutEdgeFlags[pHalfEdge->edge->id])
					{
						pCutEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_USED;
						++pHalfEdge->vstart->cut_degree;
						++pHalfEdge->twin->vstart->cut_degree;
					}
					++pPFaces->hcount;
				}
			GI_LIST_NEXT(pFace->hedges, pHalfEdge)
		}
		++mesh->patch_count;
	}
	GIDynamicQueue_destruct(&qFaces);

	/* create patches */
	mesh->patches = (GIPatch*)GI_CALLOC_ARRAY(mesh->patch_count, sizeof(GIPatch));
	pHCounts = (GIuint*)GI_MALLOC_ARRAY(mesh->patch_count, sizeof(GIuint));
	for(i=mesh->patch_count-1; i>=0; --i)
	{
		pPFaces = pPFacesList;
		pPFacesList = pPFacesList->next;
		pPatch = mesh->patches + i;
		pPatch->id = i;
		pPatch->mesh = mesh;
		pPatch->fcount = pPFaces->fcount;
		pPatch->faces = pPFaces->faces;
		GIDynamicQueue_construct(&pPatch->split_paths);
		pPatch->next = mesh->patches + ((i+1)%mesh->patch_count);
		pHCounts[i] = pPFaces->hcount;
		GI_FREE_SINGLE(pPFaces, sizeof(GIPatchFaces));
	}
	mesh->faces = mesh->patches->faces;

	/* create parameter coordinates */
	bSuccess = GICutter_create_params(cutter, 
		mesh, pFacePatchMap, pCutEdgeFlags);
	GI_FREE_ARRAY(pFacePatchMap);
	GI_FREE_ARRAY(pCutEdgeFlags);
	if(!bSuccess)
	{
		GI_FREE_ARRAY(pHCounts);
		GIMesh_destroy_cut(mesh);
		GIContext_error(cutter->context, GI_INVALID_CUT);
		return 0;
	}

	/* check for valid parameterization */
	for(i=0; i<mesh->patch_count; ++i)
	{
		pPatch = mesh->patches + i;
		if(GIPatch_valid_parameterization(pPatch))
		{
			pPatch->parameterized = GI_TRUE;
			pPatch->resolution = UINT_MAX;
			++mesh->param_patches;
		}
	}

	/* compute paths, prevent simgularities and renumerate */
	GICutter_compute_paths(cutter, mesh);
	for(i=0; i<mesh->patch_count; ++i)
	{
		pPatch = mesh->patches + i;
		if(pPatch->hcount != pHCounts[i])
		{
			/* no cylindrical patches */
			GI_FREE_ARRAY(pHCounts);
			GIMesh_destroy_cut(mesh);
			GIContext_error(cutter->context, GI_INVALID_CUT);
			return 0;
		}
		GIPatch_find_corners(pPatch);
		GIPatch_prevent_singularities(pPatch);
		GIPatch_renumerate_params(pPatch);
	}
	GI_FREE_ARRAY(pHCounts);
	mesh->cut_splits = mesh->split_hedges.size;
	mesh->varray_attribs = 0;
	return mesh->patch_count;
}

/** \internal
 *  \brief Gu's initial cut computation.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \return number of created patches (0 on error)
 *  \ingroup cutting
 */
GIint GICutter_initial_gim(GICutter *cutter, GIMesh *mesh)
{
	GIPatch *pPatch;
	GIFace *pFace;
	GIEdge *pEdge, *pEWork;
	GIHalfEdge *pHalfEdge, *pHWork;
	GIVertex *pVertex, *pRoot, *pVWork;
	GIParam *pParam;
	GIHeap qFringe;
	GIDynamicQueue qVertices;
	GIubyte *pEdgeFlags, *pFlag;
	GIdouble vec[3], vec2[3];
	GIdouble dDist, dOldDist, dMinDist;
	GIuint i, uiCount = 0, uiNodes = 0;
	GIubyte ubDegree;
	GIboolean bCenter;

	/* create patch */
	mesh->patches = pPatch = (GIPatch *)GI_CALLOC_ARRAY(sizeof(GIPatch), 1);
	pPatch->id = 0;
	mesh->patch_count = 1;
	pPatch->next = pPatch;
	pPatch->mesh = mesh;
	pPatch->fcount = mesh->fcount;
	pPatch->faces = mesh->faces;
	GIDynamicQueue_construct(&pPatch->split_paths);

	/* initialize datastructures */
	GIHeap_construct(&qFringe, mesh->ecount, lessd, -DBL_MAX, GI_FALSE);
	GIDynamicQueue_construct(&qVertices);
	pEdgeFlags = (GIubyte*)GI_CALLOC_ARRAY(mesh->ecount, sizeof(GIubyte));

	/* seed triangle */
	pFace = mesh->faces;
	GIFace_center(pFace, vec);
	GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
		pEdge = pHalfEdge->edge;
		if(pHalfEdge->twin->face)
		{
			/* no boundary edge, compute geodesic distance and add to queue */
			GIFace_center(pHalfEdge->twin->face, vec2);
			GIHeap_enqueue(&qFringe, pEdge, GIvec3d_dist(vec, vec2));
			pEdgeFlags[pEdge->id] = GI_EDGE_GREY;
			if(pHalfEdge == &pEdge->hedge[1])
				pEdgeFlags[pEdge->id] |= GI_EDGE_FACE0;
		}
		else
		{
			/* boundary edge, definitely on cut -> increase degrees of vertices */
			pEdgeFlags[pEdge->id] = GI_EDGE_USED;
			if(++pHalfEdge->vstart->cut_degree == 1)
				GIDynamicQueue_enqueue(&qVertices, pHalfEdge->vstart);
			if(++pHalfEdge->next->vstart->cut_degree == 1)
				GIDynamicQueue_enqueue(&qVertices, pHalfEdge->next->vstart);
		}
		++uiCount;
	GI_LIST_NEXT(pFace->hedges, pHalfEdge)

	/* process edge priority queue ("remove edges and triangles") */
	GIDebug(printf("remove triangles\n"));
	while(qFringe.count)
	{
		pEWork = (GIEdge*)GIHeap_dequeue(&qFringe, &dDist);
		pFlag = pEdgeFlags + pEWork->id;
		if(*pFlag & GI_EDGE_GREY)
		{
			/* edge processed */
			if(*pFlag & GI_EDGE_FACE0)
				pFace = pEWork->hedge[0].face;
			else
				pFace = pEWork->hedge[1].face;
			*pFlag = GI_EDGE_BLACK;
			bCenter = GI_FALSE;

			/* remove triangle (add edges to queue) */
			GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
				pEdge = pHalfEdge->edge;
				if(pEdge != pEWork)
				{
					/* edge not allready in queue? */
					pFlag = pEdgeFlags + pEdge->id;
					if(*pFlag == GI_EDGE_WHITE)
					{
						/* edge on boundary? */
						if(pHalfEdge->twin->face)
						{
							/* no boundary edge -> compute geodesic distance and add to queue */
							if(!bCenter)
							{
								GIFace_center(pFace, vec);
								bCenter = GI_TRUE;
							}
							GIFace_center(pHalfEdge->twin->face, vec2);
							GIHeap_enqueue(&qFringe, pEdge, dDist+GIvec3d_dist(vec, vec2));
							*pFlag = GI_EDGE_GREY;
							if(pHalfEdge == &pEdge->hedge[1])
								*pFlag |= GI_EDGE_FACE0;
						}
						else
						{
							/* boundary edge, definitely on cut -> increase degrees of vertices */
							*pFlag = GI_EDGE_USED;
							if(++pHalfEdge->vstart->cut_degree == 1)
								GIDynamicQueue_enqueue(&qVertices, pHalfEdge->vstart);
							if(++pHalfEdge->next->vstart->cut_degree == 1)
								GIDynamicQueue_enqueue(&qVertices, pHalfEdge->next->vstart);
						}
						++uiCount;
					}
					else if(*pFlag & GI_EDGE_GREY)
					{
						/* edge perhaps on cut -> increase degrees of vertices */
						*pFlag = GI_EDGE_USEABLE;
						if(++pHalfEdge->vstart->cut_degree == 1)
							GIDynamicQueue_enqueue(&qVertices, pHalfEdge->vstart);
						if(++pHalfEdge->next->vstart->cut_degree == 1)
							GIDynamicQueue_enqueue(&qVertices, pHalfEdge->next->vstart);
					}
				}
			GI_LIST_NEXT(pFace->hedges, pHalfEdge)
		}
	}
	GIHeap_destruct(&qFringe);

	/* mesh not connected */
	if(uiCount != mesh->ecount)
	{
		printf("cutting error:mesh not connected\n");
		GIDebug(printf("not connected: %d / %d\n", uiCount, mesh->ecount));
		GIContext_error(cutter->context, GI_INVALID_CUT);
		GIDynamicQueue_destruct(&qVertices);
		GI_FREE_ARRAY(pEdgeFlags);
		GIMesh_destroy_cut(mesh);
		return 0;
	}

	/* process vertex queue ("remove vertices and edges") */
	GIDebug(printf("remove vertices\n"));
	uiCount = qVertices.size;
	while(qVertices.size)
	{
		pVWork = (GIVertex*)GIDynamicQueue_dequeue(&qVertices);
		if(pVWork->cut_degree == 1)
		{
			/* find cut-edge */
			pVWork->cut_degree = 0;
			pHalfEdge = pVWork->hedge;
			pFlag = pEdgeFlags + pHalfEdge->edge->id;
			while(*pFlag != GI_EDGE_USEABLE)
			{
				pHalfEdge = pHalfEdge->twin->next;
				pFlag = pEdgeFlags + pHalfEdge->edge->id;
			}

			/* decrease degree of cut-adjacent vertex and add to queue if 1 */
			*pFlag = GI_EDGE_BLACK;
			pVertex = pHalfEdge->next->vstart;
			--pVertex->cut_degree;
			if(pVertex->cut_degree == 1)
				GIDynamicQueue_enqueue(&qVertices, pVertex);
			else if(!pVertex->cut_degree)
				--uiCount;
			--uiCount;
		}
	}

	/* create one param per vertex */
	GI_LIST_FOREACH(mesh->vertices, pVertex)
		pParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
		GI_LIST_ADD(pPatch->params, pParam);
		pParam->id = pPatch->pcount++;
		GI_VEC2_SET(pParam->params, 0.0, 0.0);
		pParam->stretch = 0.0;
		pParam->vertex = pVertex;
		pParam->cut_hedge = NULL;
		pVertex->hedge->pstart = pVertex->hedge->face ? pParam : NULL;
		for(pHalfEdge=pVertex->hedge->twin->next; 
			pHalfEdge!=pVertex->hedge; pHalfEdge=pHalfEdge->twin->next)
			pHalfEdge->pstart = pParam;
	GI_LIST_NEXT(mesh->vertices, pVertex)

	/* closed genus-0 mesh? */
	if(!uiCount)
	{
		GIdouble dCos, dCosPrev;

		/* construct cut */
		pVertex = mesh->vertices;
		pHalfEdge = pVertex->hedge;
		pHalfEdge->pstart->cut_hedge = pHalfEdge;
		pVertex->cut_degree = 2;
		pHalfEdge->next->pstart->cut_hedge = pHalfEdge->twin;
		pHalfEdge->next->vstart->cut_degree = 1;
		pPatch->params = pHalfEdge->next->pstart;
		pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_USED;

		/* create new param */
		pParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
		GI_LIST_ADD(pPatch->params, pParam);
		pParam->id = pPatch->pcount++;
		GI_VEC2_SET(pParam->params, 0.0, 0.0);
		pParam->stretch = 0.0;
		pParam->vertex = pVertex;

		/* adjust params and find second cut-edge */
		GI_VEC3_SUB(vec, pHalfEdge->next->vstart->coords, pVertex->coords);
		pHalfEdge = pHalfEdge->twin->next;
		GI_VEC3_SUB(vec2, pHalfEdge->next->vstart->coords, pVertex->coords);
		dCos = GI_VEC3_DOT(vec, vec2)/(GI_VEC3_LENGTH(vec)*GI_VEC3_LENGTH(vec2));
		do
		{
			pHalfEdge->pstart = pParam;
			pHalfEdge = pHalfEdge->twin->next;
			GI_VEC3_SUB(vec2, pHalfEdge->next->vstart->coords, pVertex->coords);
			dCosPrev = dCos;
			dCos = GI_VEC3_DOT(vec, vec2)/(GI_VEC3_LENGTH(vec)*GI_VEC3_LENGTH(vec2));
		}while(dCos < dCosPrev);
		pHalfEdge = pHalfEdge->prev->twin;

		/* extend cut and clean up */
		pParam->cut_hedge = pHalfEdge;
		pHalfEdge->next->pstart->cut_hedge = pHalfEdge->twin;
		pHalfEdge->next->vstart->cut_degree = 1;
		pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_USED;
	}
	else
	{
		/* find vertex with highest degree */
		GIDebug(printf("find root\n"));
		pRoot = mesh->vertices;
		uiCount = 0;
		GI_LIST_FOREACH(mesh->vertices, pVertex)
			if(pVertex->cut_degree == 2)
				++uiCount;
			else if(pVertex->cut_degree)
				++uiNodes;
			if(pVertex->cut_degree > pRoot->cut_degree)
				pRoot = pVertex;
		GI_LIST_NEXT(mesh->vertices, pVertex)

		/* straighten cut paths */
		if(cutter->straighten && pRoot->cut_degree > 2)
		{
			GIVertex *pVStart, *pVEnd, *pVOther, *pVPrev;
			GIHash hNodes;
			GIDynamicQueue qCutPath;
			GIdouble *vend;
			GIuint uiID;

			/* initialize datastructures */
			GIHalfEdge **pSources = (GIHalfEdge**)GI_CALLOC_ARRAY(mesh->vcount, sizeof(GIHalfEdge*));
			GIdouble *pDistances = (GIdouble*)GI_MALLOC_ARRAY(mesh->vcount, sizeof(GIdouble));
			GIdouble *pCutDistances = (GIdouble*)GI_MALLOC_ARRAY(mesh->ecount, sizeof(GIdouble));
			GIHeap_construct(&qFringe, mesh->ecount, lessd, -DBL_MAX, GI_TRUE);
			GIHash_construct(&hNodes, uiNodes, 0.0f, sizeof(GIuint), hash_uint, compare_uint, copy_uint);
			GIDynamicQueue_construct(&qCutPath);

			/* process cut nodes */
			GIDebug(printf("straighten cut:\n"));
			GIDynamicQueue_enqueue(&qVertices, pRoot);
			while(qVertices.size)
			{
				/* dequeue vertex and check if allready processed */
				pVWork = (GIVertex *)GIDynamicQueue_dequeue(&qVertices);
				if(!GIHash_insert(&hNodes, &pVWork->id, NULL))
					continue;

				/* check if on boundary and walk to "leftmost" boundary */
				pDistances[pVWork->id] = 0.0;
				ubDegree = pVWork->cut_degree;
				pHWork = pVWork->hedge;
				if(!pHWork->face)
				{
					--ubDegree;
					pHWork = pHWork->twin->next;
				}

				/* process cut paths at this cut node */
				while(ubDegree)
				{
					pEdge = pHWork->edge;
					if(pEdgeFlags[pEdge->id] == GI_EDGE_USED)
					{
						--ubDegree;
						if(!pHWork->twin->face)
						{
							pVertex = pHWork->next->vstart;
							while(pVertex->cut_degree == 2)
								pVertex = pVertex->hedge->prev->vstart;
							if(pVertex != pVWork)
								GIDynamicQueue_enqueue(&qVertices, pVertex);
							break;
						}
					}
					else if(pEdgeFlags[pEdge->id] == GI_EDGE_USEABLE)
					{
						/* init datastructures */
//						for(i=0; i<mesh->vcount; ++i)
//							pWeightScale[i] = GI_REAL_MAX;

						/* delete old cut path and find end vertex */
/*						GIDebug(printf("   delete ... "));
						--ubDegree;
						pEdgeFlags[pEdge->id] = GI_EDGE_USED;
						pVStart = pHWork->next->vstart;
						if(pVStart->cut_degree != 2)
						{
							GIDynamicQueue_enqueue(&qVertices, pVStart);
							continue;
						}
						pVPrev = pVStart;
						pHalfEdge = pVStart->hedge;
						while(pEdgeFlags[pHalfEdge->edge->id] != GI_EDGE_USEABLE)
							pHalfEdge = pHalfEdge->twin->next;
						pVertex = pHalfEdge->next->vstart;
						while(pVertex->cut_degree == 2)
						{
							pVertex->cut_degree = 0;
							pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_BLACK;
							pHalfEdge = pVertex->hedge;
							while(pEdgeFlags[pHalfEdge->edge->id] != GI_EDGE_USEABLE)
								pHalfEdge = pHalfEdge->twin->next;
							pVPrev = pVertex;
							pVertex = pHalfEdge->next->vstart;
						}
						pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_USED;
						GIDynamicQueue_enqueue(&qVertices, pVertex);
						pVEnd = pVPrev;
						vend = pVEnd->coords;
*/
						GIDebug(printf("   delete ... "));
						--ubDegree;
						pEdgeFlags[pEdge->id] = GI_EDGE_USED;
						pVStart = pHWork->next->vstart;
						if(pVStart->cut_degree != 2)
						{
							GIDynamicQueue_enqueue(&qVertices, pVStart);
							pHWork = pHWork->twin->next;
							continue;
						}
						GIHeap_clear(&qFringe);
						for(i=0; i<mesh->ecount; ++i)
							pCutDistances[i] = DBL_MAX;
						pVPrev = pVStart;
						pHalfEdge = pVStart->hedge;
						while(pEdgeFlags[pHalfEdge->edge->id] != GI_EDGE_USEABLE)
							pHalfEdge = pHalfEdge->twin->next;
						pVertex = pHalfEdge->next->vstart;
						while(pVertex->cut_degree == 2)
						{
							pVertex->cut_degree = 0;
							pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_BLACK;
							pCutDistances[pHalfEdge->edge->id] = 0.0;
							GIHeap_enqueue(&qFringe, pHalfEdge->edge, 0.0);
							GIDynamicQueue_enqueue(&qCutPath, pHalfEdge->edge);
							pHalfEdge = pVertex->hedge;
							while(pEdgeFlags[pHalfEdge->edge->id] != GI_EDGE_USEABLE)
								pHalfEdge = pHalfEdge->twin->next;
							pVPrev = pVertex;
							pVertex = pHalfEdge->next->vstart;
						}
						pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_USED;
						pCutDistances[pHalfEdge->edge->id] = 0.0;
						GIHeap_enqueue(&qFringe, pHalfEdge->edge, 0.0);
						GIDynamicQueue_enqueue(&qCutPath, pHalfEdge->edge);
						pCutDistances[pHWork->edge->id] = 0.0;
						GIHeap_enqueue(&qFringe, pHWork->edge, 0.0);
						GIDynamicQueue_enqueue(&qCutPath, pHWork->edge);
						GIDynamicQueue_enqueue(&qVertices, pVertex);
						pVEnd = pVPrev;
						vend = pVEnd->coords;

						/* compute distances of edges to old cut */
						GIDebug(printf("dijkstra ... "));
						dMinDist = DBL_MAX;
						while(qFringe.count)
						{
							pEdge = (GIEdge *)GIHeap_dequeue(&qFringe, &dOldDist);
							GI_VEC3_ADD(vec, pEdge->hedge[0].vstart->coords, 
								pEdge->hedge[1].vstart->coords);
							GI_VEC3_SCALE(vec, vec, 0.5);
							pHalfEdge = &pEdge->hedge[0];
							for(i=0; i<4; ++i)
							{
								pHalfEdge = pHalfEdge->next;
								pEWork = pHalfEdge->edge;
								if(!(pEdgeFlags[pEWork->id] & GI_EDGE_USEABLE))
								{
									GI_VEC3_ADD(vec2, pEWork->hedge[0].vstart->coords, 
										pEWork->hedge[1].vstart->coords);
									GI_VEC3_SCALE(vec2, vec2, 0.5);
									dDist = dOldDist + GIvec3d_dist(vec, vec2);
									if(dDist < pCutDistances[pEWork->id])
									{
										pCutDistances[pEWork->id] = dDist;
										GIHeap_enqueue(&qFringe, pEWork, dDist);
										if(dDist < dMinDist)
											dMinDist = dDist;
									}
								}
								if(i == 1)
									pHalfEdge = &pEdge->hedge[1];
							}
						}
//						dMinDist *= 0.5;
						while(qCutPath.size)
							pCutDistances[((GIEdge*)GIDynamicQueue_dequeue(
								&qCutPath))->id] = dMinDist;

						/* find shortest path (Dijkstra) */
						GIDebug(printf("dijkstra ... "));
						GIHeap_clear(&qFringe);
						pDistances[pVStart->id] = 0.0;
						memset(pSources, 0, mesh->vcount*sizeof(GIHalfEdge*));
						pVertex = pVStart;
						while(pVertex != pVEnd)
						{
							pHalfEdge = pVertex->hedge;
							do
							{
								/* path may not cross another cut path */
								pVOther = pHalfEdge->next->vstart;
								if(!pVOther->cut_degree)
								{
									uiID = pVOther->id;
									dDist = pDistances[pVertex->id] + pHalfEdge->edge->length /* 
										pCutDistances[pHalfEdge->edge->id]*/;
									if(!pSources[uiID] || dDist < pDistances[uiID])
									{
										pSources[uiID] = pHalfEdge;
										pDistances[uiID] = dDist;
										GIHeap_enqueue(&qFringe, pVOther, 
											dDist/*+GIvec3d_dist(pVOther->coords, vend)*/);
									}
								}
								pHalfEdge = pHalfEdge->twin->next;
							}while(pHalfEdge != pVertex->hedge);
							pVertex = (GIVertex *)GIHeap_dequeue(&qFringe, NULL);
						}

						/* retrace shortest path and assemble cut path */
						GIDebug(printf("assemble\n"));
						pHalfEdge = pSources[pVEnd->id];
						if ( pHalfEdge == NULL )
						{
							break;
						}
						pVertex = pHalfEdge->vstart;
						while(pVertex != pVStart)
						{
							pVertex->cut_degree = 2;
							pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_USED;
							pHalfEdge = pSources[pVertex->id];
							pVertex = pHalfEdge->vstart;
						}
						pVEnd->cut_degree = 2;
						pEdgeFlags[pHalfEdge->edge->id] = GI_EDGE_USED;
					}
					pHWork = pHWork->twin->next;
				}
			}

			/* clean up */
			GI_FREE_ARRAY(pSources);
			GI_FREE_ARRAY(pDistances);
			GIHeap_destruct(&qFringe);
			GIHash_destruct(&hNodes, 0);
		}
		ubDegree = pRoot->cut_degree;
		if(!pRoot->hedge->face)
			--ubDegree;

		/* construct cut */
		GIDebug(printf("construct cut: %d\n", ubDegree));
		pVertex = pRoot;
		pHalfEdge = pVertex->hedge->twin;
		while(!(pEdgeFlags[pHalfEdge->edge->id] & GI_EDGE_USEABLE))
			pHalfEdge = pHalfEdge->twin->prev;
		pHalfEdge = pHalfEdge->next;
		pPatch->params = pHalfEdge->pstart;
		while(pVertex != pRoot || ubDegree)
		{
			if(pVertex == pRoot)
				--ubDegree;

			/* visiting vertex a second time -> new param */
			pParam = pHalfEdge->pstart;
			if(pParam->cut_hedge)
			{
				pParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
				GI_LIST_ADD(pPatch->params, pParam);
				pParam->id = pPatch->pcount++;
				GI_VEC2_SET(pParam->params, 0.0, 0.0);
				pParam->stretch = 0.0;
				pParam->vertex = pVertex;
				pHalfEdge->pstart = pParam;
			}

			/* change param of relevant half edges and prevent degenerate triangles */
			while(!(pEdgeFlags[pHalfEdge->edge->id] & GI_EDGE_USEABLE))
			{
				pHalfEdge = pHalfEdge->twin->next;
				pHalfEdge->pstart = pParam;
			}

			/* extend cut */
			pParam->cut_hedge = pHalfEdge;
			pHalfEdge = pHalfEdge->next;
			pVertex = pHalfEdge->vstart;
		}
	}

	/* clean up and compute cut paths */
	GI_FREE_ARRAY(pEdgeFlags);
	GIDynamicQueue_destruct(&qVertices);
	GICutter_compute_paths(cutter, mesh);
	GIPatch_find_corners(pPatch);
	GIPatch_prevent_singularities(pPatch);
	GIPatch_renumerate_params(pPatch);
	mesh->cut_splits = mesh->split_hedges.size;
	mesh->varray_attribs = 0;
	return 1;
}

/** \internal
 *  \brief Garland's hierarchical face clustering.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \return number of created patches (0 on error)
 *  \ingroup cutting
 */
GIint GICutter_face_clustering(GICutter *cutter, GIMesh *mesh)
{
	static GIfloat zero[2] = { 0.0f, 0.0f };
	GIFaceCluster *pClusters = NULL, *pCluster;
	GIPatch *pPatch;
	GIFace *pFace;
	GIEdge *pEdge;
	GIVertex *pVRoot = NULL;
	GIint *pFacePatchMap;
	GIubyte *pCutEdgeFlags;
	GIboolean bSuccess;
	GIuint i = 0;

	/* cluster faces */
	mesh->patch_count = GIMultiresolution_face_clustering(mesh, pClusters, 
		cutter->orientation_weight, cutter->shape_weight, 
		DBL_MAX, 128, UINT_MAX);
	mesh->patches = (GIPatch*)GI_CALLOC_ARRAY(mesh->patch_count, sizeof(GIPatch));
	pFacePatchMap = (GIint*)GI_MALLOC_ARRAY(mesh->fcount, sizeof(GIint));
	pCutEdgeFlags = (GIubyte*)GI_CALLOC_ARRAY(mesh->ecount, sizeof(GIubyte));

	printf("cluster faces:%d\n", mesh->patch_count);

	printf("create patches start\n");
	/* compute patch membership and create patches */
	GI_LIST_FOREACH(pClusters, pCluster)
		pPatch = mesh->patches + i;
		pPatch->id = i;
		pPatch->mesh = mesh;
		GIDynamicQueue_construct(&pPatch->split_paths);
		pPatch->next = mesh->patches + ((i+1)%mesh->patch_count);
		if(cutter->straighten)
		{
			pPatch->fcount = pCluster->fcount;
			pPatch->faces = pCluster->faces;
		}
		GI_LIST_FOREACH(pCluster->faces, pFace)
			pFacePatchMap[pFace->id] = i;
		GI_LIST_NEXT(pCluster->faces, pFace)
		++i;
	GI_LIST_NEXT(pClusters, pCluster)
	GI_LIST_CLEAR(pClusters, sizeof(GIFaceCluster))
	
	printf("create patches end\n");

	printf("cut edges start\n");
	/* mark cut edges */
	GI_LIST_FOREACH(mesh->edges, pEdge)
		if(!pEdge->hedge[0].face || !pEdge->hedge[1].face || 
			pFacePatchMap[pEdge->hedge[0].face->id] != 
			pFacePatchMap[pEdge->hedge[1].face->id])
		{
			pCutEdgeFlags[pEdge->id] = GI_EDGE_USEABLE;
			if(++pEdge->hedge[0].vstart->cut_degree == 3)
				pVRoot = pEdge->hedge[0].vstart;
			if(++pEdge->hedge[1].vstart->cut_degree == 3)
				pVRoot = pEdge->hedge[1].vstart;
		}
		pEdge->hedge[0].pstart = pEdge->hedge[1].pstart = (GIParam*)zero;
	GI_LIST_NEXT(mesh->edges, pEdge)
	printf("cut edges end\n");

	/* straighten cut paths if neccessary */
	if(cutter->straighten)
	{
		GICutter_straighten_paths(cutter, mesh, pVRoot, pCutEdgeFlags);
		GICutter_sort_faces(cutter, mesh, pFacePatchMap, pCutEdgeFlags);
	}
	mesh->faces = mesh->patches->faces;

	/* create parameter coordinates */
	bSuccess = GICutter_create_params(cutter, 
		mesh, pFacePatchMap, pCutEdgeFlags);
	GI_FREE_ARRAY(pFacePatchMap);
	GI_FREE_ARRAY(pCutEdgeFlags);
	if(!bSuccess)
	{
		GIMesh_destroy_cut(mesh);
		GIContext_error(cutter->context, GI_INVALID_CUT);
		return 0;
	}

	/* compute paths, prevent simgularities and renumerate */
	GICutter_compute_paths(cutter, mesh);
	for(i=0; i<mesh->patch_count; ++i)
	{
		pPatch = mesh->patches + i;
		GIPatch_find_corners(pPatch);
		GIPatch_prevent_singularities(pPatch);
		GIPatch_renumerate_params(pPatch);
	}
	mesh->cut_splits = mesh->split_hedges.size;
	mesh->varray_attribs = 0;
	return mesh->patch_count;
}

/** \internal
 *  \brief Purnomo's seamless texture atlas patchification.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \return number of created patches (0 on error)
 *  \ingroup cutting
 */
GIint GICutter_seamless_atlas(GICutter *cutter, GIMesh *mesh)
{
	static GIfloat zero[2] = { 0.0f, 0.0f };
	GIFaceCluster *pClusters = NULL, *pCluster;
	GIPatch *pPatch;
	GIFace *pFace;
	GIEdge *pEdge;
	GIHalfEdge *pHalfEdge, *pHStart;
	GIVertex *pVRoot = NULL, *pVertex;
	GIint *pFacePatchMap;
	GIubyte *pCutEdgeFlags;
	GIClusterPatch *pPatches;
	GIClusterPatch *pClusterPatch;
	GIHash hMedians;
	GIboolean bSuccess;
	GIuint i = 0, j, uiClusters, uiCorners = 0;
	GIuint uiMaxFaces = mesh->fcount, uiMaxEdges = mesh->ecount;

	/* cluster faces */
	uiClusters = GIMultiresolution_face_clustering(mesh, pClusters, 
		cutter->orientation_weight, cutter->shape_weight, 
		DBL_MAX, 128, UINT_MAX);
	pFacePatchMap = (GIint*)GI_MALLOC_ARRAY(mesh->fcount, sizeof(GIint));
	pCutEdgeFlags = (GIubyte*)GI_CALLOC_ARRAY(mesh->ecount, sizeof(GIubyte));
	pPatches = (GIClusterPatch*)GI_MALLOC_ARRAY(uiClusters, sizeof(GIClusterPatch));

	/* compute patch membership and create patches */
	GI_LIST_FOREACH(pClusters, pCluster)
		GI_LIST_FOREACH(pCluster->faces, pFace)
			pFacePatchMap[pFace->id] = i;
		GI_LIST_NEXT(pCluster->faces, pFace)
		++i;
	GI_LIST_NEXT(pClusters, pCluster)
	GI_LIST_CLEAR(pClusters, sizeof(GIFaceCluster))

	/* mark cut edges */
	GI_LIST_FOREACH(mesh->edges, pEdge)
		if(!pEdge->hedge[0].face || !pEdge->hedge[1].face || 
			pFacePatchMap[pEdge->hedge[0].face->id] != 
			pFacePatchMap[pEdge->hedge[1].face->id])
		{
			pCutEdgeFlags[pEdge->id] = GI_EDGE_USEABLE;
			for(i=0; i<2; ++i)
			{
				/* remember start vertex for every patch */
				pVertex = pEdge->hedge[i].vstart;
				if(++pVertex->cut_degree == 3)
				{
					pVRoot = pVertex;
					pHalfEdge = pVertex->hedge->face ? 
						pVertex->hedge : pVertex->hedge->twin->next;
					do
					{
						if(pCutEdgeFlags[pHalfEdge->edge->id])
						{
							pClusterPatch = pPatches + 
								pFacePatchMap[pHalfEdge->face->id];
							pClusterPatch->root = pVertex;
							++pClusterPatch->bcount;
						}
						pHalfEdge = pHalfEdge->twin->next;
					}while(pHalfEdge != pVertex->hedge);
					++uiCorners;
				}
			}
		}
//		pEdge->hedge[0].pstart = pEdge->hedge[1].pstart = (GIParam*)zero;
	GI_LIST_NEXT(mesh->edges, pEdge)

	/* straighten cut paths */
	GICutter_straighten_paths(cutter, mesh, pVRoot, pCutEdgeFlags);

	/* find cluster boundaries and medians of boundaries */
	GIHash_construct(&hMedians, uiCorners+uiClusters, 0.0f, 
		sizeof(GIUIntPair), hash_uintpair, compare_uintpair, copy_uintpair);
	for(i=0; i<uiClusters; ++i)
	{
		pClusterPatch = pPatches + i;
		pClusterPatch->boundaries = (GIPatchBoundary*)GI_MALLOC_ARRAY(
			pClusterPatch->bcount, sizeof(GIPatchBoundary));
		pHStart = pClusterPatch->root->hedge;
		if(!pHStart->face)
			pHStart = pHStart->twin->next;
//		while(!pCutEdgeFlags
		for(j=0; j<pClusterPatch->bcount; ++j)
		{
			GIUIntPair pair;
			GIdouble dLength = 0.0;

			/* traverse boundary */
			pClusterPatch->boundaries[j].hedge = pHStart;
			do
			{
				dLength += pHStart->edge->length;
				pHStart = pHStart->next;
				while(!pCutEdgeFlags[pHStart->edge->id])
					pHStart = pHStart->twin->next;
			}while(pHStart->vstart->cut_degree == 2);
			pair.first = pClusterPatch->boundaries[j].hedge->vstart->id;
			pair.second = pHStart->vstart->id;

			/* set median */
			pVertex = (GIVertex*)GIHash_find(&hMedians, &pair);
			if(!pVertex)
			{
				/* find median */
				dLength *= 0.5;
				pHalfEdge = pClusterPatch->boundaries[j].hedge;
				dLength -= pHalfEdge->edge->length;
				while(dLength > 0.0)
				{
					pHalfEdge = pHalfEdge->next;
					while(!pCutEdgeFlags[pHalfEdge->edge->id])
						pHalfEdge = pHalfEdge->twin->next;
					dLength -= pHalfEdge->edge->length;
				}

				/* select median and split if neccessary */
				if(pHalfEdge->vstart->cut_degree != 2 && 
					pHalfEdge->next->vstart->cut_degree != 2)
				{
					GIHalfEdge_split(pHalfEdge, NULL, (GIPatch*)mesh, 0.5, NULL);
					if(mesh->ecount > uiMaxEdges)
					{
						uiMaxEdges = (GIint)sqrt((GIdouble)uiMaxEdges) + 1;
						pCutEdgeFlags = (GIubyte*)GI_REALLOC_ARRAY(
							pCutEdgeFlags, uiMaxEdges, sizeof(GIubyte));
					}
					pCutEdgeFlags[pHalfEdge->next->twin->next->edge->id] = GI_EDGE_USED;
				}
				else if(pHalfEdge->vstart->cut_degree != 2 || 
					fabs(dLength)<fabs(dLength+pHalfEdge->edge->length))
					pVertex = pHalfEdge->next->vstart;
				else
					pVertex = pHalfEdge->vstart;
				GIHash_insert(&hMedians, &pair, pVertex);
			}
			pClusterPatch->boundaries[j].median = pVertex;
		}
		mesh->patch_count += pClusterPatch->bcount;
	}
	GIHash_destruct(&hMedians, 0);

	/* split clusters */
	for(i=0; i<uiClusters; ++i)
	{
	}
	GI_FREE_ARRAY(pPatches);

	/* create patches and sort faces into patches */
	mesh->patches = (GIPatch*)GI_CALLOC_ARRAY(mesh->patch_count, sizeof(GIPatch));
	for(i=0; i<mesh->patch_count; ++i)
	{
		pPatch = mesh->patches + i;
		pPatch->id = i;
		pPatch->mesh = mesh;
		GIDynamicQueue_construct(&pPatch->split_paths);
		pPatch->next = mesh->patches + ((i+1)%mesh->patch_count);
	}
	if(mesh->fcount > uiMaxFaces)
	{
		pFacePatchMap = (GIint*)GI_REALLOC_ARRAY(
			pFacePatchMap, mesh->fcount, sizeof(GIint));
	}
	GICutter_sort_faces(cutter, mesh, pFacePatchMap, pCutEdgeFlags);

	/* create parameter coordinates */
	mesh->pre_cut_splits = mesh->split_hedges.size;
	bSuccess = GICutter_create_params(cutter, 
		mesh, pFacePatchMap, pCutEdgeFlags);
	GI_FREE_ARRAY(pFacePatchMap);
	GI_FREE_ARRAY(pCutEdgeFlags);
	if(!bSuccess)
	{
		GIMesh_destroy_cut(mesh);
		GIContext_error(cutter->context, GI_INVALID_CUT);
		return 0;
	}

	/* compute paths, prevent simgularities and renumerate */
	GICutter_compute_paths(cutter, mesh);
	for(i=0; i<mesh->patch_count; ++i)
	{
		pPatch = mesh->patches + i;
		GIPatch_find_corners(pPatch);
		GIPatch_prevent_singularities(pPatch);
		GIPatch_renumerate_params(pPatch);
	}
	mesh->cut_splits = mesh->split_hedges.size;
	mesh->varray_attribs = 0;
	return mesh->patch_count;
}

/** \internal
 *  \brief Catmull-Clark subdivision.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \return number of created patches (0 on error)
 *  \ingroup cutting
 */
GIint GICutter_catmull_clark(GICutter *cutter, GIMesh *mesh)
{
	static const GIdouble corners[8] = { 
		0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
	static const GIdouble d5_3 = 5.0 / 3.0;
	static const GIdouble d5_6 = 5.0 / 6.0;
	static const GIdouble d1_12 = 1.0 / 12.0;
	static const GIdouble d1_36 = 1.0 / 36.0;
	static const GIdouble d4_36 = 1.0 / 9.0;
	static const GIdouble d16_36 = 4.0 / 9.0;
	static const GIdouble d1_6 = 1.0 / 6.0;
	GIPatch *pPatch, *pPatch2;
	GIFace *pFaces = NULL;
	GIEdge *pEdge;
	GIHalfEdge *pHalfEdge, *pHStart;
	GIVertex *pVertex;
	GIParam *pParam;
	GIEdgeAngle *pEdges;
	GIHalfEdge **pFaceHedges;
	GIint *pFacePatchMap;
	GIboolean *pEdgeFlags;
	GIHalfEdge *pBoundary[4];
	GIHalfEdge *pHBStart[3];
	GIdouble *pNewCoords, *vec, *vec2;
	GIfloat *aa, *a, *a0, *a1;
	GIdouble v1[3], v2[3];
	GIuint i, j, k, l, m, uiVCount = mesh->vcount, uiECount = mesh->ecount, 
		uiFCount = mesh->fcount, uiVECount = uiVCount + uiECount, 
		uiResolution = (1<<(cutter->iterations-1)) + 1;

	/* set cut degrees of original vertices */
	vec = pNewCoords = (GIdouble*)GI_MALLOC_ARRAY(3*uiVECount, sizeof(GIdouble));
	GI_LIST_FOREACH(mesh->vertices, pVertex)
		pHalfEdge = pVertex->hedge;
		do
		{
			++pVertex->cut_degree;
			pHalfEdge = pHalfEdge->twin->next;
		}while(pHalfEdge != pVertex->hedge);

		/* compute new vertex position */
		if(pHalfEdge->face)
		{
			/* interior vertex */
			GIdouble dInvN = 1.0 / (GIdouble)pVertex->cut_degree;
			GIdouble d5_3N2 = d5_3 * dInvN * dInvN;
			GIdouble dN_53_N = ((GIdouble)pVertex->cut_degree-d5_3) * dInvN;
			pHalfEdge = pHalfEdge->twin;
			GI_VEC3_COPY(vec, pHalfEdge->vstart->coords);
			for(pHalfEdge=pHalfEdge->twin->prev; pHalfEdge->twin!=pVertex->hedge; 
				pHalfEdge=pHalfEdge->twin->prev)
			{
				GI_VEC3_ADD(vec, vec, pHalfEdge->vstart->coords);
			}
			GI_VEC3_SCALE(vec, vec, d5_3N2);
			GI_VEC3_ADD_SCALED(vec, vec, pVertex->coords, dN_53_N);
		}
		else
		{
			/* boundary vertex */
			GI_VEC3_ADD(vec, pVertex->hedge->next->vstart->coords, 
				pVertex->hedge->prev->vstart->coords);
			GI_VEC3_ADD_SCALED(vec, vec, pVertex->coords, 6.0);
			GI_VEC3_SCALE(vec, vec, 0.125);
		}
		vec += 3;
	GI_LIST_NEXT(mesh->vertices, pVertex)

	/* split longest edge with largest angle first */
	pEdges = (GIEdgeAngle*)GI_MALLOC_ARRAY(uiECount, sizeof(GIEdgeAngle));
	for(i=0,pEdge=mesh->edges; i<uiECount; ++i,pEdge=pEdge->next)
	{
		/* compute angle opposite to edge */
		pEdges[i].edge = pEdge;
		pEdges[i].angle = 0.0;
		for(j=0; j<2; ++j)
		{
			pHalfEdge = &pEdge->hedge[j];
			if(pHalfEdge->face)
			{
				GI_VEC3_SUB(v1, pHalfEdge->vstart->coords, 
					pHalfEdge->prev->vstart->coords);
				GI_VEC3_SUB(v2, pHalfEdge->next->vstart->coords, 
					pHalfEdge->prev->vstart->coords);
				pEdges[i].angle += GI_VEC3_DOT(v1, v2) / 
					(pHalfEdge->next->edge->length*pHalfEdge->prev->edge->length);
			}
		}
		if(pEdge->hedge[1].face)
			pEdges[i].angle *= 0.5;
	}
	qsort(pEdges, uiECount, sizeof(GIEdgeAngle), compare);

	/* compute position of edge vertices */
	for(i=0,vec=pNewCoords+3*uiVCount; i<uiECount; ++i,vec+=3)
	{
		pHalfEdge = &pEdges[i].edge->hedge[1];
		if(pHalfEdge->face)
		{
			GI_VEC3_ADD(vec, pHalfEdge->prev->vstart->coords, 
				pHalfEdge->twin->prev->vstart->coords);
		}
	}

	/* split original edges */
	pFaceHedges = (GIHalfEdge**)GI_CALLOC_ARRAY(uiFCount, sizeof(GIHalfEdge*));
	for(i=0; i<uiECount; ++i)
	{
		pEdge = pEdges[i].edge;
		GIHalfEdge_split(&pEdge->hedge[0], NULL, (GIPatch*)mesh, 0.5, NULL);
		pEdge->hedge[0].next->vstart->cut_degree = 3;
		if(pEdge->hedge[0].face->id < uiFCount && 
			!pFaceHedges[pEdge->hedge[0].face->id])
			pFaceHedges[pEdge->hedge[0].face->id] = pEdge->hedge[0].next;
		if(pEdge->hedge[1].face)
		{
			++pEdge->hedge[1].vstart->cut_degree;
			if(pEdge->hedge[1].face->id < uiFCount && 
				!pFaceHedges[pEdge->hedge[1].face->id])
				pFaceHedges[pEdge->hedge[1].face->id] = 
					pEdge->hedge[1].prev->twin;
		}
	}
	GI_FREE_ARRAY(pEdges);

	/* create auxiliary datastructures if neccessary */
	if(cutter->iterations > 1)
	{
		GIuint f = 6 * uiFCount, e = (uiECount<<1) + f;
		for(k=1; k<cutter->iterations; ++k)
		{
			e = (e<<1) + 6*f;
			f <<= 2;
		}
		pFacePatchMap = (GIint*)GI_MALLOC_ARRAY(f, sizeof(GIint));
		pEdgeFlags = (GIboolean*)GI_CALLOC_ARRAY(e, sizeof(GIboolean));
	}

	/* split original faces and create patches */
	mesh->patch_count = mesh->param_patches = 3 * uiFCount;
	mesh->patches = (GIPatch*)GI_CALLOC_ARRAY(mesh->patch_count, sizeof(GIPatch));
	for(i=0,j=0; i<uiFCount; ++i)
	{
		/* split center halfedge */
		pHalfEdge = pFaceHedges[i];
		GIHalfEdge_split(pHalfEdge, NULL, (GIPatch*)mesh, 0.333333, NULL);
		pHalfEdge->next->vstart->cut_degree = 3;
		pHBStart[0] = pHalfEdge;
		pHBStart[1] = pHalfEdge->twin->prev;
		pHBStart[2] = pHalfEdge->next->twin->prev;

		/* create patches */
		for(k=0; k<3; ++k,++j)
		{
			/* find boundary halfedges */
			pBoundary[0] = pHBStart[k];
			pBoundary[1] = pBoundary[0]->next;
			pBoundary[2] = pBoundary[1]->next->twin->next;
			pBoundary[3] = pBoundary[2]->next;
			for(l=1,m=0; l<4; ++l)
			{
				if(!pBoundary[l-1]->twin->face && !pBoundary[l]->twin->face)
				{
					m = l;
					break;
				}
			}

			/* create patch */
			pPatch = mesh->patches + j;
			pPatch->id = j;
			pPatch->mesh = mesh;
			GIDynamicQueue_construct(&pPatch->split_paths);
			pPatch->next = mesh->patches + ((j+1)%mesh->patch_count);
			pPatch->fcount = 2;
			pPatch->faces = pBoundary[0]->face;
			GI_LIST_REMOVE(mesh->faces, pBoundary[0]->face);
			GI_LIST_ADD(pFaces, pBoundary[0]->face);
			GI_LIST_REMOVE(mesh->faces, pBoundary[2]->face);
			GI_LIST_ADD(pFaces, pBoundary[2]->face);

			/* create params */
			for(l=0; l<4; ++l,m=(m+1)&3)
			{
				pPatch->corners[l] = pBoundary[m]->pstart = pParam = 
					(GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
				GI_LIST_ADD(pPatch->params, pParam);
				pParam->id = l;
				GI_VEC2_COPY(pParam->params, corners+(l<<1));
				pParam->stretch = 0.0;
				pParam->vertex = pBoundary[m]->vstart;
				pParam->cut_hedge = pBoundary[m];
			}
			pBoundary[1]->next->pstart = pBoundary[2]->pstart;
			pBoundary[3]->next->pstart = pBoundary[0]->pstart;
			pPatch->pcount = 4;
			pPatch->fixed_corners = GI_TRUE;
			pPatch->parameterized = GI_TRUE;
			pPatch->resolution = uiResolution;

			/* remember face splitting edges */
			if(cutter->iterations > 1)
			{
				pFacePatchMap[pBoundary[1]->next->face->id] = j;
				pFacePatchMap[pBoundary[3]->next->face->id] = j;
				pEdgeFlags[pBoundary[1]->next->edge->id] = GI_TRUE;
			}
		}
	}
	mesh->faces = mesh->patches->faces;
	mesh->pre_cut_splits = mesh->split_hedges.size;
	GI_FREE_ARRAY(pFaceHedges);

	/* save old positions and set new ones */
	mesh->old_coords = (GIdouble*)GI_MALLOC_ARRAY(3*uiVCount, sizeof(GIdouble));
	for(i=0,vec=pNewCoords,vec2=mesh->old_coords,pVertex=mesh->vertices; 
		i<uiVCount; ++i,vec+=3,vec2+=3,pVertex=pVertex->next)
	{
		GI_VEC3_COPY(vec2, pVertex->coords);
		GI_VEC3_COPY(pVertex->coords, vec);
	}
	for(; i<uiVECount; ++i,vec+=3,pVertex=pVertex->next)
	{
		if(pVertex->hedge->face)
		{
			vec2 = pVertex->coords;
			GI_VEC3_SCALE(vec2, vec2, d5_6);
			GI_VEC3_ADD_SCALED(vec2, vec2, vec, d1_12);
		}
	}
	GI_FREE_ARRAY(pNewCoords);

	/* further subdivision (now only quads) */
	if(cutter->iterations > 1)
	{
		for(k=1; k<cutter->iterations; ++k)
		{
			uiECount = mesh->ecount;
			uiVCount = mesh->vcount;

			/* split faces into four triangles and compute face points */
			for(i=0,pEdge=mesh->edges; i<uiECount; ++i,pEdge=pEdge->next)
			{
				if(pEdgeFlags[i])
				{
					pHalfEdge = &pEdge->hedge[0];
					pPatch = mesh->patches + pFacePatchMap[pHalfEdge->face->id];
					GI_VEC3_ADD(v1, pHalfEdge->prev->vstart->coords, 
						pHalfEdge->twin->prev->vstart->coords);
					GIHalfEdge_split(pHalfEdge, pPatch, pPatch, 0.5, NULL);
					vec = pHalfEdge->next->vstart->coords;
					GI_VEC3_SCALE(vec, vec, 0.5);
					GI_VEC3_ADD_SCALED(vec, vec, v1, 0.25);
					pEdgeFlags[pHalfEdge->next->edge->id] = GI_TRUE;
					pEdgeFlags[pHalfEdge->next->twin->next->edge->id] = GI_TRUE;
					pEdgeFlags[pHalfEdge->twin->prev->edge->id] = GI_TRUE;
					pFacePatchMap[pHalfEdge->next->twin->face->id] = pPatch->id;
					pFacePatchMap[pHalfEdge->twin->prev->twin->face->id] = pPatch->id;

					/* interpolate attributes (diagonal centroid is not face centroid) */
					for(j=0; j<GI_ATTRIB_COUNT; ++j)
					{
						if(mesh->aoffset[j] && mesh->asemantic[j] == GI_NONE)
						{
							for(l=0,a=aa=(GIfloat*)((GIubyte*)pHalfEdge
								->next->astart+mesh->aoffset[j]),
								a0=(GIfloat*)((GIubyte*)pHalfEdge
								->prev->astart+mesh->aoffset[j]),
								a1=(GIfloat*)((GIubyte*)pHalfEdge
								->twin->prev->astart+mesh->aoffset[j]); 
								l<mesh->asize[j]; ++l,++a,++a0,++a1)
									*a = 0.5**a + 0.25*(*a0+*a1);
							if(mesh->anorm[j])
							{
								if(mesh->asize[j] >= 3) GIvec3f_normalize(aa);
								else if(mesh->asize[j] == 2) GIvec2f_normalize(aa);
								else *aa = GI_SIGN(*aa);
							}
						}
					}
				}
			}

			/* split non-face edges and compute edge points */
			for(i=0,pEdge=mesh->edges; i<uiECount; ++i,pEdge=pEdge->next)
			{
				if(!pEdgeFlags[i])
				{
					pHalfEdge = &pEdge->hedge[0];
					pPatch = mesh->patches + pFacePatchMap[pHalfEdge->face->id];
					if(pHalfEdge->twin->face)
					{
						pPatch2 = mesh->patches + 
							pFacePatchMap[pHalfEdge->twin->face->id];
						GI_VEC3_ADD(v1, pHalfEdge->prev->vstart->coords, 
							pHalfEdge->twin->prev->vstart->coords);
					}
					else
						pPatch2 = NULL;
					GIHalfEdge_split(pHalfEdge, pPatch, pPatch2, 0.5, NULL);
					pFacePatchMap[pHalfEdge->next->twin->face->id] = pPatch->id;
					if(pPatch2)
					{
						vec = pHalfEdge->next->vstart->coords;
						GI_VEC3_SCALE(vec, vec, 0.5);
						GI_VEC3_ADD_SCALED(vec, vec, v1, 0.25);
						pFacePatchMap[pHalfEdge->twin->prev->twin->face->id] = pPatch2->id;
					}
				}
			}

			/* compute new positions for original vertices */
			for(i=0,pVertex=mesh->vertices; i<uiVCount; 
				++i,pVertex=pVertex->next)
			{
				pHalfEdge = pVertex->hedge;
				vec = pVertex->coords;
				if(pHalfEdge->face)
				{
					/* interior vertex */
					GIdouble dInvN, d1_N2, d4_N2, dN_3_N;
					GIuint uiDegree = 0;
					GI_VEC3_SET(v1, 0.0, 0.0, 0.0);
					GI_VEC3_SET(v2, 0.0, 0.0, 0.0);
					pHalfEdge = pHStart = pEdgeFlags[pHalfEdge->edge->id] ? 
						pHalfEdge->twin : pHalfEdge->prev;
					do
					{
						++uiDegree;
						GI_VEC3_SUB(v1, v1, pHalfEdge->vstart->coords);
						pHalfEdge = pHalfEdge->twin->prev;
						GI_VEC3_ADD(v2, v2, pHalfEdge->vstart->coords);
						pHalfEdge = pHalfEdge->twin->prev;
					}while(pHalfEdge != pHStart);
					dInvN = 1.0 / (GIdouble)uiDegree;
					d1_N2 = dInvN * dInvN;
					d4_N2 = 4.0 * d1_N2;
					dN_3_N = dInvN * (GIdouble)(uiDegree-3);
					GI_VEC3_SCALE(vec, vec, dN_3_N);
					GI_VEC3_ADD_SCALED(vec, vec, v1, d1_N2);
					GI_VEC3_ADD_SCALED(vec, vec, v2, d4_N2);
				}
				else
				{
					/* boundary vertex */
					GI_VEC3_ADD(v1, pHalfEdge->next->vstart->coords, 
						pHalfEdge->prev->vstart->coords);
					GI_VEC3_SCALE(vec, vec, 0.5);
					GI_VEC3_ADD_SCALED(vec, vec, v1, 0.25);
				}
			}
		}
		GI_FREE_ARRAY(pFacePatchMap);

		if(0)
		{
			GIdouble *pLimitCoords;
			pLimitCoords = (GIdouble*)GI_MALLOC_ARRAY(3*(uiVCount+
				(mesh->fcount>>3)), sizeof(GIdouble));

			/* compute limit positions for vertex-vertices */
			for(i=0,vec=pLimitCoords,pVertex=mesh->vertices; 
				i<uiVCount; ++i,vec+=3,pVertex=pVertex->next)
			{
				pHalfEdge = pVertex->hedge;
				if(pHalfEdge->face)
				{
					/* interior vertex */
					GIdouble dN2_NN_5, d4_NN_5, d1_NN_5;
					GIuint uiDegree = 0;
					GI_VEC3_SET(v1, 0.0, 0.0, 0.0);
					GI_VEC3_SET(v2, 0.0, 0.0, 0.0);
					pHalfEdge = pHStart = pEdgeFlags[pHalfEdge->edge->id] ? 
						pHalfEdge->twin : pHalfEdge->prev;
					do
					{
						++uiDegree;
						GI_VEC3_ADD(v1, v1, pHalfEdge->vstart->coords);
						pHalfEdge = pHalfEdge->twin->prev;
						GI_VEC3_ADD(v2, v2, pHalfEdge->vstart->coords);
						pHalfEdge = pHalfEdge->twin->prev;
					}while(pHalfEdge != pHStart);
					d1_NN_5 = 1.0 / (GIdouble)(uiDegree*(uiDegree+5));
					d4_NN_5 = 4.0 * d1_NN_5;
					dN2_NN_5 = (GIdouble)(uiDegree*uiDegree) * d1_NN_5;
					GI_VEC3_SCALE(vec, pVertex->coords, dN2_NN_5);
					GI_VEC3_ADD_SCALED(vec, vec, v1, d1_NN_5);
					GI_VEC3_ADD_SCALED(vec, vec, v2, d4_NN_5);
				}
				else
				{
					/* boundary vertex */
					GI_VEC3_ADD(vec, pHalfEdge->next->vstart->coords, 
						pHalfEdge->prev->vstart->coords);
					GI_VEC3_ADD_SCALED(vec, vec, pVertex->coords, 4.0);
					GI_VEC3_SCALE(vec, vec, d1_6);
				}
			}

			/* compute limit positions for face-vertices */
			uiVCount += mesh->fcount >> 3;
			for(; i<uiVCount; ++i,vec+=3,pVertex=pVertex->next)
			{
				pHalfEdge = pVertex->hedge;
				GI_VEC3_SET(v1, 0.0, 0.0, 0.0);
				GI_VEC3_SET(v2, 0.0, 0.0, 0.0);
				pHalfEdge = pHStart = pEdgeFlags[pHalfEdge->edge->id] ? 
					pHalfEdge->twin : pHalfEdge->prev;
				do
				{
					GI_VEC3_ADD(v1, v1, pHalfEdge->vstart->coords);
					pHalfEdge = pHalfEdge->twin->prev;
					GI_VEC3_ADD(v2, v2, pHalfEdge->vstart->coords);
					pHalfEdge = pHalfEdge->twin->prev;
				}while(pHalfEdge != pHStart);
				GI_VEC3_SCALE(vec, pVertex->coords, d16_36);
				GI_VEC3_ADD_SCALED(vec, vec, v1, d1_36);
				GI_VEC3_ADD_SCALED(vec, vec, v2, d4_36);
			}

			/* compute limit positions for edge-vertices */
			for(; pVertex!=mesh->vertices; pVertex=pVertex->next)
			{
				vec = pVertex->coords;
				pHalfEdge = pVertex->hedge;
				if(pHalfEdge->face)
				{
					/* interior vertex */
					GI_VEC3_SET(v1, 0.0, 0.0, 0.0);
					GI_VEC3_SET(v2, 0.0, 0.0, 0.0);
					pHalfEdge = pHalfEdge->twin;
					do
					{
						GI_VEC3_ADD(v1, v1, 
							pHalfEdge->prev->twin->prev->vstart->coords);
						GI_VEC3_ADD(v2, v2, pHalfEdge->vstart->coords);
						pHalfEdge = pHalfEdge->twin->prev;
					}while(pHalfEdge->twin != pVertex->hedge);
					GI_VEC3_SCALE(vec, pVertex->coords, d16_36);
					GI_VEC3_ADD_SCALED(vec, vec, v1, d1_36);
					GI_VEC3_ADD_SCALED(vec, vec, v2, d4_36);
				}
				else
				{
					/* boundary vertex */
					GI_VEC3_SCALE(vec, vec, 4.0);
					GI_VEC3_ADD(vec, vec, pHalfEdge->next->vstart->coords);
					GI_VEC3_ADD(vec, vec, pHalfEdge->prev->vstart->coords);
					GI_VEC3_SCALE(vec, vec, d1_6);
				}
			}

			/* set limit positions for vertex-vertices and face-vertices */
			for(i=0,vec=pLimitCoords; i<uiVCount; 
				++i,vec+=3,pVertex=pVertex->next)
			{
				GI_VEC3_COPY(pVertex->coords, vec);
			}
			GI_FREE_ARRAY(pLimitCoords);
		}
		GI_FREE_ARRAY(pEdgeFlags);
	}

	/* compute new edge lengths */
	GI_LIST_FOREACH(mesh->edges, pEdge)
		pEdge->length = GIvec3d_dist(pEdge->hedge[0].vstart->coords, 
			pEdge->hedge[1].vstart->coords);
	GI_LIST_NEXT(mesh->edges, pEdge)
	for(i=0; i<mesh->patch_count; ++i)
	{
		pPatch = mesh->patches + i;
		pParam = pPatch->corners[0];
		for(j=0; j<4; ++j)
			for(; pParam!=pPatch->corners[(j+1)&3]; 
				pParam=pParam->cut_hedge->next->pstart)
				pPatch->side_lengths[j] += pParam->cut_hedge->edge->length;
	}

	/* compute paths */
	--uiResolution;
	GICutter_compute_paths(cutter, mesh);
	for(i=0; i<mesh->patch_count; ++i)
	{
		GICutPath *pPath;
		pPatch = mesh->patches + i;
		GI_LIST_FOREACH(pPatch->paths, pPath)
			pPath->glength = uiResolution;
		GI_LIST_NEXT(pPatch->paths, pPath)
		GIPatch_renumerate_params(pPatch);
	}
	mesh->cut_splits = mesh->split_hedges.size;
	mesh->varray_attribs = 0;
	return mesh->patch_count;
}

/** \internal
 *  \brief Create parameter coordinates from patchification.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \param patch_ids per-face patch IDs
 *  \param cut_flags per-edge cut flags
 *  \retval GI_TRUE if created successfully
 *  \retval GI_FALSE on error
 *  \ingroup cutting
 */
GIboolean GICutter_create_params(GICutter *cutter, GIMesh *mesh, 
								 GIint *patch_ids, GIubyte *cut_flags)
{
	GIPatch *pPatch;
	GIHalfEdge *pHalfEdge, *pHStart;
	GIVertex *pVertex, dummyVertex;
	GIParam *pParam, dummyParam;
	GIParam **pRoots = (GIParam**)GI_MALLOC_ARRAY(
		mesh->patch_count, sizeof(GIParam*));
	GIint i;

	/* initialize auxiliary data structures */
	dummyVertex.cut_degree = 0;
	dummyParam.vertex = &dummyVertex;
	pRoots = (GIParam**)GI_MALLOC_ARRAY(mesh->patch_count, sizeof(GIParam*));
	for(i=0; i<mesh->patch_count; ++i)
		pRoots[i] = &dummyParam;

	/* traverse vertices and create params */
	GI_LIST_FOREACH(mesh->vertices, pVertex)
		pHStart = pVertex->hedge;
		if(pVertex->cut_degree)
			while(!cut_flags[pHStart->edge->id])
				pHStart = pHStart->prev->twin;
		pHalfEdge = pHStart->face ? pHStart : pHStart->prev->twin;
		pPatch = mesh->patches + patch_ids[pHalfEdge->face->id];
		pParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
		GI_LIST_ADD(pPatch->params, pParam);
		pParam->id = pPatch->pcount++;
		GI_VEC2_COPY(pParam->params, (GIfloat*)pHalfEdge->pstart);
		pParam->stretch = 0.0;
		pParam->vertex = pVertex;
		pParam->cut_hedge = cut_flags[pHalfEdge->edge->id] ? pHalfEdge : NULL;
		pHalfEdge->pstart = pParam;
		pHalfEdge = pHalfEdge->prev->twin;

		/* traverse half edges counter-clockwise */
		while(pHalfEdge != pHStart)
		{
			if(cut_flags[pHalfEdge->edge->id])
			{
				pPatch = mesh->patches + patch_ids[pHalfEdge->face->id];
				pParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
				GI_LIST_ADD(pPatch->params, pParam);
				pParam->id = pPatch->pcount++;
				GI_VEC2_COPY(pParam->params, (GIfloat*)pHalfEdge->pstart);
				pParam->stretch = 0.0;
				pParam->vertex = pVertex;
				pParam->cut_hedge = pHalfEdge;
			}
			pHalfEdge->pstart = pParam;
			pHalfEdge = pHalfEdge->prev->twin;
		}
		if(!pHalfEdge->face)
			pHalfEdge->pstart = NULL;

		/* no loose ends in multi-chart */
		if(mesh->patch_count > 1 && pVertex->cut_degree == 1)
			return GI_FALSE;

		/* remember vertex with highest cut degree for every patch */
		if(mesh->patch_count == 1)
		{
			if(pVertex->cut_degree > pRoots[0]->vertex->cut_degree)
				pRoots[0] = pVertex->hedge->pstart;
		}
		else if(pVertex->cut_degree)
		{
			pHalfEdge = pVertex->hedge->face ? pVertex->hedge : 
				pVertex->hedge->twin->next;
			do
			{
				if(cut_flags[pHalfEdge->edge->id] && pVertex->cut_degree > 
				   pRoots[patch_ids[pHalfEdge->face->id]]->vertex->cut_degree)
					pRoots[patch_ids[pHalfEdge->face->id]] = pHalfEdge->pstart;
				pHalfEdge = pHalfEdge->twin->next;
			}while(pHalfEdge != pVertex->hedge);
		}
	GI_LIST_NEXT(mesh->vertices, pVertex)

	/* set roots and check if cut existing */
	if(!pRoots[0]->vertex->cut_degree)
		return GI_FALSE;
	for(i=0; i<mesh->patch_count; ++i)
		mesh->patches[i].params = pRoots[i];
	GI_FREE_ARRAY(pRoots);
	return GI_TRUE;
}

/** \internal
 *  \brief Straighten cut paths.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \param root vertex at cut intersection
 *  \param cut_flags per-edge cut flags
 *  \ingroup cutting
 */
void GICutter_straighten_paths(GICutter *cutter, GIMesh *mesh, 
							   GIVertex *root, GIubyte *cut_flags)
{
	GIHalfEdge *pHalfEdge;
	GIVertex *pVEnd, *pVStart, *pVertex;
	GIDynamicQueue qCutPaths;
	GIHeap qFringe;
	GIHash hPathInfos;
	GIPathInfo *pInfo;
	GIdouble *vend;
	GIdouble dDist, dOldDist;
	if(!root)
		return;

	/* intialize data structures */
	GIHeap_construct(&qFringe, mesh->vcount, lessd, -DBL_MAX, GI_TRUE);
	GIHash_construct(&hPathInfos, (GIuint)sqrt((GIdouble)mesh->vcount), 
		0.0f, sizeof(GIuint), hash_uint, compare_uint, copy_uint);
	GIDynamicQueue_construct(&qCutPaths);
	pHalfEdge = root->hedge;
	do
	{
		if(cut_flags[pHalfEdge->edge->id] == GI_EDGE_USEABLE)
			GIDynamicQueue_enqueue(&qCutPaths, pHalfEdge);
		pHalfEdge = pHalfEdge->twin->next;
	}while(pHalfEdge != root->hedge);

	/* traverse all cut paths */
	while(qCutPaths.size)
	{
		pHalfEdge = (GIHalfEdge*)GIDynamicQueue_dequeue(&qCutPaths);
		if(cut_flags[pHalfEdge->edge->id] == GI_EDGE_USED)
			continue;
		else if(!pHalfEdge->face || !pHalfEdge->twin->face)
		{
			/* boundary definitely on cut */
			pVEnd = pHalfEdge->next->vstart;
			while(pVEnd->cut_degree == 2)
			{
				cut_flags[pHalfEdge->edge->id] = GI_EDGE_USED;
				pHalfEdge = pHalfEdge->next;
				while(!cut_flags[pHalfEdge->edge->id])
					pHalfEdge = pHalfEdge->twin->next;
				pVEnd = pHalfEdge->next->vstart;
			}
			cut_flags[pHalfEdge->edge->id] = GI_EDGE_USED;
		}
		else
		{
			/* remove old cut path and find end vertex */
			pVStart = pHalfEdge->vstart;
			pVEnd = pHalfEdge->next->vstart;
			while(pVEnd->cut_degree == 2)
			{
				cut_flags[pHalfEdge->edge->id] = GI_NONE;
				pVEnd->cut_degree = 0;
				pHalfEdge = pHalfEdge->next;
				while(!cut_flags[pHalfEdge->edge->id])
					pHalfEdge = pHalfEdge->twin->next;
				pVEnd = pHalfEdge->next->vstart;
			}
			cut_flags[pHalfEdge->edge->id] = GI_NONE;
			vend = pVEnd->coords;

			/* find new cut path by A* */
			GIHeap_clear(&qFringe);
			GIHash_clear(&hPathInfos, sizeof(GIPathInfo));
			GIHash_insert(&hPathInfos, &pVStart->id, 
				GI_CALLOC_SINGLE(sizeof(GIPathInfo)));
			pVertex = pVStart;
			dOldDist = 0.0;
			while(pVertex != pVEnd)
			{
				pHalfEdge = pVertex->hedge;
				do
				{
					GIVertex *pVOther = pHalfEdge->next->vstart;
					if(!pVOther->cut_degree || pVOther == pVEnd)
					{
						dDist = dOldDist + pHalfEdge->edge->length;
						pInfo = (GIPathInfo*)GIHash_find(&hPathInfos, &pVOther->id);
						if(!pInfo || dDist < pInfo->distance)
						{
							if(!pInfo)
								GIHash_insert(&hPathInfos, &pVOther->id, pInfo=
									(GIPathInfo*)GI_MALLOC_SINGLE(sizeof(GIPathInfo)));
							pInfo->source = pHalfEdge;
							pInfo->distance = dDist;
							GIHeap_enqueue(&qFringe, pVOther, 
								dDist+GIvec3d_dist(pVOther->coords, vend));
						}
					}
					pHalfEdge = pHalfEdge->twin->next;
				}while(pHalfEdge != pVertex->hedge);
				pVertex = (GIVertex *)GIHeap_dequeue(&qFringe, &dOldDist);
			}

			/* retrace and assemble new cut path */
			pHalfEdge = ((GIPathInfo*)GIHash_find(
				&hPathInfos, &pVEnd->id))->source;
			pVertex = pHalfEdge->vstart;
			while(pVertex != pVStart)
			{
				pVertex->cut_degree = 2;
				cut_flags[pHalfEdge->edge->id] = GI_EDGE_USED;
				pHalfEdge = ((GIPathInfo*)GIHash_find(
					&hPathInfos, &pVertex->id))->source;
				pVertex = pHalfEdge->vstart;
			}
			cut_flags[pHalfEdge->edge->id] = GI_EDGE_USED;
		}

		/* enqueue new found cut paths (if not already traversed) */
		pHalfEdge = pVEnd->hedge;
		do
		{
			if(cut_flags[pHalfEdge->edge->id] == GI_EDGE_USEABLE)
				GIDynamicQueue_enqueue(&qCutPaths, pHalfEdge);
			pHalfEdge = pHalfEdge->twin->next;
		}while(pHalfEdge != pVEnd->hedge);
	}

	/* clean up */
	GIHeap_destruct(&qFringe);
	GIHash_destruct(&hPathInfos, sizeof(GIPathInfo));
	GIDynamicQueue_destruct(&qCutPaths);
}

/** \internal
 *  \brief Sort faces into patches.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \param patch_ids per-face patch IDs
 *  \param cut_flags per-edge cut flags
 *  \ingroup cutting
 */
void GICutter_sort_faces(GICutter *cutter, GIMesh *mesh, 
						 GIint *patch_ids, GIubyte *cut_flags)
{
	GIPatch *pPatch;
	GIFace *pFace, *pFaces = NULL;
	GIHalfEdge *pHalfEdge;
	GIDynamicQueue qFaces;
	GIuint i;

	/* sort faces into patches by region growing */
	GIDynamicQueue_construct(&qFaces);
	for(i=0; i<mesh->fcount; ++i)
		patch_ids[i] = -1;
	for(i=0; mesh->faces; ++i)
	{
		pPatch = mesh->patches + i;
		pPatch->faces = mesh->faces;
		patch_ids[mesh->faces->id] = i;
		GIDynamicQueue_enqueue(&qFaces, mesh->faces);
		while(qFaces.size)
		{
			pFace = (GIFace*)GIDynamicQueue_dequeue(&qFaces);
			GI_LIST_REMOVE(mesh->faces, pFace);
			GI_LIST_ADD(pFaces, pFace);
			++pPatch->fcount;

			/* look at neighbouring faces */
			GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
				if(!cut_flags[pHalfEdge->edge->id] && 
					patch_ids[pHalfEdge->twin->face->id] == -1)
				{
					patch_ids[pHalfEdge->twin->face->id] = i;
					GIDynamicQueue_enqueue(&qFaces, pHalfEdge->twin->face);
				}
			GI_LIST_NEXT(pFace->hedges, pHalfEdge)
		}
	}
	GIDynamicQueue_destruct(&qFaces);
}

/** \internal
 *  \brief Compute cut paths.
 *  \param cutter cutter to use
 *  \param mesh mesh to work on
 *  \ingroup cutting
 */
void GICutter_compute_paths(GICutter *cutter, GIMesh *mesh)
{
	GIPatch *pPatch;
	GICutPath *pPath, *pTwin;
	GIHalfEdge *pHalfEdge, *pCutHedge;
	GIVertex *pVertex, *pVRoot;
	GIParam *pParam;
	GIHash hCutPaths;
	GIuint i;

	/* auxiliary datastructures */
	GIHash_construct(&hCutPaths, 32*mesh->patch_count, 0.0f, 
		sizeof(GIHalfEdge*), hash_uint, compare_uint, copy_uint);

	/* traverse patches and compute path information */
	for(i=0; i<mesh->patch_count; ++i)
	{
		/* clear paths */
		pPatch = mesh->patches + i;
		pPatch->hcount = 0;
		pPatch->hlength = 0.0;
		pPatch->path_count = 0;
		pPatch->groups = 0;
		GI_LIST_CLEAR_PERSISTENT(pPatch->paths, sizeof(GICutPath));
		GIDynamicQueue_clear(&pPatch->split_paths);

		/* accumulate path information */
		pParam = pPatch->params;
		pVRoot = pParam->vertex;
		pPath = NULL;
		do
		{
			pCutHedge = pParam->cut_hedge;
			pVertex = pParam->vertex;
			if(pVertex->cut_degree != 2 || 
				(pVertex->flags & GI_VERTEX_EXACT_BIT) || pVertex == pVRoot)
			{
				/* add length and insert path in hash */
				if(pPath)
				{
					pPatch->hlength += pPath->elength;
					if(pHalfEdge->twin->face)
						GIHash_insert(&hCutPaths, &pHalfEdge->edge->id, pPath);
				}

				/* cut hedge allready in hash */
				pPath = (GICutPath*)GI_MALLOC_PERSISTENT(sizeof(GICutPath));
				GI_LIST_ADD(pPatch->paths, pPath);
				pPath->id = pPatch->path_count++;
				pPath->patch = pPatch;
				pPath->glength = 0;
				pPath->pstart = pParam;
				pTwin = (GICutPath*)GIHash_find(&hCutPaths, &pCutHedge->edge->id);
				if(pTwin)
				{
					/* copy data */
					if(pTwin->patch == pPatch)
						pPath->group = pTwin->group;
					else
						pPath->group = pPatch->groups++;
					pPath->twin = pTwin;
					pTwin->twin = pPath;
					pPath->elength = pTwin->elength;
					pPatch->hlength += pPath->elength;
					pPath = NULL;
				}
				else
				{
					/* new path */
					pPath->group = pPatch->groups++;
					pPath->elength = pCutHedge->edge->length;
					pPath->twin = NULL;
				}
			}
			else if(pPath)
				pPath->elength += pCutHedge->edge->length;
			++pPatch->hcount;
			pHalfEdge = pCutHedge;
			pParam = pCutHedge->next->pstart;
		}while(pParam != pPatch->params);

		/* finish last path */
		if(pPath)
		{
			pPatch->hlength += pPath->elength;
			if(pHalfEdge->twin->face)
				GIHash_insert(&hCutPaths, &pHalfEdge->edge->id, pPath);
		}
	}

	/* clean up */
	GIHash_destruct(&hCutPaths, 0);
}

/** \internal
 *  \brief Patch destructor.
 *  \param patch patch to destruct
 *  \ingroup cutting
 */
void GIPatch_destruct(GIPatch *patch)
{
	/* clear lists */
	GI_LIST_CLEAR_PERSISTENT(patch->paths, sizeof(GICutPath));
	GI_LIST_CLEAR_PERSISTENT(patch->params, sizeof(GIParam));
	GIDynamicQueue_destruct(&patch->split_paths);
}

/** \internal
 *  \brief Prevent degenerate regions in parameter space by splitting edges.
 *  \param patch patch to work on
 *  \ingroup cutting
 */
void GIPatch_prevent_singularities(GIPatch *patch)
{
	GIHalfEdge *pHalfEdge, *pCutHedge;
	GIParam *pParam;

	/* prevent degenerate regions (not only triangles) */
	pHalfEdge = patch->params->cut_hedge->next;
	do
	{
		pParam = pHalfEdge->pstart;
		pCutHedge = pParam->cut_hedge;
		while(pHalfEdge != pCutHedge)
		{
			if(pHalfEdge->next->vstart->cut_degree)
				GIHalfEdge_split(pHalfEdge, patch, patch, 0.5, NULL);
			pHalfEdge = pHalfEdge->twin->next;
		}
		pHalfEdge = pCutHedge->next;
	}while(pParam != patch->params);
}

/** \internal
 *  \brief Numerate parameter coords so that all interior params come first
 *  \param patch patch to work on
 *  \ingroup cutting
 */
void GIPatch_renumerate_params(GIPatch *patch)
{
	GIParam *pParam;
	GIuint i = 0, b = patch->pcount - patch->hcount;

	/* numerate */
	GI_LIST_FOREACH(patch->params, pParam)
		if(pParam->cut_hedge)
			pParam->id = b++;
		else
			pParam->id = i++;
	GI_LIST_NEXT(patch->params, pParam)
}

/** \internal
 *  \brief Search for vertices with corner flag set.
 *  \param patch patch to work on
 *  \retval GI_TRUE if patch has 4 valid corners
 *  \retval GI_FALSE else
 *  \ingroup cutting
 */
GIboolean GIPatch_find_corners(GIPatch *patch)
{
	GICutPath *pPath, *pStartPath;
	GIHalfEdge *pHalfEdge;
	GIParam *pParam = patch->params;
	GIdouble dLength;
	GIint iCorners = 0;
	if(patch->corners[0])
		return GI_TRUE;
	memset(patch->side_lengths, 0, 4*sizeof(GIdouble));

	/* find corners */
	do
	{
		if(pParam->vertex->flags & GI_VERTEX_CORNER_BIT)
		{
			if(iCorners < 4)
				patch->corners[iCorners] = pParam;
			++iCorners;
		}
		pParam = pParam->cut_hedge->next->pstart;
	}while(pParam != patch->params);

	/* too less or two many corners */
	if(iCorners != 4)
	{
		memset(patch->corners, 0, 4*sizeof(GIParam*));
		return GI_FALSE;
	}

	/* compute side lengths and split paths if neccessary */
	iCorners = pParam==patch->corners[0] ? 1 : 0;
	pStartPath = patch->paths;
	GI_LIST_FOREACH(patch->paths, pPath)
		dLength = 0.0;
		do
		{
			pHalfEdge = pParam->cut_hedge;
			pParam = pHalfEdge->next->pstart;
			dLength += pHalfEdge->edge->length;
			if(iCorners)
				patch->side_lengths[iCorners-1] += pHalfEdge->edge->length;
			if(pParam == patch->corners[iCorners])
			{
				if(GIPatch_split_path(patch, pPath, pParam))
				{
					pPath->next->elength = pPath->elength - dLength;
					pPath->elength = dLength;
					if(pPath->twin)
					{
						if(pPath->twin->patch != patch)
						{
							GIPatch_split_path(pPath->twin->patch, 
								pPath->twin, pParam->cut_hedge->twin->next->pstart);
							pPath->twin->twin = pPath->next;
							pPath->twin = pPath->twin->next;
						}
						pPath->twin->elength = pPath->elength;
						pPath->next->twin->elength = pPath->next->elength;
					}
				}
				if(!iCorners)
					pStartPath = pPath->next;
				if(++iCorners == 4)
				{
					pPath = patch->paths->prev;
					break;
				}
			}
		}while(pParam != pPath->next->pstart);
	GI_LIST_NEXT(patch->paths, pPath)

	/* cut starts at lower left corner */
	patch->side_lengths[3] = patch->hlength - patch->side_lengths[0] - 
		patch->side_lengths[1] - patch->side_lengths[2];
	patch->fixed_corners = GI_TRUE;
	patch->params = patch->corners[0];
	patch->paths = pStartPath;
	return GI_TRUE;
}

/** \internal
 *  \brief Check if patch has bijective mapping onto the unit square.
 *  \param patch patch to check
 *  \retval GI_TRUE if patch has valid parameterization
 *  \retval GI_FALSE else
 *  \ingroup cutting
 */
GIboolean GIPatch_valid_parameterization(GIPatch *patch)
{
	static const GIdouble dEPS = 1e-3;
	static const GIdouble dEPS2 = 1e-6;
	static const GIdouble corner[8] = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
	GIParam *pParam = patch->params;
	GIdouble dDist, dOldDist;
	GIint i, j, c, iInc = 1, iMod = 0;
	memset(patch->corners, 0, 4*sizeof(GIParam*));
	if(!pParam)
		return GI_FALSE;

	/* find lower left corner */
	do
	{
		dDist = GIvec2d_dist_sqr(pParam->params, (GIvec2d)corner);
		if(patch->corners[0] && dDist > dOldDist)
			break;
		if(dDist < dEPS2)
			patch->corners[0] = pParam;
		dOldDist = dDist;
		pParam = pParam->cut_hedge->next->pstart;
	}while(pParam != patch->params);
	if(!patch->corners[0])
		return GI_FALSE;

	/* which direction? */
	if(GIvec2d_dist_sqr(pParam->params, ((GIvec2d)corner)+2) > 
		GIvec2d_dist_sqr(pParam->params, ((GIvec2d)corner)+6))
	{
		iInc = -1;
		iMod = 1;
	}
	i = (4+iInc) & 3;
	c = (i+iMod) & 1;
	j = i << 1;

	/* find other corners and check boundary consistency */
	dOldDist = GIvec2d_dist_sqr(pParam->params, ((GIvec2d)corner)+j);
	while(pParam != patch->corners[0])
	{
		dDist = GIvec2d_dist_sqr(pParam->params, ((GIvec2d)corner)+j);
		if(dDist > dOldDist)
		{
			if(patch->corners[i])
			{
				i = (i+4+iInc) & 3;
				c = (i+iMod) & 1;
				j = i << 1;
				dOldDist = GIvec2d_dist_sqr(pParam->params, ((GIvec2d)corner)+j);
				continue;
			}
			else
				break;
		}
		if(fabs(pParam->params[c]-corner[j+c]) > dEPS)
			break;
		if(dDist < dEPS2)
			patch->corners[i] = pParam;
		dOldDist = dDist;
		pParam = pParam->cut_hedge->next->pstart;
	}

	/* error encountered? */
	if(pParam != patch->corners[0])
	{
		memset(patch->corners, 0, 4*sizeof(GIParam*));
		return GI_FALSE;
	}
	patch->params = patch->corners[0];
	return GI_TRUE;
}

/** \internal
 *  \brief Split cut path.
 *  \param patch patch path belongs to
 *  \param path path to split
 *  \param param param to split at
 *  \retval GI_TRUE if path split
 *  \retval GI_FALSE if split point is allready a cut node
 *  \ingroup cutting
 */
GIboolean GIPatch_split_path(GIPatch *patch, GICutPath *path, GIParam *param)
{
	GICutPath *pNew, *pTwin = path->twin, *pNew2 = NULL;
	if(param == path->pstart || param == path->next->pstart)
		return GI_FALSE;

	/* create path */
	pNew = (GICutPath*)GI_MALLOC_PERSISTENT(sizeof(GICutPath));
	GI_LIST_INSERT(patch->paths, path->next, pNew);
	pNew->id = patch->path_count++;
	pNew->patch = patch;
	pNew->group = patch->groups++;
	pNew->pstart = param;
	pNew->twin = pTwin;

	/* split twin if same patch */
	if(pTwin && pTwin->patch == patch)
	{
		pNew2 = (GICutPath*)GI_MALLOC_PERSISTENT(sizeof(GICutPath));
		GI_LIST_INSERT(patch->paths, pTwin->next, pNew2);
		pNew2->id = patch->path_count++;
		pNew2->patch = patch;
		pNew2->group = path->group;
		pNew2->pstart = param->cut_hedge->twin->next->pstart;
		pNew2->twin = path;
		path->twin = pNew2;
		pTwin->twin = pNew;
		pTwin->group = pNew->group;
	}
	GIDynamicQueue_push(&patch->split_paths, path);
	return GI_TRUE;
}

/** \internal
 *  \brief Revert all path splits.
 *  \param patch patch to work on
 *  \param count number of splits to reverse or -1 for all
 *  \ingroup cutting
 */
void GIPatch_revert_splits(GIPatch *patch, GIint count)
{
	GICutPath *pPath, *pPath2, *pTwin;
	if(count < 0 || count > patch->split_paths.size)
		count = patch->split_paths.size;

	/* reverse splits */
	while(count--)
	{
		pPath = (GICutPath*)GIDynamicQueue_pop(&patch->split_paths);
		pPath2 = pPath->next;
		pTwin = pPath2->twin;

		/* delete path */
		pPath->elength += pPath2->elength;
		pPath->glength += pPath2->glength;
		GI_LIST_DELETE_PERSISTENT(patch->paths, pPath2, sizeof(GICutPath));
		--patch->path_count;
		--patch->groups;

		/* reverse twin if same patch */
		if(pTwin && pTwin->patch == patch)
		{
			pPath2 = pTwin->next;
			GI_LIST_DELETE_PERSISTENT(patch->paths, pPath2, sizeof(GICutPath));
			--patch->path_count;
			pTwin->group = pPath->group;
			pTwin->elength = pPath->elength;
			pTwin->glength = pPath->glength;
			pPath->twin = pTwin;
			pTwin->twin = pPath;
		}
	}
}

/** \internal
 *  \brief Compute stretch of patch.
 *  \param patch patch to work on
 *  \param metric stretch metric to use
 *  \param param_stretches GI_TRUE to compute param stretches GI_FALSE else
 *  \param init_stretches GI_FALSE if param stretches all 0, GI_TRUE else
 *  \return face with maximal stretch
 *  \ingroup cutting
 */
GIFace* GIPatch_compute_stretch(GIPatch *patch, GIuint metric, 
							   GIboolean param_stretches, 
							   GIboolean init_stretches)
{
	GIFace *pFace = patch->faces, *pFEnd = patch->next->faces, *pMaxFace;
	GIHalfEdge *pHalfEdge;
	GIParam *pParam;
	GIdouble dWeight = patch->mesh->context->parameterizer.area_weight;
	GIdouble dA2D, dA3D, dStretch, dS2A;
	GIdouble dStretchSum = 0.0, dA2DSum = 0.0, dA3DSum = 0.0;
	GIdouble dMaxFaceStretch = 0.0;
	GIdouble *pParamA3DSums = NULL;
	GIvoid *pArgs = NULL;

	/* error checking and initialization */
	switch(metric)
	{
	case GI_RMS_GEOMETRIC_STRETCH:
	case GI_COMBINED_STRETCH:
		if(param_stretches)
			pParamA3DSums = (GIdouble*)GI_CALLOC_ARRAY(patch->pcount, sizeof(GIdouble));
		pArgs = &dWeight;
		break;
	case GI_MAX_GEOMETRIC_STRETCH:
		patch->stretch[GI_MAX_GEOMETRIC_STRETCH-GI_STRETCH_BASE] = 0.0;
		break;
	default:
		return NULL;
	}
	if(param_stretches)
	{
		patch->min_param_stretch = DBL_MAX;
		patch->max_param_stretch = 0.0;
		if(init_stretches)
		{
			GI_LIST_FOREACH(patch->params, pParam)
				pParam->stretch = 0.0;
			GI_LIST_NEXT(patch->params, pParam)
		}
	}

	/* compute stretch for each face */
	pMaxFace = pFace;
	do
	{
		/* compute stretch and areas */
		dStretch = GIFace_stretch(pFace, metric, pArgs, &dA2D);
		if(!patch->surface_area)
		{
			dA3D = GIFace_area(pFace);
			dA3DSum += dA3D;
		}
		if(!patch->param_area)
			dA2DSum += dA2D;

		/* accumulate values */
		switch(metric)
		{
		/* L2-stretch or combined energy */
		case GI_RMS_GEOMETRIC_STRETCH:
		case GI_COMBINED_STRETCH:
			if(patch->surface_area)
				dA3D = GIFace_area(pFace);
			dS2A = dStretch * dA3D;
			dStretchSum += dS2A;
			if(param_stretches)
			{
				GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
					pParam = pHalfEdge->pstart;
					pParam->stretch += dS2A;
					if(pParamA3DSums)
						pParamA3DSums[pParam->id] += dA3D;
				GI_LIST_NEXT(pFace->hedges, pHalfEdge)
			}
			break;

		/* Linf-stretch */
		case GI_MAX_GEOMETRIC_STRETCH:
			if(dStretch > patch->stretch[GI_MAX_GEOMETRIC_STRETCH-GI_STRETCH_BASE])
				patch->stretch[GI_MAX_GEOMETRIC_STRETCH-GI_STRETCH_BASE] = dStretch;
			if(param_stretches)
			{
				GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
					pParam = pHalfEdge->pstart;
					if(dStretch > pParam->stretch)
						pParam->stretch = dStretch;
					if(dStretch > patch->max_param_stretch)
						patch->max_param_stretch = dStretch;
					if(dStretch < patch->min_param_stretch)
						patch->min_param_stretch = dStretch;
				GI_LIST_NEXT(pFace->hedges, pHalfEdge)
			}
		}

		/* face with maximum stretch */
		if(dStretch > dMaxFaceStretch)
		{
			pMaxFace = pFace;
			dMaxFaceStretch = dStretch;
		}
		pFace = pFace->next;
	}while(pFace != pFEnd);

	/* compute overall stretch and param stretch */
	if(!patch->surface_area)
		patch->surface_area = dA3DSum;
	if(!patch->param_area)
		patch->param_area = dA2DSum;
	if(metric == GI_RMS_GEOMETRIC_STRETCH || metric == GI_COMBINED_STRETCH)
	{
		patch->stretch[metric-GI_STRETCH_BASE] = 
			sqrt(dStretchSum/patch->surface_area);
		if(param_stretches)
		{
			GI_LIST_FOREACH(patch->params, pParam)
				pParam->stretch = sqrt(pParam->stretch/pParamA3DSums[pParam->id]);
				if(pParam->stretch > patch->max_param_stretch)
					patch->max_param_stretch = pParam->stretch;
				if(pParam->stretch < patch->min_param_stretch)
					patch->min_param_stretch = pParam->stretch;
			GI_LIST_NEXT(patch->params, pParam)
		}
	}

	/* clean up */
	if(param_stretches && patch->param_metric != metric)
	{
		patch->param_metric = metric;
		if(patch->mesh->patch_count == 1)
		{
			if(patch->min_param_stretch < patch->mesh->min_param_stretch)
				patch->mesh->min_param_stretch = patch->min_param_stretch;
			if(patch->max_param_stretch > patch->mesh->max_param_stretch)
				patch->mesh->max_param_stretch = patch->max_param_stretch;
			patch->mesh->param_metric = metric;
		}
		else
			patch->mesh->param_metric = 0;
	}
	if(pParamA3DSums)
		GI_FREE_ARRAY(pParamA3DSums);
	return pMaxFace;
}
