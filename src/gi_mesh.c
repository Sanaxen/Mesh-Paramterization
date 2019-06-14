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
 *  \brief Implementation of structures and functions for mesh handling.
 */

#include "gi_mesh.h"
#include "gi_cutter.h"
#include "gi_context.h"
#include "gi_memory.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
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


/** \internal
 *  \brief Information about half edge split.
 */
typedef struct _GISplitInfo
{
	GIHalfEdge	*hedge;							/**< Split half edge. */
	GIPatch		*patch;							/**< Patch half edge belongs to. */
	GIPatch		*twin_patch;					/**< Patch twin half edge belongs to (if any). */
} GISplitInfo;

/** \internal
 *  \brief Hash packed attributes.
 *  \param a attribute package to hash
 *  \param size size of hash table
 *  \return computed hash value
 */
static GIuint hash_attribs(const GIvoid *a, GIuint size)
{
/*	const GIubyte *c = (const GIubyte*)a + sizeof(GIuint);
	const GIubyte *end = c + *(const GIuint*)a;
	GIuint uiHash = 2166136261U;
	while(c != end)
	{
		uiHash ^= (GIuint)(*(c++));
		uiHash *= 16777619;
	}
	return uiHash % size;
*/	GIuint i = 0;
	const GIuint *p = (const GIuint*)a;
	p += *p / sizeof(GIuint);
	for(; p!=a; --p)
		i ^= *p;
	return i % size;
}

/** \internal
 *  \brief Copy packed attributes.
 *  \param d destination
 *  \param s source
 */
static void copy_attribs(GIvoid *d, const GIvoid *s)
{
	memcpy(d, s, *(const GIuint*)s+sizeof(GIuint));
}

/** \internal
 *  \brief Compare packed attributes.
 *  \param a first attribute package
 *  \param b second attribute package
 *  \retval GI_TRUE if packages are equal
 *  \retval GI_FALSE if packages are not equal
 */
static GIboolean compare_attribs(const GIvoid *a, const GIvoid *b)
{
	return !memcmp((const GIuint*)a+1, (const GIuint*)b+1, *(const GIuint*)a);
}

/** \internal
 *  \brief Compare unsigned ints for qsort and bsearch.
 *  \param a first value
 *  \param b second value
 *  \return negative value if a < b, positive value if a > b and 0 if equal
 */
static int compare(const void *a, const void *b)
{
	return *(const GIuint*)a - *(const GIuint*)b;
}

/** Create new mesh object.
 *  \return id of new mesh
 *  \ingroup mesh
 */
GIuint GIAPIENTRY giGenMesh()
{
	GIContext *pContext = GIContext_current();
	GIint a;

	/* create mesh and add to hash */
	GIMesh *pMesh = (GIMesh*)GI_CALLOC_SINGLE(sizeof(GIMesh));
	pMesh->id = pContext->next_mid++;
	pMesh->context = pContext;
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
		pMesh->aoffset[a] = -1;
	pMesh->genus = -1;
	GIHash_insert(&pContext->mesh_hash, &pMesh->id, pMesh);
	return pMesh->id;
}

/** Create new nesh objects.
 *  \param n number of meshes to create
 *  \param meshes array to store ids of newly created meshes
 *  \ingroup mesh
 */
void GIAPIENTRY giGenMeshes(GIsizei n, GIuint *meshes)
{
	GIint i;

	/* create meshes */
	for(i=0; i<n; ++i)
		meshes[i] = giGenMesh();
}

/** Bind mesh to use for further processing.
 *  \param mesh id of mesh to bind
 *  \ingroup mesh
 */
void GIAPIENTRY giBindMesh(GIuint mesh)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh;
	if(pContext->mesh && pContext->mesh->id == mesh)
		return;
	if(!mesh)
	{
		pContext->mesh = NULL;
		return;
	}

	/* search and bind mesh */
	pMesh = (GIMesh*)GIHash_find(&pContext->mesh_hash, &mesh);
	if(!pMesh)
		GIContext_error(pContext, GI_INVALID_ID);
	else
		pContext->mesh = pMesh;
}

/** Delete mesh and all associated data.
 *  \param mesh id of mesh to delete
 *  \ingroup mesh
 */
void GIAPIENTRY giDeleteMesh(GIuint mesh)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh;

	/* find mesh */
	pMesh = (GIMesh*)GIHash_remove(&pContext->mesh_hash, &mesh);
	if(!pMesh)
	{
		GIContext_error(pContext, GI_INVALID_ID);
		return;
	}

	/* remove and clean up */
	if(pMesh == pContext->mesh)
		pContext->mesh = NULL;
	if(mesh == pContext->next_mid-1)
		--pContext->next_mid;
	GIMesh_destruct(pMesh);
	GI_FREE_SINGLE(pMesh, sizeof(GIMesh));
}

/** Delete mesh objects.
 *  \param n number of meshes to delete
 *  \param meshes ids of meshes to delete
 *  \ingroup mesh
 */
void GIAPIENTRY giDeleteMeshes(GIsizei n, const GIuint *meshes)
{
	GIint i;

	/* delete meshes */
	for(i=0; i<n; ++i)
		giDeleteMesh(meshes[i]);
}

/** Create mesh from indices into current attribute arrays.
 *  This function creates a new mesh by indexing into the currently set and 
 *  enabled attribute arrays. WARNING: All previous mesh data of the current 
 *  bound mesh object will be deleted.
 *  \param start minimal value in index array
 *  \param end maximal value in index array
 *  \param count size of index array
 *  \param indices array of indices into attribute arrays
 *  \ingroup mesh
 */
void GIAPIENTRY giIndexedMesh(GIuint start, GIuint end, 
							  GIsizei count, const GIuint *indices)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIuint *pSubset[GI_SUBSET_COUNT];
	GIsizei iNumVertices = indices ? ((GIsizei)end-(GIsizei)start+1) : count;
	GIint i, j, a;
	GIuint idx, bidx = start;
	GIuint uiPosAttrib = pContext->semantic[GI_POSITION_ATTRIB-GI_SEMANTIC_BASE];
	GIuint uiParamAttrib = pContext->semantic[GI_PARAM_ATTRIB-GI_SEMANTIC_BASE];
	GIuint uiStretchAttrib = pContext->semantic[GI_PARAM_STRETCH_ATTRIB-GI_SEMANTIC_BASE];
	GIdouble dNormSqr;
	GIFace *pFace;
	GIEdge *pEdge;
	GIHalfEdge *pHalfEdge, *pHBoundary, hedge[3];
	GIVertex *pVertex, *pVertex1, *pVertex2, *pRoot = NULL;
	GIAttribute *pAttribute;
	GIParam *pParam = NULL;
	GIVertex *pCorners[3];
	GIUIntPair pair;
	GIHash hVectorVertexMap, hVertexEdgeMap, hAttribMap;
	GIVertex **pIndexVertexMap;
	GIAttribute **pIndexAttributeMap;
	GIuint *pFaceCount;
	GIvoid *pKey;
	GIfloat *pPackedAttribs;
	GIboolean bAttributes = GI_FALSE, bParams = pContext->attrib_enabled[uiParamAttrib];
	GIboolean bManifold = GI_TRUE;

	N_Cylinder errorEdge;
	FILE* fp = fopen("NonManifold.obj", "w");
	/* error checking */
	count -= count % 3;
	if(!pMesh || !pContext->attrib_enabled[uiPosAttrib] || 
		pContext->attrib_size[uiPosAttrib] < 3 || 
		(pContext->attrib_enabled[uiParamAttrib] && 
		pContext->attrib_size[uiParamAttrib] < 2) || 
		pContext->attrib_enabled[uiStretchAttrib] || 
		!count || iNumVertices < 3)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	GIDebug(printf("create mesh ... "));

	/* initialize mesh */
	GIMesh_destruct(pMesh);
	GIDynamicQueue_construct(&pMesh->split_hedges);
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
	{
		pMesh->asemantic[a] = pContext->attrib_semantic[a];
		if(a == uiPosAttrib)
			pMesh->asize[a] = 3;
		else if(a == uiParamAttrib)
			pMesh->asize[a] = 2;
		else if(a == uiStretchAttrib)
			pMesh->asize[a] = 1;
		else if(pContext->attrib_enabled[a])
		{
			pMesh->aoffset[a] = sizeof(GIAttribute) + pMesh->attrib_size;
			pMesh->asize[a] = pContext->attrib_size[a];
			pMesh->anorm[a] = pContext->attrib_normalized[a];
			pMesh->attrib_size += pMesh->asize[a] * sizeof(GIfloat);
			bAttributes = GI_TRUE;
		}
	}
	memcpy(pMesh->semantic, pContext->semantic, GI_SEMANTIC_COUNT*sizeof(GIuint));

	/* initialize subsets */
	for(a=0; a<GI_SUBSET_COUNT; ++a)
	{
		if(pContext->subset_enabled[a] && pContext->subset[a])
		{
			if(pContext->subset_sorted[a])
				pSubset[a] = (GIuint *)pContext->subset[a];
			else
			{
				GIsizei iCount = pContext->subset_count[a];
				pSubset[a] = (GIuint*)GI_MALLOC_ARRAY(iCount, sizeof(GIuint));
				memcpy(pSubset[a], pContext->subset[a], iCount*sizeof(GIuint));
				qsort(pSubset[a], iCount, sizeof(GIuint), compare);
			}
		}
		else
			pSubset[a] = NULL;
	}

	/* create hash tables and other maps */
	GIHash_construct(&hVectorVertexMap, iNumVertices, 0.0f, 
		3*sizeof(GIfloat), hash_vec3f, compare_vec3f, copy_vec3f);
	GIHash_construct(&hVertexEdgeMap, 3*iNumVertices, 0.0f, 
		sizeof(GIUIntPair), hash_uintpair, compare_uintpair, copy_uintpair);
	if(indices)
		pIndexVertexMap = (GIVertex**)GI_CALLOC_ARRAY(
			iNumVertices, sizeof(GIVertex*));
	if(bAttributes)
	{
		if(indices)
			pIndexAttributeMap = (GIAttribute**)GI_CALLOC_ARRAY(
				iNumVertices, sizeof(GIAttribute*));
		GIHash_construct(&hAttribMap, iNumVertices, 0.0f, sizeof(GIuint)+
			pMesh->attrib_size, hash_attribs, compare_attribs, copy_attribs);
		pKey = GI_MALLOC_SINGLE(sizeof(GIuint)+pMesh->attrib_size);
		*((GIuint*)pKey) = pMesh->attrib_size;
		pPackedAttribs = (GIfloat*)((GIuint*)pKey+1);
	}
	pFaceCount = (GIuint*)GI_CALLOC_ARRAY(iNumVertices, sizeof(GIuint));

	/* create faces (and other data respectively) */
	GI_VEC3_SET(pMesh->aabb_min, DBL_MAX, DBL_MAX, DBL_MAX);
	GI_VEC3_SET(pMesh->aabb_max, -DBL_MAX, -DBL_MAX, -DBL_MAX);
	pMesh->radius = 0.0;
	for(i=0; i<count; i+=3)
	{
		pFace = (GIFace*)GI_MALLOC_PERSISTENT(sizeof(GIFace));
		GI_LIST_ADD(pMesh->faces, pFace);
		pFace->id = pMesh->fcount++;
		pFace->hedges = NULL;
		for(j=0; j<3; ++j, ++bidx)
		{
			/* create half edge */
			pHalfEdge = hedge + j;
			GI_LIST_ADD(pFace->hedges, pHalfEdge);
			pHalfEdge->face = pFace;
			pHalfEdge->edge = NULL;
			pHalfEdge->twin = NULL;

			/* vertex index already visited? */
			if(indices)
			{
				bidx = indices[i+j];
				idx = bidx - start;
				pVertex = pIndexVertexMap[idx];
			}
			if(!pVertex || !indices)
			{
				/* vertex with same coordinates already existing? */
				const GIfloat *fvec = pContext->attrib_pointer[uiPosAttrib] + 
					bidx*pContext->attrib_stride[uiPosAttrib];
				pVertex = (GIVertex *)GIHash_find(&hVectorVertexMap, fvec);
				if(!pVertex)
				{
					/* create new vertex */
					pVertex = (GIVertex*)GI_MALLOC_PERSISTENT(sizeof(GIVertex));
					GI_LIST_ADD(pMesh->vertices, pVertex);
					pVertex->id = pMesh->vcount++;
					pVertex->hedge = NULL;
					pVertex->flags = 0;
					pVertex->cut_degree = 0;
					GI_VEC3_COPY(pVertex->coords, fvec);
					GIHash_insert(&hVectorVertexMap, fvec, pVertex);

					/* check subsets */
					for(a=0; a<GI_SUBSET_COUNT; ++a)
						if(pSubset[a] && bsearch(&bidx, pSubset[a], 
							pContext->subset_count[a], sizeof(GIuint), compare))
							pVertex->flags |= 1 << a;

					/* analyse geometry */
					GI_VEC3_MIN(pMesh->aabb_min, pMesh->aabb_min, pVertex->coords);
					GI_VEC3_MAX(pMesh->aabb_max, pMesh->aabb_max, pVertex->coords);
					dNormSqr = GI_VEC3_LENGTH_SQR(pVertex->coords);
					if(dNormSqr > pMesh->radius)
						pMesh->radius = dNormSqr;
				}
				if(indices)
					pIndexVertexMap[idx] = pVertex;
			}
			pHalfEdge->vstart = pVertex;
			pCorners[j] = pVertex;
			++pFaceCount[pVertex->id];

			/* attribute index already visited? */
			if(bAttributes)
			{
				if(indices)
					pAttribute = pIndexAttributeMap[idx];
				if(!pAttribute || !indices)
				{
					for(a=0; a<GI_ATTRIB_COUNT; ++a)
					{
						if(pMesh->aoffset[a] >= 0)
						{
							memcpy(pPackedAttribs, pContext->attrib_pointer[a]+
								bidx*pContext->attrib_stride[a], 
								pContext->attrib_size[a]*sizeof(GIfloat));
							pPackedAttribs += pContext->attrib_size[a];
						}
					}
					pPackedAttribs = (GIfloat*)((GIuint*)pKey+1);
					pAttribute = (GIAttribute*)GIHash_find(&hAttribMap, pKey);
					if(!pAttribute)
					{
						pAttribute = (GIAttribute*)GI_MALLOC_PERSISTENT(
							sizeof(GIAttribute)+pMesh->attrib_size);
						GI_LIST_ADD(pMesh->attributes, pAttribute);
						pAttribute->id = pMesh->acount++;
						memcpy((GIbyte*)pAttribute+sizeof(GIAttribute), 
							pPackedAttribs, pMesh->attrib_size);
						GIHash_insert(&hAttribMap, pKey, pAttribute);
					}
					if(indices)
						pIndexAttributeMap[idx] = pAttribute;
				}
				pHalfEdge->astart = pAttribute;
			}
			else
				pHalfEdge->astart = NULL;

			//set params
			pHalfEdge->pstart = bParams ? ((GIParam*)(
				pContext->attrib_pointer[uiParamAttrib]+bidx*
				pContext->attrib_stride[uiParamAttrib])) : NULL;
		}

		/* create edges */
		for(j=0; j<3; ++j)
		{
			/* vertex pair already visited? */
			pVertex1 = pCorners[j];
			pVertex2 = pCorners[(j+1)%3];
			pair.first = pVertex1->id;
			pair.second = pVertex2->id;
			pEdge = (GIEdge*)GIHash_find(&hVertexEdgeMap, &pair);
			if(pEdge)
			{
				if (pEdge->hedge[0].twin)
				{
					double d = (pMesh->aabb_max[0] - pMesh->aabb_min[0])*0.01;
					int color[3] = { 255,0,0 };
					double p1[3] = { (double)pVertex1->coords[0] + d, (double)pVertex1->coords[1], (double)pVertex1->coords[2] };
					double p2[3] = { (double)pVertex2->coords[0], (double)pVertex2->coords[1] - d, (double)pVertex2->coords[2] };
					Cylinder cy(p1, p2, 0.5*d);
					cy.setColor(color);
					errorEdge.Add(cy);

					printf("v %f %f %f\n", (double)pVertex1->coords[0], (double)pVertex1->coords[1], (double)pVertex1->coords[2]);
					printf("#id %d\n", pVertex1->id);
					printf("v %f %f %f\n", (double)pVertex2->coords[0], (double)pVertex2->coords[1], (double)pVertex2->coords[2]);
					printf("#id %d\n", pVertex2->id);

					fprintf(fp, "v %f %f %f 255 0 0\n", (double)pVertex1->coords[0], (double)pVertex1->coords[1], (double)pVertex1->coords[2]);
					fprintf(fp, "#id %d\n", pVertex1->id);
					fprintf(fp, "v %f %f %f 255 0 0\n", (double)pVertex2->coords[0], (double)pVertex2->coords[1], (double)pVertex2->coords[2]);
					fprintf(fp, "#id %d\n", pVertex2->id);
					bManifold = GI_FALSE;
					//break;
				}
				pHalfEdge = &pEdge->hedge[1];
				*pHalfEdge = hedge[j];
				pHalfEdge->twin = &pEdge->hedge[0];
				pHalfEdge->twin->twin = pHalfEdge;
			}
			else
			{
				pEdge = (GIEdge*)GI_MALLOC_PERSISTENT(sizeof(GIEdge));
				GI_LIST_ADD(pMesh->edges, pEdge);
				pEdge->id = pMesh->ecount++;
				pHalfEdge = &pEdge->hedge[0];
				*pHalfEdge = hedge[j];
				pEdge->hedge[1].edge = pEdge;
				pEdge->length = GIvec3d_dist(pHalfEdge->vstart->coords, 
					pHalfEdge->next->vstart->coords);
				pMesh->mean_edge += pEdge->length;
				GIHash_insert(&hVertexEdgeMap, &pair, pEdge);
			}
			pHalfEdge->edge = pEdge;
			if(!pHalfEdge->vstart->hedge)
				pHalfEdge->vstart->hedge = pHalfEdge;
			GI_LIST_REMOVE(pFace->hedges, &hedge[j]);
			GI_LIST_ADD(pFace->hedges, pHalfEdge);
		}
	}
	pMesh->radius = sqrt(pMesh->radius);
	pMesh->mean_edge /= (GIdouble)pMesh->ecount;

	/* complete mesh */
	if(bManifold)
	{
		/* create boundary halfedges */
		GI_LIST_FOREACH(pMesh->edges, pEdge)
			if(!pEdge->hedge[0].twin)
			{
				pHalfEdge = &pEdge->hedge[1];
				pHBoundary = NULL;
				while(pHalfEdge != pHBoundary)
				{
					/* create halfedge and find next one */
					GI_LIST_ADD(pHBoundary, pHalfEdge);
					pHalfEdge->face = NULL;
					pHalfEdge->astart = NULL;
					pHalfEdge->pstart = NULL;
					pHalfEdge->twin = &pHalfEdge->edge->hedge[0];
					pHalfEdge->twin->twin = pHalfEdge;
					pHalfEdge->vstart = pHalfEdge->twin->next->vstart;
					pHalfEdge->vstart->hedge = pHalfEdge;
					for(pHalfEdge=pHalfEdge->twin->prev; 
						pHalfEdge->twin && pHalfEdge->edge!=pHBoundary->edge; 
						pHalfEdge=pHalfEdge->twin->prev) ;
					pHalfEdge = &pHalfEdge->edge->hedge[1];
				}
			}
		GI_LIST_NEXT(pMesh->edges, pEdge)

		/* check if manifold */
		GI_LIST_FOREACH(pMesh->vertices, pVertex)
			i = pFaceCount[pVertex->id];
			pHalfEdge = pVertex->hedge;
			if(!pHalfEdge->face)
				pHalfEdge = pHalfEdge->twin->next;
			do
			{
				--i;
				pHalfEdge = pHalfEdge->twin->next;
			}while(pHalfEdge != pVertex->hedge);
			if(i)
			{
				double d = (pMesh->aabb_max[0] - pMesh->aabb_min[0])*0.01;
				int color[3] = { 255,0,0 };
				double p1[3] = { (double)pVertex->coords[0]-d, (double)pVertex->coords[1], (double)pVertex->coords[2] };
				double p2[3] = { (double)pVertex->hedge->vstart->coords[0], (double)pVertex->hedge->vstart->coords[1] + d, (double)pVertex->hedge->vstart->coords[2] };
				Cylinder cy(p1, p2, 0.5*d);
				cy.setColor(color);
				errorEdge.Add(cy);

				printf("v %f %f %f\n", (double)pVertex->coords[0], (double)pVertex->coords[1], (double)pVertex->coords[2]);
				printf("#id %d\n", pVertex->id);
				printf("v %f %f %f\n", (double)pVertex->hedge->vstart->coords[0], (double)pVertex->hedge->vstart->coords[1], (double)pVertex->hedge->vstart->coords[2]);
				printf("#id %d\n", pVertex->hedge->vstart->id);

				fprintf(fp, "v %f %f %f 255 0 0\n", (double)pVertex->coords[0], (double)pVertex->coords[1], (double)pVertex->coords[2]);
				fprintf(fp, "#id %d\n", pVertex->id);
				fprintf(fp, "v %f %f %f 255 0 0\n", (double)pVertex->hedge->vstart->coords[0], (double)pVertex->hedge->vstart->coords[1], (double)pVertex->hedge->vstart->coords[2]);
				fprintf(fp, "#id %d\n", pVertex->hedge->vstart->id);
				printf("Manifold Error.\n");
				bManifold = GI_FALSE;
				//break;
			}
		GI_LIST_NEXT(pMesh->vertices, pVertex)
	}

	/* clean up subsets */
	for(a=0; a<GI_SUBSET_COUNT; ++a)
		if(pSubset[a] && !pContext->subset_sorted[a])
			GI_FREE_ARRAY(pSubset[a]);

	/* clean up */
	GIHash_destruct(&hVectorVertexMap, 0);
	GIHash_destruct(&hVertexEdgeMap, 0);
	if(indices)
		GI_FREE_ARRAY(pIndexVertexMap);
	if(bAttributes)
	{
		if(indices)
			GI_FREE_ARRAY(pIndexAttributeMap);
		GIHash_destruct(&hAttribMap, 0);
		GI_FREE_SINGLE(pKey, sizeof(GIuint)+pMesh->attrib_size);
	}
	GI_FREE_ARRAY(pFaceCount);
	GIDebug(printf("finished\n"));
	fclose(fp);
	if (bManifold) _unlink("NonManifold.obj");
	else{
		errorEdge.put("NonManifold_edge.obj", (pMesh->aabb_max[0] - pMesh->aabb_min[0])*0.01);

		//int iNumIndices = 0;
		//int iNumVertices = 0;
		//unsigned int *pIndices = NULL;
		//GLfloat *pVertices = NULL;
		//GLfloat *pUV = NULL;
		//GLfloat *pNormals = NULL;

		//putMesh("aaa", &iNumIndices, &iNumVertices, &pIndices, &pVertices, &pUV, &pNormals);
	}

	/* check for manifold and set params */
	if(!bManifold || (bParams && 
		!GICutter_from_params(&pContext->cutter, pMesh)))
	{
		GIMesh_destruct(pMesh);
		GIContext_error(pContext, GI_INVALID_MESH);
	}
}

/** Create mesh from current attribute arrays.
 *  This function creates a new mesh by walking through the currently set and 
 *  enabled attribute arrays. WARNING: All previous mesh data of the current 
 *  bound mesh object will be deleted.
 *  \param first index to start at
 *  \param count number of vertices to process
 *  \ingroup mesh
 */
void GIAPIENTRY giNonIndexedMesh(GIint first, GIsizei count)
{
	giIndexedMesh(first, first+count-1, count, NULL);
}

/** Copy mesh from existing mesh
 *  \param mesh mesh to copy from
 *  \ingroup mesh
 */
void GIAPIENTRY giCopyMesh(GIuint mesh)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh, *pSource;
	GIPatch *pSPatch, *pDPatch = NULL;
	GICutPath *pSPath, *pDPath;
	GIFace *pSFace, *pDFace;
	GIEdge *pSEdge, *pDEdge;
	GIHalfEdge *pSHalfEdge, *pDHalfEdge, *pHBoundary;
	GIVertex *pSVertex, *pDVertex;
	GIAttribute *pSAttribute, *pDAttribute;
	GIParam *pSParam, *pDParam;
	GIPatch **pIndexPatchMap;
	GIFace **pIndexFaceMap;
	GIEdge **pIndexEdgeMap;
	GIVertex **pIndexVertexMap;
	GIAttribute **pIndexAttributeMap;
	GIParam **pIndexParamMap;
	GICutPath ***pIndexPathMap;
	GIQueueNode *pQNode;
	GISplitInfo *pSSplit, *pDSplit;
	GIuint i, a, uiFaces = 0, uiNextPatch = 0;

	/* error checking */
	if(!pMesh)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	pSource = (GIMesh*)GIHash_find(&pContext->mesh_hash, &mesh);
	if(!pSource)
	{
		GIContext_error(pContext, GI_INVALID_ID);
		return;
	}
	pSPatch = pSource->patches;
	GIMesh_destruct(pMesh);

	/* copy data */
	pMesh->fcount = pSource->fcount;
	pMesh->ecount = pSource->ecount;
	pMesh->vcount = pSource->vcount;
	pMesh->acount = pSource->acount;
	pMesh->pcount = pSource->pcount;
	pMesh->attrib_size = pSource->attrib_size;
	pMesh->varray_size = pSource->varray_size;
	pMesh->varray_attribs = pSource->varray_attribs;
	pMesh->genus = pSource->genus;
	GI_VEC3_COPY(pMesh->aabb_min, pSource->aabb_min);
	GI_VEC3_COPY(pMesh->aabb_max, pSource->aabb_max);
	pMesh->radius = pSource->radius;
	pMesh->mean_edge = pSource->mean_edge;
	pMesh->surface_area = pSource->surface_area;
	pMesh->param_area = pSource->param_area;
	memcpy(pMesh->stretch, pSource->stretch, GI_STRETCH_COUNT*sizeof(GIdouble));
	pMesh->min_param_stretch = pSource->min_param_stretch;
	pMesh->max_param_stretch = pSource->max_param_stretch;
	pMesh->param_metric = pSource->param_metric;
	GIDynamicQueue_construct(&pMesh->split_hedges);
	pMesh->patch_count = pSource->patch_count;
	pMesh->param_patches = pSource->param_patches;
	pMesh->resolution = pSource->resolution;
	pMesh->cut_splits = pSource->cut_splits;
	pMesh->active_patch = NULL;
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
	{
		pMesh->aoffset[a] = pSource->aoffset[a];
		pMesh->asize[a] = pSource->asize[a];
		pMesh->anorm[a] = pSource->anorm[a];
		pMesh->asemantic[a] = pSource->asemantic[a];
	}
	memcpy(pMesh->semantic, pSource->semantic, GI_SEMANTIC_COUNT*sizeof(GIuint));

	/* auxiliary datastructures */
	if(pSPatch)
	{
		uiNextPatch = pSPatch->fcount;
		pIndexPatchMap = (GIPatch**)GI_MALLOC_ARRAY(pSource->patch_count, sizeof(GIPatch*));
		pIndexPathMap = (GICutPath***)GI_MALLOC_ARRAY(pSource->patch_count, sizeof(GICutPath**));
		pMesh->patches = (GIPatch*)GI_MALLOC_ARRAY(pSource->patch_count, sizeof(GIPatch));
		pMesh->varray_patch = pSource->varray_patch ? 
			(pMesh->patches+pSource->varray_patch->id) : NULL;
	}
	else
		pMesh->patches = pMesh->varray_patch = NULL;
	pIndexFaceMap = (GIFace**)GI_MALLOC_ARRAY(pSource->fcount, sizeof(GIFace*));
	pIndexEdgeMap = (GIEdge**)GI_CALLOC_ARRAY(pSource->ecount, sizeof(GIEdge*));
	pIndexVertexMap = (GIVertex**)GI_CALLOC_ARRAY(pSource->ecount, sizeof(GIVertex*));
	if(pSource->attributes)
		pIndexAttributeMap = (GIAttribute**)GI_CALLOC_ARRAY(pSource->acount, sizeof(GIAttribute*));

	/* copy faces (and other data respectively) */
	GI_LIST_FOREACH(pSource->faces, pSFace)
		pDFace = (GIFace*)GI_MALLOC_PERSISTENT(sizeof(GIFace));
		GI_LIST_ADD(pMesh->faces, pDFace);
		pDFace->id = pSFace->id;
		pDFace->hedges = NULL;
		pIndexFaceMap[pSFace->id] = pDFace;

		/* create new patch and copy data if neccessary */
		if(pSPatch && !pDPatch)
		{
			pDPatch = pMesh->patches + pSPatch->id;
			pDPatch->id = pSPatch->id;
			pDPatch->mesh = pMesh;
			pDPatch->fcount = pSPatch->fcount;
			pDPatch->pcount = pSPatch->pcount;
			pDPatch->hcount = pSPatch->hcount;
			pDPatch->path_count = pSPatch->path_count;
			pDPatch->faces = pDFace;
			pDPatch->params = NULL;
			pDPatch->paths = NULL;
			pDPatch->groups = pSPatch->groups;
			pDPatch->hlength = pSPatch->hlength;
			memcpy(pDPatch->side_lengths, pSPatch->side_lengths, 4*sizeof(GIdouble));
			pDPatch->resolution = pSPatch->resolution;
			pDPatch->surface_area = pSPatch->surface_area;
			pDPatch->param_area = pSPatch->param_area;
			memcpy(pDPatch->stretch, pSPatch->stretch, GI_STRETCH_COUNT*sizeof(GIdouble));
			pDPatch->min_param_stretch = pSPatch->min_param_stretch;
			pDPatch->max_param_stretch = pSPatch->max_param_stretch;
			pDPatch->param_metric = pSPatch->param_metric;
			pDPatch->parameterized = pSPatch->parameterized;
			pDPatch->fixed_corners = pSPatch->fixed_corners;
			GIDynamicQueue_construct(&pDPatch->split_paths);
			pDPatch->next = pMesh->patches + ((pDPatch->id+1)%pMesh->patch_count);
			pIndexPatchMap[pSPatch->id] = pDPatch;
			pIndexParamMap = (GIParam**)GI_CALLOC_ARRAY(pSPatch->pcount, sizeof(GIParam*));
			pIndexPathMap[pDPatch->id] = (GICutPath**)GI_CALLOC_ARRAY(pSPatch->path_count, sizeof(GICutPath*));
		}

		/* create half edges and related data */
		GI_LIST_FOREACH(pSFace->hedges, pSHalfEdge)
			/* copy edge and connect to half edge */
			pSEdge = pSHalfEdge->edge;
			pDEdge = pIndexEdgeMap[pSEdge->id];
			if(!pDEdge)
			{
				pDEdge = (GIEdge*)GI_MALLOC_PERSISTENT(sizeof(GIEdge));
				GI_LIST_ADD(pMesh->edges, pDEdge);
				pDEdge->id = pSEdge->id;
				pDEdge->length = pSEdge->length;
				pDHalfEdge = &pDEdge->hedge[0];
				pDEdge->hedge[1].edge = pDEdge;
				pIndexEdgeMap[pSEdge->id] = pDEdge;
			}
			else
			{
				pDHalfEdge = &pDEdge->hedge[1];
				pDHalfEdge->twin = &pDEdge->hedge[0];
				pDHalfEdge->twin->twin = pDHalfEdge;
			}
			pDHalfEdge->edge = pDEdge;
			pDHalfEdge->face = pDFace;
			GI_LIST_ADD(pDFace->hedges, pDHalfEdge);

			/* copy vertex and connect to half edge */
			pSVertex = pSHalfEdge->vstart;
			pDVertex = pIndexVertexMap[pSVertex->id];
			if(!pDVertex)
			{
				pDVertex = (GIVertex*)GI_MALLOC_PERSISTENT(sizeof(GIVertex));
				GI_LIST_ADD(pMesh->vertices, pDVertex);
				pDVertex->id = pSVertex->id;
				GI_VEC3_COPY(pDVertex->coords, pSVertex->coords);
				pDVertex->flags = pSVertex->flags;
				pDVertex->cut_degree = pSVertex->cut_degree;
				pDVertex->hedge = pDHalfEdge;
				pIndexVertexMap[pSVertex->id] = pDVertex;
			}
			pDHalfEdge->vstart = pDVertex;

			/* copy attribute and connect to half edge */
			if(pSource->attributes)
			{
				pSAttribute = pSHalfEdge->astart;
				pDAttribute = pIndexAttributeMap[pSAttribute->id];
				if(!pDAttribute)
				{
					pDAttribute = (GIAttribute*)GI_MALLOC_PERSISTENT(sizeof(GIAttribute)+pSource->attrib_size);
					GI_LIST_ADD(pMesh->attributes, pDAttribute);
					pDAttribute->id = pSAttribute->id;
					memcpy((GIbyte*)pDAttribute+sizeof(GIAttribute), 
						(GIbyte*)pSAttribute+sizeof(GIAttribute), pSource->attrib_size);
					pIndexAttributeMap[pSAttribute->id] = pDAttribute;
				}
				pDHalfEdge->astart = pDAttribute;
			}
			else
				pDHalfEdge->astart = NULL;

			/* copy parameter coordinates and connect to half edge */
			if(pDPatch)
			{
				pSParam = pSHalfEdge->pstart;
				pDParam = pIndexParamMap[pSParam->id];
				if(!pDParam)
				{
					pDParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
					GI_LIST_ADD(pDPatch->params, pDParam);
					pDParam->id = pSParam->id;
					GI_VEC2_COPY(pDParam->params, pSParam->params);
					pDParam->stretch = pSParam->stretch;
					pDParam->vertex = pDVertex;
					pDParam->cut_hedge = NULL;
					pIndexParamMap[pSParam->id] = pDParam;
				}
				pDHalfEdge->pstart = pDParam;
				if(pSParam->cut_hedge == pSHalfEdge)
					pDParam->cut_hedge = pDHalfEdge;
			}
			else
				pDHalfEdge->pstart = NULL;
		GI_LIST_NEXT(pSFace->hedges, pSHalfEdge)

		/* finish patch and create paths */
		if(++uiFaces == uiNextPatch)
		{
			/* set root and copy corners */
			pDPatch->params = pIndexParamMap[pSPatch->params->id];
			if(pSPatch->corners[0])
				for(i=0; i<4; ++i)
					pDPatch->corners[i] = pIndexParamMap[pSPatch->corners[i]->id];
			else
				memset(pDPatch->corners, 0, 4*sizeof(GIParam*));

			/* create paths */
			i = pDPatch->id;
			GI_LIST_FOREACH(pSPatch->paths, pSPath)
				pDPath = (GICutPath*)GI_MALLOC_PERSISTENT(sizeof(GICutPath));
				GI_LIST_ADD(pDPatch->paths, pDPath);
				pDPath->id = pSPath->id;
				pDPath->patch = pDPatch;
				pDPath->group = pSPath->group;
				pDPath->elength = pSPath->elength;
				pDPath->glength = pSPath->glength;
				pDPath->pstart = pIndexParamMap[pSPath->pstart->id];
				pDPath->twin = NULL;
				pIndexPathMap[i][pSPath->id] = pDPath;
			GI_LIST_NEXT(pSPatch->paths, pSPath)

			/* copy split stack */
			pQNode = pSPatch->split_paths.head;
			while(pQNode)
			{
				GIDynamicQueue_enqueue(&pDPatch->split_paths, 
					pIndexPathMap[i][((GICutPath*)pQNode->data)->id]);
				pQNode = pQNode->next;
			}

			/* next patch */
			pDPatch = NULL;
			pSPatch = pSPatch->next;
			uiNextPatch += pSPatch->fcount;
			GI_FREE_ARRAY(pIndexParamMap);
		}
	GI_LIST_NEXT(pSource->faces, pSFace)

	/* copy boundaries */
	GI_LIST_FOREACH(pSource->edges, pSEdge)
		if(!pSEdge->hedge[1].face)
		{
			pHBoundary = NULL;
			pSHalfEdge = &pSEdge->hedge[1];
			pDHalfEdge = &pIndexEdgeMap[pSEdge->id]->hedge[1];
			while(pDHalfEdge != pHBoundary)
			{
				/* create halfedge and find next one */
				GI_LIST_ADD(pHBoundary, pDHalfEdge);
				pDHalfEdge->face = NULL;
				pDHalfEdge->astart = NULL;
				pDHalfEdge->pstart = NULL;
				pDHalfEdge->twin = &pDHalfEdge->edge->hedge[0];
				pDHalfEdge->twin->twin = pDHalfEdge;
				pDHalfEdge->vstart = pDHalfEdge->twin->next->vstart;
				pDHalfEdge->vstart->hedge = pDHalfEdge;
				pSHalfEdge = pSHalfEdge->next;
				pDHalfEdge = &pIndexEdgeMap[pSHalfEdge->edge->id]->hedge[1];
			}
		}
	GI_LIST_NEXT(pSource->edges, pSEdge)

	/* copy split stack */
	pQNode = pSource->split_hedges.head;
	while(pQNode)
	{
		pSSplit = (GISplitInfo*)pQNode->data;
		pDSplit = (GISplitInfo*)GI_MALLOC_PERSISTENT(sizeof(GISplitInfo));
		pDSplit->hedge = GIFace_halfedge_at(pIndexFaceMap[
			pSSplit->hedge->face->id], GIHalfEdge_index(pSSplit->hedge));
		if(pSSplit->patch)
			pDSplit->patch = pIndexPatchMap[pSSplit->patch->id];
		if(pSSplit->twin_patch)
			pDSplit->twin_patch = pIndexPatchMap[pSSplit->twin_patch->id];
		GIDynamicQueue_enqueue(&pMesh->split_hedges, pDSplit);
		pQNode = pQNode->next;
	}

	/* connect path twins */
	for(i=0; i<pMesh->patch_count; ++i)
	{
		pSPatch = pSource->patches + i;
		pDPatch = pMesh->patches + i;
		GI_LIST_FOREACH(pSPatch->paths, pSPath)
			if(pSPath->twin)
				pIndexPathMap[i][pSPath->id]->twin = 
					pIndexPathMap[pSPath->twin->patch->id][pSPath->twin->id];
		GI_LIST_NEXT(pSPatch->paths, pSPath)
	}

	/* clean up */
	GI_FREE_ARRAY(pIndexFaceMap);
	GI_FREE_ARRAY(pIndexEdgeMap);
	GI_FREE_ARRAY(pIndexVertexMap);
	if(pSource->attributes)
		GI_FREE_ARRAY(pIndexAttributeMap);
	if(pSource->patches)
	{
		for(i=0; i<pMesh->patch_count; ++i)
			GI_FREE_ARRAY(pIndexPathMap[i]);
		GI_FREE_ARRAY(pIndexPatchMap);
		GI_FREE_ARRAY(pIndexPathMap);
	}
}

/** Recompute per param stretch values.
 *  \ingroup mesh
 */
void GIAPIENTRY giComputeParamStretch(GIenum metric)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIPatch *pPatch;

	/* error checking */
	if(!pMesh || !pMesh->patches)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	if(metric < GI_STRETCH_BASE || metric > GI_STRETCH_END)
	{
		GIContext_error(pContext, GI_INVALID_ENUM);
		return;
	}
	pPatch = pMesh->active_patch;

	/* compute if neccessary */
	if(pPatch && pPatch->parameterized)
	{
		if(pPatch->param_metric != metric)
		{
			GIPatch_compute_stretch(pPatch, metric, GI_TRUE, GI_TRUE);
			pMesh->param_metric = pMesh->patch_count>1 ? 0 : pPatch->param_metric;
		}
	}
	else if(!pPatch && pMesh->param_patches == pMesh->patch_count)
	{
		if(pMesh->param_metric != metric)
			GIMesh_compute_stretch(pMesh, metric, GI_TRUE);
	}
	else
		GIContext_error(pContext, GI_INVALID_OPERATION);
}

/** Set active patch.
 *  \param patch patch to set active
 *  \ingroup mesh
 */
void GIAPIENTRY giMeshActivePatch(GIint patch)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;

	/* error checking */
	if(!pMesh || !pMesh->patches)
		GIContext_error(pContext, GI_INVALID_OPERATION);
	else if(patch < 0)
		pMesh->active_patch = NULL;
	else if(patch > pMesh->patch_count)
		GIContext_error(pContext, GI_INVALID_VALUE);
	else
		pMesh->active_patch = pMesh->patches + patch;
}

/** Retrieve boolean mesh or patch property.
 *  \param pname property to query
 *  \param params address of variable to store value
 *  \ingroup mesh
 */
void GIAPIENTRY giGetMeshbv(GIenum pname, GIboolean *params)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIPatch *pPatch;

	/* error checking */
	if(!pMesh)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	pPatch = pMesh->patch_count>1 ? pMesh->active_patch : NULL;

	/* select property and get value */
	switch(pname)
	{
	case GI_HAS_PARAMS:
		*params = pPatch ? pPatch->parameterized : (pMesh->param_patches>0);
		break;
	case GI_HAS_CUT:
		*params = pMesh->patches != NULL;
		break;
	default:
		GIContext_error(pContext, GI_INVALID_ENUM);
	}
}

/** Retrieve integer mesh or patch property.
 *  \param pname property to query
 *  \param params address of variable to store value
 *  \ingroup mesh
 */
void GIAPIENTRY giGetMeshiv(GIenum pname, GIint *params)
{
	/* error checking */
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIPatch *pPatch;
	if(!pMesh)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	pPatch = pMesh->patch_count>1 ? pMesh->active_patch : NULL;

	/* select property and get value */
	switch(pname)
	{
	case GI_PATCH_COUNT:
		if(pPatch)
			GIContext_error(pContext, GI_INVALID_ENUM);
		else
			*params = pMesh->patch_count;
		break;
	case GI_FACE_COUNT:
		*params = pPatch ? pPatch->fcount : pMesh->fcount;
		break;
	case GI_EDGE_COUNT:
		*params = pPatch ? (pPatch->fcount+pPatch->pcount-1) : pMesh->ecount;
		break;
	case GI_VERTEX_COUNT:
		*params = pPatch ? pPatch->pcount : pMesh->vcount;
		break;
	case GI_PARAM_COUNT:
		if(pPatch)
			*params = pPatch->pcount;
		else
		{
			GIuint i;
			*params = 0;
			if(pMesh->param_patches)
				for(i=0; i<pMesh->patch_count; ++i)
					*params += pMesh->patches[i].pcount;
		}
		break;
	case GI_PARAM_STRETCH_METRIC:
		*params = pPatch ? pPatch->param_metric : pMesh->param_metric;
		break;
	case GI_PARAM_RESOLUTION:
		*params = pPatch ? pPatch->resolution : pMesh->resolution;
		break;
	case GI_TOPOLOGICAL_SIDEBAND_LENGTH:
		if(pMesh->patch_count == 1)
			pPatch = pMesh->patches;
		if(pPatch)
			*params = pPatch->path_count;
		else
			GIContext_error(pContext, GI_INVALID_ENUM);
		break;
	case GI_TOPOLOGICAL_SIDEBAND:
		if(pMesh->patch_count == 1)
			pPatch = pMesh->patches;
		if(pPatch)
		{
			GICutPath *pPath;
			GIParam *pParam = pPatch->params;
			GIuint c = 0;
			GI_LIST_FOREACH(pPatch->paths, pPath)
				if(pMesh->patch_count > 1)
					params[c] = pPath->twin ? pPath->twin->patch->id : -1;
				else
					params[c] = pPath->group;
				do
				{
					++params[c+1];
					pParam = pParam->cut_hedge->next->pstart;
				}while(pParam != pPath->next->pstart);
				params[c+2] = pPath->glength;
				c += 3;
			GI_LIST_NEXT(pPatch->paths, pPath)
		}
		else
			GIContext_error(pContext, GI_INVALID_ENUM);
		break;
	case GI_ACTIVE_PATCH:
		*params = pMesh->active_patch ? pPatch->id : GI_ALL_PATCHES;
		break;
	case GI_POSITION_ATTRIB:
	case GI_PARAM_ATTRIB:
	case GI_PARAM_STRETCH_ATTRIB:
		*params = pMesh->semantic[pname-GI_SEMANTIC_BASE];
		break;
	default:
		{
			GIboolean bBool = GI_FALSE;
			giGetMeshbv(pname, &bBool);
			*params = bBool;
		}
	}
}

/** Retrieve floating point mesh property.
 *  \param pname property to query
 *  \param params address of variable to store value
 *  \ingroup mesh
 */
void GIAPIENTRY giGetMeshfv(GIenum pname, GIfloat *params)
{
	/* error checking */
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIPatch *pPatch;
	if(!pMesh)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	pPatch = pMesh->patch_count>1 ? pMesh->active_patch : NULL;

	/* select property and get value */
	switch(pname)
	{
	case GI_AABB_MIN:
		if(pPatch)
			GIContext_error(pContext, GI_INVALID_ENUM);
		else
			GI_VEC3_COPY(params, pMesh->aabb_min);
		break;
	case GI_AABB_MAX:
		if(pPatch)
			GIContext_error(pContext, GI_INVALID_ENUM);
		else
			GI_VEC3_COPY(params, pMesh->aabb_max);
		break;
	case GI_RADIUS:
		if(pPatch)
			GIContext_error(pContext, GI_INVALID_ENUM);
		else
			*params = pMesh->radius;
		break;
	case GI_MIN_PARAM_STRETCH:
		if((pPatch && !pPatch->param_metric) || (!pPatch && !pMesh->param_metric))
		{
			*params = 0.0;
			GIContext_error(pContext, GI_INVALID_OPERATION);
		}
		else if(pPatch)
			*params = pPatch->min_param_stretch;
		else
			*params = pMesh->min_param_stretch;
		break;
	case GI_MAX_PARAM_STRETCH:
		if((pPatch && !pPatch->param_metric) || (!pPatch && !pMesh->param_metric))
		{
			*params = 0.0;
			GIContext_error(pContext, GI_INVALID_OPERATION);
		}
		else if(pPatch)
			*params = pPatch->max_param_stretch;
		else
			*params = pMesh->max_param_stretch;
		break;
	case GI_MAX_GEOMETRIC_STRETCH:
	case GI_RMS_GEOMETRIC_STRETCH:
	case GI_COMBINED_STRETCH:
		if(pPatch && pPatch->parameterized)
		{
			if(!pPatch->stretch[pname-GI_STRETCH_BASE])
				GIPatch_compute_stretch(pPatch, pname, GI_FALSE, GI_FALSE);
			*params = pPatch->stretch[pname-GI_STRETCH_BASE] * 
				sqrt(pPatch->param_area/pPatch->surface_area);
		}
		else if(!pPatch && pMesh->param_patches)
		{
			if(!pMesh->stretch[pname-GI_STRETCH_BASE])
				GIMesh_compute_stretch(pMesh, pname, GI_FALSE);
			*params = pMesh->stretch[pname-GI_STRETCH_BASE] * 
				sqrt(pMesh->param_area/pMesh->surface_area);
		}
		else
			GIContext_error(pContext, GI_INVALID_OPERATION);
		break;
	default:
		{
			GIint iInt = 0;
			giGetMeshiv(pname, &iInt);
			*params = iInt;
		}
	}
}

/** Retrieve boolean mesh or patch attribute property.
 *  \param attrib attribute to query
 *  \param pname property to query
 *  \param params address of variable to store value
 *  \ingroup mesh
 */
void GIAPIENTRY giGetMeshAttribbv(GIuint attrib, GIenum pname, GIboolean *params)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIPatch *pPatch;

	/* error checking */
	if(!pMesh)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	else if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}
	pPatch = pMesh->patch_count>1 ? pMesh->active_patch : NULL;

	/* select property and get value */
	switch(pname)
	{
	case GI_HAS_ATTRIB:
		switch(pMesh->asemantic[attrib])
		{
		case GI_POSITION_ATTRIB:
			*params = GI_TRUE;
			break;
		case GI_PARAM_ATTRIB:
			*params = pPatch ? pPatch->parameterized : (pMesh->param_patches>0);
			break;
		case GI_PARAM_STRETCH_ATTRIB:
			*params = (pPatch && pPatch->param_metric) || 
				(!pPatch && pMesh->param_metric);
			break;
		default:
			*params = (pMesh->aoffset[attrib] >= 0);
		}
		break;
	case GI_ATTRIB_NORMALIZED:
		*params = pMesh->anorm[attrib];
		break;
	default:
		GIContext_error(pContext, GI_INVALID_ENUM);
	}
}

/** Retrieve integer mesh or patch attribute property.
 *  \param attrib attribute to query
 *  \param pname property to query
 *  \param param address of variable to store value
 *  \ingroup mesh
 */
void GIAPIENTRY giGetMeshAttribiv(GIuint attrib, GIenum pname, GIint *param)
{
	/* error checking */
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIPatch *pPatch;
	if(!pMesh)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	else if(attrib >= GI_ATTRIB_COUNT)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}
	pPatch = pMesh->patch_count>1 ? pMesh->active_patch : NULL;

	/* select property and get value */
	switch(pname)
	{
	case GI_ATTRIB_SIZE:
		*param = pMesh->asize[attrib];
		break;
	case GI_ATTRIB_SEMANTIC:
		*param = pMesh->asemantic[attrib];
		break;
	default:
		{
			GIboolean bBool = GI_FALSE;
			giGetMeshbv(pname, &bBool);
			*param = bBool;
		}
	}
}

/** Extract mesh data as indexed vertex arrays.
 *  \param vcount address to store number of vertices at
 *  \param icount address to store number of indices at
 *  \param indices array to fill or NULL if just querying sizes
 *  \ingroup mesh
 */
void GIAPIENTRY giGetIndexedMesh(GIuint *vcount, GIuint *icount, 
								 GIuint *indices)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIPatch *pPatch;
	GIFace *pFace, *pFEnd;
	GIHalfEdge *pHalfEdge;
	GIuint *pIndex = indices;
	GIuint uiAttribs = 0, uiNumAttribs = 0, uiAttrib[GI_ATTRIB_COUNT];
	GIfloat *dst, *fsrc;
	GIuint i, j, a, c, uiAttribSize = 0, idx, vidx = 0;
	GIboolean bParams, bStretch;
	GIvoid *pKey;
	GIfloat *pPackedAttribs;
	GIHash hAttribIndex;

	/* error checking */
	*vcount = *icount = 0;
	if(!pMesh || !pMesh->vertices)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	pPatch = pMesh->patch_count>1 ? pMesh->active_patch : NULL;

	/* one patch or whole mesh */
	if(pMesh->active_patch)
	{
		pFace = pMesh->active_patch->faces;
		pFEnd = pMesh->active_patch->next->faces;
		bParams = pMesh->active_patch->parameterized;
		bStretch = bParams && pMesh->active_patch->param_metric;
	}
	else
	{
		pFace = pFEnd = pMesh->faces;
		bParams = pMesh->param_patches > 0;
		bStretch = bParams && pMesh->param_metric;
	}

	/* error checking */
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
	{
		if(pContext->attrib_enabled[a])
		{
			if(pContext->attrib_semantic[a] != pMesh->asemantic[a] || 
				pContext->attrib_size[a] != pMesh->asize[a] || 
				(pMesh->asemantic[a] == GI_PARAM_ATTRIB && !bParams) || 
				(pMesh->asemantic[a] == GI_PARAM_STRETCH_ATTRIB && !bStretch) || 
				(pMesh->asemantic[a] == GI_NONE && pMesh->aoffset[a] < 0))
			{
				GIContext_error(pContext, GI_INVALID_OPERATION);
				return;
			}
			uiAttrib[uiNumAttribs++] = a;
			uiAttribSize += pMesh->asize[a] * sizeof(GIfloat);
			uiAttribs |= 1 << a;
		}
	}

	/* query size only */
	*icount = 3 * (pPatch ? pPatch->fcount : pMesh->fcount);
	if(!indices)
	{
		GIboolean bPos = uiAttribs & (1<<pMesh->semantic[
			GI_POSITION_ATTRIB-GI_SEMANTIC_BASE]);
		GIboolean bParam = uiAttribs & (1<<pMesh->semantic[GI_PARAM_ATTRIB-
			GI_SEMANTIC_BASE]) || uiAttribs & (1<<pMesh->semantic[
			GI_PARAM_STRETCH_ATTRIB-GI_SEMANTIC_BASE]);
		if(uiAttribs == pMesh->varray_attribs && pPatch == pMesh->varray_patch)
			*vcount = pMesh->varray_size;
		else if(uiNumAttribs == 1)
		{
			if(bPos)
				*vcount = pMesh->varray_size = pPatch ? pPatch->pcount : pMesh->vcount;
			else if(bParam && pPatch)
				*vcount = pMesh->varray_size = pPatch->pcount;
		}
		else if(uiNumAttribs == 2 && bPos)
		{
			if(bParam && pPatch)
				*vcount = pMesh->varray_size = pPatch->pcount;
			else if(!bParam && !pPatch)
				*vcount = pMesh->varray_size = pMesh->acount;
		}
	}
	pMesh->varray_attribs = uiAttribs;
	pMesh->varray_patch = pPatch;
	if(*vcount)
		return;

	/* create hash table */
	GIHash_construct(&hAttribIndex, pPatch ? pPatch->pcount : 
		pMesh->vcount, 0.0f, sizeof(GIuint)+uiAttribSize, 
		hash_attribs, compare_attribs, copy_attribs);
	pKey = GI_MALLOC_SINGLE(sizeof(GIuint)+uiAttribSize);
	*(GIuint*)pKey = uiAttribSize;

	/* fill arrays */
	do
	{
		for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,pHalfEdge=pHalfEdge->next)
		{
			/* pack attributes */
			pPackedAttribs = (GIfloat*)((GIuint*)pKey+1);
			for(j=0; j<uiNumAttribs; ++j)
			{
				a = uiAttrib[j];
				switch(pMesh->asemantic[a])
				{
				case GI_POSITION_ATTRIB:
					GI_VEC3_COPY(pPackedAttribs, pHalfEdge->vstart->coords);
					break;
				case GI_PARAM_ATTRIB:
					GI_VEC2_COPY(pPackedAttribs, pHalfEdge->pstart->params);
					break;
				case GI_PARAM_STRETCH_ATTRIB:
					*pPackedAttribs = pHalfEdge->pstart->stretch;
					break;
				default:
					fsrc = (GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]);
					for(c=0; c<pMesh->asize[a]; ++c)
						pPackedAttribs[c] = fsrc[c];
				}
				pPackedAttribs += pMesh->asize[a];
			}

			/* Compute index and copy data if neccessary */
			idx = (GIuint)GIHash_find(&hAttribIndex, pKey);
			if(!idx)
			{
				idx = ++*vcount;
				GIHash_insert(&hAttribIndex, pKey, (GIvoid*)idx);
			}
			if(indices)
			{
				*(pIndex++) = --idx;
				pPackedAttribs = (GIfloat*)((GIuint*)pKey+1);
				for(j=0; j<uiNumAttribs; ++j)
				{
					a = uiAttrib[j];
					dst = pContext->attrib_pointer[a] + 
						pContext->attrib_stride[a]*idx;
					for(c=0; c<pMesh->asize[a]; ++c)
						dst[c] = *(pPackedAttribs++);
				}
			}
		}
		pFace = pFace->next;
	}while(pFace != pFEnd);

	/* clean up */
	GIHash_destruct(&hAttribIndex, 0);
	GI_FREE_SINGLE(pKey, sizeof(GIuint)+uiAttribSize);
	pMesh->varray_size = *vcount;
}

/** Extract mesh data as non-indexed vertex arrays.
 *  \param vcount address to store number of vertices at
 *  \ingroup mesh
 */
void GIAPIENTRY giGetNonIndexedMesh(GIuint *vcount)
{
	GIContext *pContext = GIContext_current();
	GIMesh *pMesh = pContext->mesh;
	GIFace *pFace, *pFEnd;
	GIHalfEdge *pHalfEdge;
	GIuint uiNumAttribs = 0, uiAttrib[GI_ATTRIB_COUNT];
	GIfloat *dst, *fsrc;
	GIuint i, j, a, c, idx = 0;
	GIboolean bParams, bStretch;

	/* error checking */
	*vcount = 0;
	if(!pMesh || !pMesh->vertices)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}

	/* one patch or whole mesh */
	if(pMesh->active_patch)
	{
		pFace = pMesh->active_patch->faces;
		pFEnd = pMesh->active_patch->next->faces;
		bParams = pMesh->active_patch->parameterized;
		bStretch = bParams && pMesh->active_patch->param_metric;
	}
	else
	{
		pFace = pFEnd = pMesh->faces;
		bParams = pMesh->param_patches > 0;
		bStretch = bParams && pMesh->param_metric;
	}

	/* error checking */
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
	{
		if(pContext->attrib_enabled[a])
		{
			if(pContext->attrib_semantic[a] != pMesh->asemantic[a] || 
				pContext->attrib_size[a] != pMesh->asize[a] || 
				(pMesh->asemantic[a] == GI_PARAM_ATTRIB && !bParams) || 
				(pMesh->asemantic[a] == GI_PARAM_STRETCH_ATTRIB && !bStretch) || 
				(pMesh->asemantic[a] == GI_NONE && pMesh->aoffset[a] < 0))
			{
				GIContext_error(pContext, GI_INVALID_OPERATION);
				return;
			}
			uiAttrib[uiNumAttribs++] = a;
		}
	}

	/* fill arrays */
	do
	{
		for(i=0,pHalfEdge=pFace->hedges; i<3; ++i,++idx,pHalfEdge=pHalfEdge->next)
		{
			for(j=0; j<uiNumAttribs; ++j)
			{
				a = uiAttrib[j];
				dst = pContext->attrib_pointer[a] + 
					pContext->attrib_stride[a]*idx;
				switch(pMesh->asemantic[a])
				{
				case GI_POSITION_ATTRIB:
					GI_VEC3_COPY(dst, pHalfEdge->vstart->coords);
					break;
				case GI_PARAM_ATTRIB:
					GI_VEC2_COPY(dst, pHalfEdge->pstart->params);
					break;
				case GI_PARAM_STRETCH_ATTRIB:
					*dst = pHalfEdge->pstart->stretch;
					break;
				default:
					fsrc = (GIfloat*)((GIbyte*)pHalfEdge->astart+pMesh->aoffset[a]);
					for(c=0; c<pMesh->asize[a]; ++c)
						dst[c] = fsrc[c];
				}
			}
		}
		pFace = pFace->next;
	}while(pFace != pFEnd);
	*vcount = 3 * (pMesh->active_patch ? 
		pMesh->active_patch->fcount : pMesh->fcount);
}

/** \internal
 *  \brief Mesh destructor.
 *  \param mesh mesh to destroy
 *  \ingroup mesh
 */
void GIMesh_destruct(GIMesh *mesh)
{
	GIuint a;

	/* clear lists */
	GIDynamicQueue_destruct(&mesh->split_hedges);
	GI_LIST_CLEAR_PERSISTENT(mesh->faces, sizeof(GIFace));
	GI_LIST_CLEAR_PERSISTENT(mesh->edges, sizeof(GIEdge));
	GI_LIST_CLEAR_PERSISTENT(mesh->vertices, sizeof(GIVertex));
	GI_LIST_CLEAR_PERSISTENT(mesh->attributes, sizeof(GIAttribute)+mesh->attrib_size);

	/* reset properties */
	mesh->fcount = mesh->ecount = mesh->vcount = mesh->acount = 0;
	mesh->attrib_size = 0;
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
	{
		mesh->aoffset[a] = -1;
		mesh->asize[a] = 0;
		mesh->anorm[a] = GI_FALSE;
		mesh->asemantic[a] = GI_NONE;
	}
	memset(mesh->semantic, 0, GI_SEMANTIC_COUNT*sizeof(GIuint));
	mesh->varray_size = mesh->varray_attribs = 0;
	mesh->genus = -1;
	GI_VEC3_SET(mesh->aabb_min, 0.0, 0.0, 0.0);
	GI_VEC3_SET(mesh->aabb_max, 0.0, 0.0, 0.0);
	mesh->radius = 0.0;
	mesh->mean_edge = 0.0;

	/* destroy patches */
	GIMesh_destroy_cut(mesh);
}

/** \internal
 *  \brief Delete mesh cut.
 *  \param mesh mesh to work on
 *  \ingroup mesh
 */
void GIMesh_destroy_cut(GIMesh *mesh)
{
	GIPatch *pPatch = mesh->patches;
	GIFace *pFace;
	GIHalfEdge *pHalfEdge;
	GIuint i;
	if(!mesh->patch_count)
		return;

	/* delete splits made after cut creation */
	GIMesh_revert_splits(mesh, mesh->split_hedges.size-mesh->pre_cut_splits);

	/* destroy patches */
	for(i=0; i<mesh->patch_count; ++i,++pPatch)
		GIPatch_destruct(pPatch);
	GI_FREE_ARRAY(mesh->patches);
	mesh->patches = mesh->active_patch = mesh->varray_patch = NULL;
	mesh->varray_attribs = 0;

	/* reset connections */
	GI_LIST_FOREACH(mesh->faces, pFace)
		GI_LIST_FOREACH(pFace->hedges, pHalfEdge)
			pHalfEdge->pstart = NULL;
			pHalfEdge->vstart->cut_degree = 0;
		GI_LIST_NEXT(pFace->hedges, pHalfEdge)
	GI_LIST_NEXT(mesh->faces, pFace)

	/* reset parameterization */
	memset(mesh->stretch, 0, GI_STRETCH_COUNT*sizeof(GIdouble));
	mesh->param_metric = 0;
	mesh->surface_area = mesh->param_area = 0.0;
	mesh->min_param_stretch = mesh->max_param_stretch = 0.0;
	mesh->patch_count = mesh->param_patches = 0;
	mesh->cut_splits = mesh->pre_cut_splits = 0;
	mesh->resolution = 0;

	/* delete rest of the splits */
	GIMesh_revert_splits(mesh, -1);

	/* reset coordinates if neccessary */
	if(mesh->old_coords)
	{
		GIEdge *pEdge;
		GIVertex *pVertex = mesh->vertices;
		GIdouble *vec = mesh->old_coords;
		for(i=0; i<mesh->vcount; ++i,vec+=3,pVertex=pVertex->next)
		{
			GI_VEC3_COPY(pVertex->coords, vec);
		}
		GI_LIST_FOREACH(mesh->edges, pEdge)
			pEdge->length = GIvec3d_dist(pEdge->hedge[0].vstart->coords, 
				pEdge->hedge[1].vstart->coords);
		GI_LIST_NEXT(mesh->edges, pEdge)
		GI_FREE_ARRAY(mesh->old_coords);
		mesh->old_coords = NULL;
	}
}

/** \internal
 *  \brief Reverse all half edge splits.
 *  \param mesh mesh to work on
 *  \param count number of splits to revert or -1 for all
 *  \ingroup mesh
 */
void GIMesh_revert_splits(GIMesh *mesh, GIint count)
{
	GIPatch *pPatch, *pTwinPatch;
	GIEdge *pEDel;
	GIHalfEdge *pHalfEdge, *pHTwin, *pHDel;
	GIVertex *pVDel;
	GIParam *pPDel1, *pPDel2;
	GIAttribute *pADel1, *pADel2;
	GISplitInfo *pSplit;
	if(count < 0 || count > mesh->split_hedges.size)
		count = mesh->split_hedges.size;
	if(count)
		mesh->varray_attribs = 0;

	/* reverse all splits */
	while(count--)
	{
		/* extract information */
		pSplit = (GISplitInfo*)GIDynamicQueue_pop(&mesh->split_hedges);
		pHalfEdge = pSplit->hedge;
		pPatch = pSplit->patch;
		pTwinPatch = pSplit->twin_patch;
		pHTwin = pHalfEdge->twin;

		/* delete new things */
		pHDel = pHalfEdge->face ? pHalfEdge->next->twin->next : pHalfEdge->next;
		pEDel = pHDel->edge;
		pVDel = pHDel->vstart;
		pPDel1 = pHDel->pstart;
		pADel1 = pHDel->astart;
		if(pPDel1 && pPDel1->cut_hedge == pHDel)
			--pPatch->hcount;
		GIHalfEdge_revert_half_split(pHalfEdge, pHDel, mesh, pPatch, GI_TRUE);

		/* delete new things on other side */
		pHDel = pHTwin->face ? pHTwin->prev->twin->prev : pHTwin->prev;
		assert(pHDel->edge==pEDel);
		pPDel2 = pHTwin->pstart;
		pADel2 = pHTwin->astart;
		pHTwin->vstart = pHDel->vstart;
		pHTwin->astart = pHDel->astart;
		pHTwin->pstart = pHDel->pstart;
		if(pHTwin->vstart->hedge == pHDel)
			pHTwin->vstart->hedge = pHTwin;
		if(pHTwin->pstart && pHTwin->pstart->cut_hedge == pHDel)
			pHTwin->pstart->cut_hedge = pHTwin;
		if(pPDel2 && pPDel2->cut_hedge == pHTwin)
			--pTwinPatch->hcount;
		GIHalfEdge_revert_half_split(pHTwin, pHDel, mesh, pTwinPatch, GI_FALSE);

		/* delete vertex */
		GI_LIST_DELETE_PERSISTENT(mesh->vertices, pVDel, sizeof(GIVertex));
		--mesh->vcount;
		if(pADel1)
		{
			GI_LIST_DELETE_PERSISTENT(mesh->attributes, pADel1, 
				sizeof(GIAttribute)+mesh->attrib_size);
			--mesh->acount;
		}
		if(pPDel1)
		{
			GI_LIST_DELETE_PERSISTENT(pPatch->params, pPDel1, sizeof(GIParam));
			--pPatch->pcount;
		}
		if(pADel2 && pADel1 != pADel2)
		{
			GI_LIST_DELETE_PERSISTENT(mesh->attributes, pADel2, 
				sizeof(GIAttribute)+mesh->attrib_size);
			--mesh->acount;
		}
		if(pPDel2 && pPDel1 != pPDel2)
		{
			GI_LIST_DELETE_PERSISTENT(pTwinPatch->params, pPDel2, sizeof(GIParam));
			--pTwinPatch->pcount;
		}

		/* delete edge and split record */
		GI_LIST_DELETE_PERSISTENT(mesh->edges, pEDel, sizeof(GIEdge));
		--mesh->ecount;
		pHalfEdge->edge->length = GIvec3d_dist(pHalfEdge->vstart->coords, 
			pHalfEdge->next->vstart->coords);
		GI_FREE_PERSISTENT(pSplit, sizeof(GISplitInfo));
	}
}

/** \internal
 *  \brief Compute genus of mesh.
 *  \param mesh mesh to work on
 *  \return genus
 *  \ingroup mesh
 */
GIint GIMesh_genus(GIMesh *mesh)
{
	GIHalfEdge *pHalfEdge;
	GIVertex *pVertex;
	GIboolean *pVisited;
	GIuint B = 0;
	if(mesh->genus >= 0 || !mesh->vertices)
		return mesh->genus;
	pVisited = (GIboolean*)GI_CALLOC_ARRAY(mesh->vcount, sizeof(GIboolean));

	/* count boundary loops */
	GI_LIST_FOREACH(mesh->vertices, pVertex)
		if(!pVisited[pVertex->id])
		{
			/* vertex on boundary? */
			if(pVertex->hedge->face)
				pVisited[pVertex->id] = GI_TRUE;
			else
			{
				/* walk along boundary */
				GI_LIST_FOREACH(pVertex->hedge, pHalfEdge)
					pVisited[pHalfEdge->vstart->id] = GI_TRUE;
				GI_LIST_NEXT(pVertex->hedge, pHalfEdge)
				++B;
			}
		}
	GI_LIST_NEXT(mesh->vertices, pVertex)

	/* Euler-Poincare formula */
	mesh->genus = (2-B-mesh->vcount-mesh->fcount+mesh->ecount) >> 1;
	return mesh->genus;
}

/** \internal
 *  \brief Compute stretch of mesh.
 *  \param mesh mesh to work on
 *  \param metric stretch metric to use
 *  \param param_stretches GI_TRUE to compute param stretches GI_FALSE else
 *  \ingroup mesh
 */
void GIMesh_compute_stretch(GIMesh *mesh, GIuint metric, 
							GIboolean param_stretches)
{
	GIPatch *pPatch;
	GIuint m = metric - GI_STRETCH_BASE;

	/* error checking and initialization */
	if(!mesh->patches)
		return;
	mesh->stretch[m] = mesh->surface_area = mesh->param_area = 0.0;
	if(param_stretches)
	{
		mesh->min_param_stretch = DBL_MAX;
		mesh->max_param_stretch = 0.0;
	}

	/* process patches */
	GI_LIST_FOREACH(mesh->patches, pPatch)
		/* compute patch stretch and accumulate areas */
		if(!pPatch->stretch[m] || 
			(param_stretches && pPatch->param_metric != metric))
			GIPatch_compute_stretch(pPatch, metric, 
				param_stretches, GI_TRUE);
		mesh->surface_area += pPatch->surface_area;
		mesh->param_area += pPatch->param_area;

		/* compute overall stretch */
		switch(metric)
		{
			case GI_RMS_GEOMETRIC_STRETCH:
			case GI_COMBINED_STRETCH:
				mesh->stretch[m] += pPatch->stretch[m] * 
					pPatch->stretch[m] * pPatch->surface_area;
				break;
			case GI_MAX_GEOMETRIC_STRETCH:
				if(pPatch->stretch[m] > mesh->stretch[m])
					mesh->stretch[m] = pPatch->stretch[m];
		}

		/* record extrema */
		if(param_stretches)
		{
			if(pPatch->min_param_stretch < mesh->min_param_stretch)
				mesh->min_param_stretch = pPatch->min_param_stretch;
			if(pPatch->max_param_stretch > mesh->max_param_stretch)
				mesh->max_param_stretch = pPatch->max_param_stretch;
		}
	GI_LIST_NEXT(mesh->patches, pPatch)

	/* scale down */
	if(metric != GI_MAX_GEOMETRIC_STRETCH)
		mesh->stretch[m] = sqrt(mesh->stretch[m]/mesh->surface_area);
	if(param_stretches)
		mesh->param_metric = metric;
}

/** \internal
 *  \brief Find half edge by index in face.
 *  \param face face to look at
 *  \param i index of half edge in face
 *  \return half edge at specified position
 *  \ingroup mesh
 */
GIHalfEdge* GIFace_halfedge_at(GIFace *face, GIuint i)
{
	/* search index */
	GIHalfEdge *pHalfEdge = face->hedges;
	for(; i; --i)
		pHalfEdge = pHalfEdge->next;
	return pHalfEdge;
}

/** \internal
 *  \brief Compute center of face.
 *  \param face face to compute center of
 *  \param center vector to store center coordinates
 *  \ingroup mesh
 */
void GIFace_center(GIFace *face, GIdouble *center)
{
	/* compute average of vertices */
	GI_VEC3_ADD(center, face->hedges->vstart->coords, 
		face->hedges->next->vstart->coords);
	GI_VEC3_ADD(center, center, face->hedges->prev->vstart->coords);
	GI_VEC3_SCALE(center, center, 0.33333333333333333333);
}

/** \internal
 *  \brief Compute area of triangular face.
 *  \param face face to compute area of
 *  \return area of face;
 *  \ingroup mesh
 */
GIdouble GIFace_area(GIFace *face)
{
	GIHalfEdge *pHalfEdge = face->hedges;
	GIdouble v12[3], v13[3], n[3];
	GIdouble *p1 = pHalfEdge->vstart->coords;
	GIdouble *p2 = pHalfEdge->next->vstart->coords;
	GIdouble *p3 = pHalfEdge->prev->vstart->coords;

	/* compute area */
	GI_VEC3_SUB(v12, p2, p1);
	GI_VEC3_SUB(v13, p3, p1);
	GI_VEC3_CROSS(n, v12, v13);
	return 0.5 * GI_VEC3_LENGTH(n);
}

/** \internal
 *  \brief Compute stretch of triangular face.
 *  \param face face to work on
 *  \param metric stretch metric to use
 *  \param args additional arguments
 *  \param area_2d address to store parameter space area at or NULL if not needed
 *  \return stretch of face
 *  \ingroup mesh
 */
GIdouble GIFace_stretch(GIFace *face, GIenum metric, GIvoid *args, GIdouble *area_2d)
{
	GIHalfEdge *pHalfEdge = face->hedges;
	GIdouble *p1, *p2, *p3, *q1, *q2, *q3;
#if OPENGI_SSE >= 2
	__m128d XMM0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6, XMM7;
	GIdouble dStretch, dWeight = args ? *(GIdouble*)args : 0.0;
#else
	GIdouble Fu[3], Fv[3];
	GIdouble u32, u13, u21, v23, v31, v12, dOne2A;
	GIdouble E, F, G, EG;
	GIdouble dA2D;
#endif

	/* get coordinates */
	p1 = pHalfEdge->pstart->params;
	q1 = pHalfEdge->vstart->coords;
	pHalfEdge = pHalfEdge->next;
	p2 = pHalfEdge->pstart->params;
	q2 = pHalfEdge->vstart->coords;
	pHalfEdge = pHalfEdge->next;
	p3 = pHalfEdge->pstart->params;
	q3 = pHalfEdge->vstart->coords;

#if OPENGI_SSE >= 2
	/* compute auxiliary values and area */
	XMM0 = _mm_load_sd(&p2[0]);
	XMM0 = _mm_loadh_pd(XMM0, &p2[1]);
	XMM1 = _mm_load_sd(&p3[0]);
	XMM1 = _mm_loadh_pd(XMM1, &p3[1]);
	XMM2 = _mm_load_sd(&p1[0]);
	XMM2 = _mm_loadh_pd(XMM2, &p1[1]);
	XMM3 = XMM0;
	XMM0 = _mm_sub_pd(XMM0, XMM1);
	XMM1 = _mm_sub_pd(XMM1, XMM2);
	XMM2 = _mm_sub_pd(XMM2, XMM3);
	XMM3 = _mm_set_sd(-0.0);
	XMM0 = _mm_xor_pd(XMM0, XMM3);
	XMM1 = _mm_xor_pd(XMM1, XMM3);
	XMM2 = _mm_xor_pd(XMM2, XMM3);
	XMM4 = _mm_shuffle_pd(XMM1, XMM1, _MM_SHUFFLE2(0, 1));
	XMM4 = _mm_mul_pd(XMM4, XMM2);
#if OPENGI_SSE >= 3
	XMM4 = _mm_hsub_pd(XMM4, XMM4);
#else
	XMM3 = _mm_shuffle_pd(XMM4, XMM4, _MM_SHUFFLE2(1, 1));
	XMM4 = _mm_shuffle_pd(XMM4, XMM4, _MM_SHUFFLE2(0, 0));
	XMM4 = _mm_sub_pd(XMM4, XMM3);
#endif
	if(area_2d)
	{
		XMM3 = _mm_set_sd(0.5);
		XMM3 = _mm_mul_sd(XMM3, XMM4);
		_mm_store_sd(area_2d, XMM3);
	}
	XMM3 = _mm_set1_pd(1.0);
	XMM3 = _mm_div_pd(XMM3, XMM4);

	/* compute partial derivatives (XMM0 = 32, XMM1 = 13, XMM2 = 21, XMM3 = 1/A) */
	XMM5 = _mm_set1_pd(q1[0]);
	XMM5 = _mm_mul_pd(XMM5, XMM0);
	XMM4 = _mm_set1_pd(q2[0]);
	XMM4 = _mm_mul_pd(XMM4, XMM1);
	XMM5 = _mm_add_pd(XMM5, XMM4);
	XMM4 = _mm_set1_pd(q3[0]);
	XMM4 = _mm_mul_pd(XMM4, XMM2);
	XMM5 = _mm_add_pd(XMM5, XMM4);
	XMM5 = _mm_mul_pd(XMM5, XMM3);
	XMM6 = _mm_set1_pd(q1[1]);
	XMM6 = _mm_mul_pd(XMM6, XMM0);
	XMM4 = _mm_set1_pd(q2[1]);
	XMM4 = _mm_mul_pd(XMM4, XMM1);
	XMM6 = _mm_add_pd(XMM6, XMM4);
	XMM4 = _mm_set1_pd(q3[1]);
	XMM4 = _mm_mul_pd(XMM4, XMM2);
	XMM6 = _mm_add_pd(XMM6, XMM4);
	XMM6 = _mm_mul_pd(XMM6, XMM3);
	XMM7 = _mm_set1_pd(q1[2]);
	XMM7 = _mm_mul_pd(XMM7, XMM0);
	XMM4 = _mm_set1_pd(q2[2]);
	XMM4 = _mm_mul_pd(XMM4, XMM1);
	XMM7 = _mm_add_pd(XMM7, XMM4);
	XMM4 = _mm_set1_pd(q3[2]);
	XMM4 = _mm_mul_pd(XMM4, XMM2);
	XMM7 = _mm_add_pd(XMM7, XMM4);
	XMM7 = _mm_mul_pd(XMM7, XMM3);

	/* compute tangents (XMM5 = Fx[0], XMM6 = Fx[1], XMM7 = Fx[2]) */
	XMM0 = _mm_mul_pd(XMM5, XMM5);
	XMM1 = _mm_mul_pd(XMM6, XMM6);
	XMM2 = _mm_mul_pd(XMM7, XMM7);
	XMM0 = _mm_add_pd(XMM0, XMM1);
	XMM0 = _mm_add_pd(XMM0, XMM2);
	if(metric != GI_RMS_GEOMETRIC_STRETCH)
	{
		XMM1 = _mm_shuffle_pd(XMM5, XMM5, _MM_SHUFFLE2(0, 1));
		XMM2 = _mm_shuffle_pd(XMM6, XMM6, _MM_SHUFFLE2(0, 1));
		XMM3 = _mm_shuffle_pd(XMM7, XMM7, _MM_SHUFFLE2(0, 1));
		XMM1 = _mm_mul_pd(XMM1, XMM5);
		XMM2 = _mm_mul_pd(XMM2, XMM6);
		XMM3 = _mm_mul_pd(XMM3, XMM7);
		XMM1 = _mm_add_pd(XMM1, XMM2);
		XMM1 = _mm_add_pd(XMM1, XMM3);
	}

	/* compute stretch values (XMM0 = (E,G), XMM1 = (F,F)) */
	switch(metric)
	{
	/* L2-stretch */
	case GI_RMS_GEOMETRIC_STRETCH:
#if OPENGI_SSE >= 3
		XMM0 = _mm_hadd_pd(XMM0, XMM0);
#else
		XMM7 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(1, 1));
		XMM0 = _mm_add_sd(XMM0, XMM7);
#endif
		XMM7 = _mm_set_sd(0.5);
		XMM0 = _mm_mul_sd(XMM0, XMM7);
		_mm_store_sd(&dStretch, XMM0);
		return dStretch;

	/* combined energy */
	case GI_COMBINED_STRETCH:
#if OPENGI_SSE >= 3
		XMM2 = _mm_hsub_pd(XMM0, XMM0);
		XMM0 = _mm_hadd_pd(XMM0, XMM0);
#else
		XMM3 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(0, 1));
		XMM2 = _mm_sub_pd(XMM0, XMM3);
		XMM0 = _mm_add_pd(XMM0, XMM3);
#endif
		XMM3 = _mm_set1_pd(0.5);
		XMM4 = _mm_set1_pd(4.0);
		XMM2 = _mm_mul_pd(XMM2, XMM2);
		XMM1 = _mm_mul_pd(XMM1, XMM1);
		XMM1 = _mm_mul_pd(XMM1, XMM4);
		XMM1 = _mm_add_pd(XMM1, XMM2);
		XMM1 = _mm_sqrt_pd(XMM1);
		XMM2 = _mm_set_sd(-0.0);
		XMM1 = _mm_xor_pd(XMM1, XMM2);
		XMM0 = _mm_add_pd(XMM0, XMM1);
		XMM0 = _mm_mul_pd(XMM0, XMM3);
		XMM1 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(1, 1));
		XMM2 = _mm_div_sd(XMM1, XMM0);
		XMM0 = _mm_mul_sd(XMM0, XMM1);
		XMM0 = _mm_shuffle_pd(XMM0, XMM2, _MM_SHUFFLE2(0, 0));
		XMM0 = _mm_sqrt_pd(XMM0);
		XMM1 = _mm_set1_pd(1.0);
		XMM1 = _mm_div_pd(XMM1, XMM0);
		XMM0 = _mm_add_pd(XMM0, XMM1);
		XMM1 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(1, 1));
		if(fabs(dWeight-1.0) > 1e-4)
		{
			_mm_store_sd(&dStretch, XMM0);
			dStretch = pow(dStretch, dWeight);
			XMM0 = _mm_set_sd(dStretch);
		}
		XMM0 = _mm_mul_sd(XMM0, XMM1);
		_mm_store_sd(&dStretch, XMM0);
		return dStretch;

	/* Linf-stretch */
	case GI_MAX_GEOMETRIC_STRETCH:
#if OPENGI_SSE >= 3
		XMM2 = _mm_hsub_pd(XMM0, XMM0);
		XMM0 = _mm_hadd_pd(XMM0, XMM0);
#else
		XMM3 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(0, 1));
		XMM2 = _mm_sub_pd(XMM0, XMM3);
		XMM0 = _mm_add_pd(XMM0, XMM3);
#endif
		XMM3 = _mm_set_sd(0.5);
		XMM4 = _mm_set_sd(4.0);
		XMM2 = _mm_mul_sd(XMM2, XMM2);
		XMM1 = _mm_mul_sd(XMM1, XMM1);
		XMM1 = _mm_mul_sd(XMM1, XMM4);
		XMM1 = _mm_add_sd(XMM1, XMM2);
		XMM1 = _mm_sqrt_sd(XMM1, XMM1);
		XMM0 = _mm_add_sd(XMM0, XMM1);
		XMM0 = _mm_mul_sd(XMM0, XMM3);
		XMM0 = _mm_sqrt_sd(XMM0, XMM0);
		_mm_store_sd(&dStretch, XMM0);
		return dStretch;
#else
	/* compute auxiliary values and area */
	u32 = p3[0] - p2[0];
	u13 = p1[0] - p3[0];
	u21 = p2[0] - p1[0];
	v23 = p2[1] - p3[1];
	v31 = p3[1] - p1[1];
	v12 = p1[1] - p2[1];
	dA2D = (u21*v31 - u13*v12);
	dOne2A = 1.0 / dA2D;
	if(area_2d)
		*area_2d = 0.5 * dA2D;
//	if(dA2D < 0.0)
//		return DBL_MAX;

	/* compute partial derivatives and metric tensor */
	Fu[0] = (v23*q1[0]+v31*q2[0]+v12*q3[0]) * dOne2A;
	Fu[1] = (v23*q1[1]+v31*q2[1]+v12*q3[1]) * dOne2A;
	Fu[2] = (v23*q1[2]+v31*q2[2]+v12*q3[2]) * dOne2A;
	Fv[0] = (u32*q1[0]+u13*q2[0]+u21*q3[0]) * dOne2A;
	Fv[1] = (u32*q1[1]+u13*q2[1]+u21*q3[1]) * dOne2A;
	Fv[2] = (u32*q1[2]+u13*q2[2]+u21*q3[2]) * dOne2A;
	E = GI_VEC3_DOT(Fu, Fu);
	G = GI_VEC3_DOT(Fv, Fv);

	/* compute stretch values */
	switch(metric)
	{
	/* L2-stretch */
	case GI_RMS_GEOMETRIC_STRETCH:
		return 0.5 * (E+G);

	/* combined energy */
	case GI_COMBINED_STRETCH:
		{
			GIdouble dSqrt, dRatio, dDet, d1, d2, dWeight = *(GIdouble*)args;
			F = GI_VEC3_DOT(Fu, Fv);
			EG = E - G;
			dSqrt = sqrt(EG*EG+4.0*F*F);
			d1 = 0.5 * (E+G+dSqrt);
			d2 = 0.5 * (E+G-dSqrt);
			dRatio = sqrt(d1/d2);
			dDet = sqrt(d1*d2);
			if(fabs(dWeight-1.0) > 1e-4)
				return (dRatio+1.0/dRatio) * pow(dDet+1.0/dDet, dWeight);
			return (dRatio+1.0/dRatio) * (dDet+1.0/dDet);
		}

	/* Linf-stretch */
	case GI_MAX_GEOMETRIC_STRETCH:
		F = GI_VEC3_DOT(Fu, Fv);
		EG = E - G;
		return sqrt(0.5*(E+G+sqrt(EG*EG+4.0*F*F)));
#endif

	/* invalid enum */
	default:
		return 0.0;
	}
}

/** \internal
 *  \brief Find index of half edge in half edge's face.
 *  \param hedge half edge to look for
 *  \return index of half edge in face
 *  \ingroup mesh
 */
GIuint GIHalfEdge_index(GIHalfEdge *hedge)
{
	/* search half edge */
	GIHalfEdge *pHalfEdge;
	GIuint i = 0;
	for(pHalfEdge=hedge->face->hedges; 
		pHalfEdge!=hedge; pHalfEdge=pHalfEdge->next)
		++i;
	return i;
}

/** \internal
 *  \brief Split half edge and its opposite and incident faces.
 *  \param hedge half edge to split
 *  \param patch patch half edge belongs to
 *  \param twin_patch patch twin half edge belongs to
 *  \param f interpolation factor
 *  \param params params of new vertex or NULL to interpolate
 *  \ingroup mesh
 */
void GIHalfEdge_split(GIHalfEdge *hedge, struct _GIPatch *patch, 
					  struct _GIPatch *twin_patch, 
					  GIdouble f, GIdouble *params)
{
	GIMesh *pMesh = patch ? patch->mesh : (GIMesh*)twin_patch;
	GIEdge *pEdge1 = hedge->edge, *pEdge2;
	GIHalfEdge *pHNew1, *pHNew2, *pHTwin = hedge->twin;
	GIVertex *pVertex;
	GIParam *pParam = NULL;
	GIAttribute *pAttribute = NULL;
	GISplitInfo *pSplit;
	GIdouble vec[3];
	GIdouble dOneF = 1.0 - f;
	GIboolean bCut = (hedge->pstart && hedge->pstart->cut_hedge == hedge) || 
		(pHTwin->pstart && pHTwin->pstart->cut_hedge == pHTwin);
	if(!patch)
		twin_patch = NULL;

	/* create and connect new edge */
	pEdge2 = (GIEdge*)GI_MALLOC_PERSISTENT(sizeof(GIEdge));
	GI_LIST_ADD(pMesh->edges, pEdge2);
	pEdge2->id = pMesh->ecount++;
	pEdge2->length = dOneF * pEdge1->length;
	pEdge1->length *= f;
	pHNew1 = &pEdge2->hedge[hedge->face ? 0 : 1];
	pHNew2 = &pEdge2->hedge[hedge->face ? 1 : 0];
	pHNew1->edge = pEdge2;
	pHNew2->edge = pEdge2;
	pHNew1->twin = pHNew2;
	pHNew2->twin = pHNew1;

	/* create center vertex */
	pVertex = (GIVertex*)GI_MALLOC_PERSISTENT(sizeof(GIVertex));
	GI_LIST_ADD(pMesh->vertices, pVertex);
	pVertex->id = pMesh->vcount++;
	GI_VEC3_SCALE(pVertex->coords, hedge->vstart->coords, dOneF);
	GI_VEC3_ADD_SCALED(pVertex->coords, 
		pVertex->coords, hedge->next->vstart->coords, f);
	if(pMesh->anorm[pMesh->semantic[GI_POSITION_ATTRIB-GI_SEMANTIC_BASE]])
		GIvec3d_normalize(pVertex->coords);
	pVertex->flags = 0;
	pVertex->cut_degree = bCut ? 2 : 0;

	/* create center attribute */
	if(hedge->astart)
	{
		pAttribute = GIAttribute_create_interpolated(
			hedge->astart, hedge->next->astart, dOneF, pMesh);
		GI_LIST_ADD(pMesh->attributes, pAttribute);
		pAttribute->id = pMesh->acount++;
	}

	/* create center param */
	if(hedge->pstart)
	{
		pParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
		GI_LIST_ADD(patch->params, pParam);
		pParam->id = patch->pcount++;
		if(params)
		{
			GI_VEC2_COPY(pParam->params, params);
		}
		else
		{
			GI_VEC2_SCALE(pParam->params, hedge->pstart->params, dOneF);
			GI_VEC2_ADD_SCALED(pParam->params, pParam->params, 
				hedge->next->pstart->params, f);
		}
		pParam->stretch = 0.0;
		pParam->vertex = pVertex;
		if(bCut)
		{
			pParam->cut_hedge = pHNew1;
			++patch->hcount;
		}
		else
			pParam->cut_hedge = NULL;
	}

	/* split one side of half edge */
	GIHalfEdge_half_split(hedge, pHNew1, pMesh, patch, GI_TRUE);
	pHNew1->vstart = pVertex;
	pHNew1->astart = pAttribute;
	pHNew1->pstart = pParam;
	if(hedge->face)
	{
		GIHalfEdge *pHMid = hedge->next;
		pHMid->vstart = pVertex;
		pHMid->astart = pAttribute;
		pHMid->pstart = pParam;
		GI_VEC3_SUB(vec, pHMid->twin->vstart->coords, pVertex->coords);
		pHMid->edge->length = GI_VEC3_LENGTH(vec);
	}

	/* other center attribute? */
	if(!pHTwin->astart)
		pAttribute = NULL;
	else if(hedge->astart != pHTwin->next->astart || 
		pHNew1->next->astart != pHTwin->astart)
	{
		pAttribute = GIAttribute_create_interpolated(
			pHTwin->astart, pHTwin->next->astart, f, pMesh);
		GI_LIST_ADD(pMesh->attributes, pAttribute);
		pAttribute->id = pMesh->acount++;
	}

	/* other center param? */
	if(!pHTwin->pstart)
		pParam = NULL;
	else if(bCut)
	{
		pParam = (GIParam*)GI_MALLOC_PERSISTENT(sizeof(GIParam));
		GI_LIST_ADD(twin_patch->params, pParam);
		pParam->id = twin_patch->pcount++;
		GI_VEC2_SCALE(pParam->params, pHTwin->pstart->params, f);
		GI_VEC2_ADD_SCALED(pParam->params, pParam->params, 
			pHTwin->next->pstart->params, dOneF);
		pParam->stretch = 0.0;
		pParam->vertex = pVertex;
		pParam->cut_hedge = pHTwin;
		++twin_patch->hcount;
	}

	/* split opposite side of half edge */
	GIHalfEdge_half_split(pHTwin, pHNew2, pMesh, twin_patch, GI_FALSE);
	pHNew2->vstart = pHTwin->vstart;
	pHNew2->astart = pHTwin->astart;
	pHNew2->pstart = pHTwin->pstart;
	pHTwin->vstart = pVertex;
	pHTwin->astart = pAttribute;
	pHTwin->pstart = pParam;
	if(pHTwin->face)
	{
		GIHalfEdge *pHMid = pHNew2->next;
		pHMid->vstart = pVertex;
		pHMid->astart = pAttribute;
		pHMid->pstart = pParam;
		GI_VEC3_SUB(vec, pHMid->twin->vstart->coords, pVertex->coords);
		pHMid->edge->length = GI_VEC3_LENGTH(vec);
		pVertex->hedge = pHNew1;
	}
	else
		pVertex->hedge = pHTwin;
	if(pHNew2->vstart->hedge == pHTwin)
		pHNew2->vstart->hedge = pHNew2;
	if(pHNew2->pstart && pHNew2->pstart->cut_hedge == pHTwin)
		pHNew2->pstart->cut_hedge = pHNew2;

	/* save split */
	pSplit = (GISplitInfo*)GI_MALLOC_PERSISTENT(sizeof(GISplitInfo));
	pSplit->hedge = hedge;
	pSplit->patch = patch;
	pSplit->twin_patch = patch ? twin_patch : NULL;
	GIDynamicQueue_push(&pMesh->split_hedges, pSplit);
	pMesh->varray_attribs = 0;
}

/** \internal
 *  \brief Split one side of half edge.
 *  \param hedge half edge to split
 *  \param hnew newly created half edge
 *  \param mesh mesh to work on
 *  \param patch patch half edge belongs to, if any
 *  \param hedge_first \c GI_TRUE if \a hedge before \a hnew, \c GI_FALSE else
 *  \ingroup mesh
 */
void GIHalfEdge_half_split(GIHalfEdge *hedge, GIHalfEdge *hnew, 
						   GIMesh *mesh, struct _GIPatch *patch, 
						   GIboolean hedge_first)
{
	GIFace *pFace = hedge->face, *pFNew;

	/* boundary halfedge? */
	if(pFace)
	{
		/* create edge and connect half edges */
		GIHalfEdge *pHMov = hedge_first ? hedge->next : hedge->prev, 
			*pHIns = pHMov->next, *pHNew1, *pHNew2, *pHTmp;
		GIEdge *pENew = (GIEdge*)GI_MALLOC_PERSISTENT(sizeof(GIEdge));
		GI_LIST_ADD(mesh->edges, pENew);
		pENew->id = mesh->ecount++;
		pHNew1 = &pENew->hedge[0];
		pHNew2 = &pENew->hedge[1];
		pHNew1->edge = pENew;
		pHNew2->edge = pENew;
		pHNew1->twin = pHNew2;
		pHNew2->twin = pHNew1;
		pHNew2->vstart = hedge->prev->vstart;
		pHNew2->astart = hedge->prev->astart;
		pHNew2->pstart = hedge->prev->pstart;
		if(!hedge_first)
		{
			GI_SWAP(pHNew1, pHNew2, pHTmp);
		}

		/* make existing face first face */
		GI_LIST_REMOVE(pFace->hedges, pHMov);
		GI_LIST_INSERT(pFace->hedges, pHIns, pHNew1);
		pHNew1->face = pFace;

		/* create second face */
		pFNew = (GIFace*)GI_MALLOC_PERSISTENT(sizeof(GIFace));
		if(patch)
		{
			GI_LIST_INSERT(mesh->faces, patch->next->faces, pFNew);
			++patch->fcount;
		}
		else
		{
			GI_LIST_ADD(mesh->faces, pFNew);
		}
		pFNew->id = mesh->fcount++;
		pFNew->hedges = NULL;
		GI_LIST_ADD(pFNew->hedges, pHMov);
		if(hedge_first)
		{
			GI_LIST_ADD(pFNew->hedges, pHNew2);
			GI_LIST_ADD(pFNew->hedges, hnew);
		}
		else
		{
			GI_LIST_ADD(pFNew->hedges, hnew);
			GI_LIST_ADD(pFNew->hedges, pHNew2);
		}
		pHMov->face = pFNew;
		pHNew2->face = pFNew;
		hnew->face = pFNew;
	}
	else
	{
		GIHalfEdge *pHNext = hedge_first ? hedge->next : hedge;
		GI_LIST_ADD(pHNext, hnew);
		hnew->face = NULL;
	}
}

/** \internal
 *  \brief Revert one side of half edge split.
 *  \param hedge half edge to split
 *  \param hdel half edge to delete
 *  \param mesh mesh to work on
 *  \param patch patch half edge belongs to, if any
 *  \param hedge_first \c GI_TRUE if \a hedge before \a hdel, \c GI_FALSE else
 *  \ingroup mesh
 */
void GIHalfEdge_revert_half_split(GIHalfEdge *hedge, GIHalfEdge *hdel, 
								  GIMesh *mesh, struct _GIPatch *patch, 
								  GIboolean hedge_first)
{
	GIFace *pFace = hedge->face;

	/* boundary halfedge? */
	if(pFace)
	{
		/* delete new things */
		GIHalfEdge *pHMov = hedge_first ? hdel->next : hdel->prev, 
			*pHDel = hedge_first ? hedge->next : hedge->prev, *pHIns = pHDel->next;
		GIFace *pFDel = pHDel->twin->face;
		GIEdge *pEDel = pHDel->edge;
		GI_LIST_REMOVE(pFDel->hedges, pHMov);
		GI_LIST_REMOVE(pFace->hedges, pHDel);
		GI_LIST_INSERT(pFace->hedges, pHIns, pHMov);
		pHMov->face = pFace;
		GI_LIST_DELETE_PERSISTENT(mesh->faces, pFDel, sizeof(GIFace));
		GI_LIST_DELETE_PERSISTENT(mesh->edges, pEDel, sizeof(GIEdge));
		--mesh->fcount;
		--mesh->ecount;
		if(patch)
			--patch->fcount;
	}
	else
	{
		GI_LIST_REMOVE(hdel->prev, hdel);
	}
}

/** \internal
 *  \brief Create new attribute as interpolation of two attributes.
 *  \param attrib1 first attribute
 *  \param attrib2 second attribute
 *  \param f weight of first attribute in [0,1]
 *  \param mesh mesh the attributes belong to
 *  \return pointer to newly created attribute
 *  \ingroup mesh
 */
GIAttribute* GIAttribute_create_interpolated(GIAttribute *attrib1, 
											 GIAttribute *attrib2, 
											 float f, GIMesh *mesh)
{
	GIAttribute *pAttribute;
	GIfloat *vec, *vec1 = (GIfloat*)((GIbyte*)attrib1+sizeof(GIAttribute)), 
		*vec2 = (GIfloat*)((GIbyte*)attrib2+sizeof(GIAttribute));
	GIfloat fOneF = 1.0f - f;
	GIuint i, iFloatCount = mesh->attrib_size / sizeof(GIfloat);

	/* create attribute and interpolate data (renormalize normals) */
	pAttribute = (GIAttribute*)GI_MALLOC_PERSISTENT(sizeof(GIAttribute)+mesh->attrib_size);
	vec = (GIfloat*)((GIbyte*)pAttribute+sizeof(GIAttribute));
	for(i=0; i<iFloatCount; ++i)
		vec[i] = f*vec1[i] + fOneF*vec2[i];
	for(i=0; i<GI_ATTRIB_COUNT; ++i)
	{
		if(mesh->anorm[i])
		{
			if(mesh->asize[i] >= 3)
				GIvec3f_normalize((GIfloat*)((GIbyte*)pAttribute+mesh->aoffset[i]));
			else if(mesh->asize[i] == 2)
				GIvec2f_normalize((GIfloat*)((GIbyte*)pAttribute+mesh->aoffset[i]));
			else
				*(GIfloat*)((GIbyte*)pAttribute+mesh->aoffset[i]) = GI_SIGN(*(GIfloat*)
					((GIbyte*)pAttribute+mesh->aoffset[i])<0.0f);
		}
	}
	return pAttribute;
}

/** \internal
 *  \brief Compute stretch of parameter coordinate.
 *  \param param param to work on
 *  \param metric stretch metric to use
 *  \param args additional arguments
 *  \param face_areas hash of face areas or NULL
 *  \return stretch of param
 *  \ingroup mesh
 */
GIdouble GIParam_stretch(GIParam *param, GIenum metric, 
						 GIvoid *args, GIHash *face_areas)
{
	GIHalfEdge *pHalfEdge = param->cut_hedge, *pHEnd = NULL;
	GIdouble dStretch, dA3D, dResult = 0.0;
	if(!pHalfEdge)
		pHalfEdge = pHEnd = param->vertex->hedge;

	/* average face stretches */
	do
	{
		dStretch = GIFace_stretch(pHalfEdge->face, metric, args, NULL);
		switch(metric)
		{
		/* L2-stretch or combined energy */
		case GI_RMS_GEOMETRIC_STRETCH:
		case GI_COMBINED_STRETCH:
			if(face_areas)
				dA3D = *(GIdouble*)GIHash_find(face_areas, pHalfEdge->face);
			else
				dA3D = GIFace_area(pHalfEdge->face);
			dResult += dStretch * dA3D;
			break;

		/* Linf-stretch */
		case GI_MAX_GEOMETRIC_STRETCH:
			if(dStretch > dResult)
				dResult = dStretch;
		}
		pHalfEdge = pHalfEdge->prev->twin;
	}while(pHalfEdge != pHEnd && pHalfEdge->pstart == param);
	return dResult;
}
