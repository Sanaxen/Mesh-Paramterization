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
 *  \brief Implementation of structures and functions for multiresolution modeling.
 */

#include "gi_multiresolution.h"
#include "gi_math.h"
#include "gi_container.h"
#include "gi_numerics.h"


GIuint GIMultiresolution_face_clustering(GIMesh *mesh, 
										 GIFaceCluster *clusters, 
										 GIdouble orientation, 
										 GIdouble shape, 
										 GIdouble max_error, 
										 GIuint min_clusters, 
										 GIuint max_clusters)
{
	GIFaceCluster *pCluster, *pCluster2;
	GIClusterBoundary *pBoundary, *pBoundary2;
	GIClusterMerge *pMerge;
	GIFace *pFace;
	GIEdge *pEdge;
	GIHalfEdge *pHedge1, *pHedge2, *pHedge3;
	GIFaceCluster **pFaceClusterMap;
	GIHeap qFringe;
	GIHash hNeighbours;
	GIdouble dError;
	GIuint i = 0, uiClusters = mesh->fcount;
	GIdouble v12[3], v13[3], q[10];

	printf("create initial clusters start\n");
	/* create initial clusters */
	GI_LIST_CLEAR(clusters, sizeof(GIFaceCluster));
	pFaceClusterMap = (GIFaceCluster**)GI_MALLOC_ARRAY(
		mesh->fcount, sizeof(GIFaceCluster*));
	GI_LIST_FOREACH(mesh->faces, pFace)
		/* create cluster for face */
		pFaceClusterMap[pFace->id] = pCluster = 
			(GIFaceCluster*)GI_MALLOC_SINGLE(sizeof(GIFaceCluster));
		GI_LIST_ADD(clusters, pCluster);
		pCluster->fcount = 1;
		pCluster->faces = pFace;
		pCluster->boundaries = NULL;

		/* compute area and perimeter */
		pHedge1 = pFace->hedges;
		pHedge2 = pFace->hedges->next;
		pHedge3 = pFace->hedges->prev;
		GI_VEC3_SUB(v12, pHedge2->vstart->coords, pHedge1->vstart->coords);
		GI_VEC3_SUB(v13, pHedge3->vstart->coords, pHedge1->vstart->coords);
		GI_VEC3_CROSS(pCluster->normal, v12, v13);
		pCluster->area = 0.5 * GI_VEC3_LENGTH(pCluster->normal);
		pCluster->perimeter = pHedge1->edge->length + 
			pHedge2->edge->length + pHedge3->edge->length;
		GIvec3d_normalize(pCluster->normal);

		/* compute quadrics */
		GI_QUADRIC_CONSTRUCT(pCluster->P, pHedge1->vstart->coords, 1.0);
		GI_QUADRIC_CONSTRUCT(q, pHedge2->vstart->coords, 1.0);
		GI_QUADRIC_ADD(pCluster->P, pCluster->P, q);
		GI_QUADRIC_CONSTRUCT(q, pHedge3->vstart->coords, 1.0);
		GI_QUADRIC_ADD(pCluster->P, pCluster->P, q);
		GI_QUADRIC_CONSTRUCT(pCluster->R, pCluster->normal, -1.0);
		GI_QUADRIC_SCALE(pCluster->R, pCluster->R, pCluster->area);
	GI_LIST_NEXT(mesh->faces, pFace)
	printf("create initial clusters end\n");

	printf("create cluster boundaries start\n");
	/* create intial merge operations and cluster boundaries */
	GIHeap_construct(&qFringe, mesh->ecount, lessd, -DBL_MAX, GI_TRUE);
	GI_LIST_FOREACH(mesh->edges, pEdge)
		pMerge = (GIClusterMerge*)GI_MALLOC_SINGLE(sizeof(GIClusterMerge));
		pMerge->length = pEdge->length;
		for(i=0; i<2; ++i)
		{
			pMerge->boundary[i].twin = &pMerge->boundary[1-i];
			pMerge->boundary[i].merge = pMerge;
			if(pEdge->hedge[i].face)
			{
				pMerge->boundary[i].cluster = 
					pFaceClusterMap[pEdge->hedge[i].face->id];
				GI_LIST_ADD(pMerge->boundary[i].cluster
					->boundaries, &pMerge->boundary[i]);
			}
			else
				pMerge->boundary[i].cluster = NULL;
		}
		GIHeap_enqueue(&qFringe, pMerge, 
			GIClusterMerge_error(pMerge, orientation, shape));
	GI_LIST_NEXT(mesh->edges, pEdge)
	printf("create cluster boundaries end\n");

	printf("cluster merging start\n");
	/* greedy cluster merging */
	GIHash_construct(&hNeighbours, 64, 0.0f, sizeof(GIFaceCluster*), 
		hash_pointer, compare_pointer, copy_pointer);
	while(qFringe.count && uiClusters > min_clusters)
	{
		printf("  qFringe.count:%d uiClusters:%d %d\n", qFringe.count, uiClusters, min_clusters);

		pMerge = (GIClusterMerge*)GIHeap_dequeue(&qFringe, &dError);
		if(dError > max_error && uiClusters <= max_clusters)
			break;

		/* merge cluster data */
		pCluster = pMerge->boundary[0].cluster;
		pCluster2 = pMerge->boundary[1].cluster;
		GI_QUADRIC_ADD(pCluster->P, pCluster->P, pCluster2->P);
		GI_QUADRIC_ADD(pCluster->R, pCluster->P, pCluster2->R);
		pCluster->area += pCluster2->area;
		pCluster->perimeter += pCluster2->perimeter - 2.0*pMerge->length;
		GI_VEC3_SCALE(pCluster->normal, pCluster->normal, 
			(GIdouble)(pCluster->fcount));
		GI_VEC3_ADD_SCALED(pCluster->normal, pCluster->normal, 
			pCluster2->normal, (GIdouble)(pCluster2->fcount));
		GIvec3d_normalize(pCluster->normal);

		/* merge face sublists */
		pCluster->fcount += pCluster2->fcount;
		if(pCluster2->next == pCluster)
			pCluster->faces = pCluster2->faces;
		else if(pCluster->next != pCluster2)
		{
			pFace = pCluster2->next->faces->prev;
			pCluster2->next->faces->prev = pCluster2->faces->prev;
			pCluster2->faces->prev->next = pCluster2->next->faces;
			pCluster->next->faces->prev->next = pCluster2->faces;
			pCluster2->faces->prev = pCluster->next->faces->prev;
			pCluster->next->faces->prev = pFace;
			pFace->next = pCluster->next->faces;
		}
		GI_LIST_DELETE(clusters, pCluster2, sizeof(GIFaceCluster));
		--uiClusters;

		/* remove cluster boundary */
		if(pCluster->boundaries == &pMerge->boundary[0])
			pCluster->boundaries = pCluster->boundaries->next;
		pMerge->boundary[0].next->prev = pMerge->boundary[1].prev;
		pMerge->boundary[1].prev->next = pMerge->boundary[0].next;
		pMerge->boundary[0].prev->next = pMerge->boundary[1].next;
		pMerge->boundary[1].next->prev = pMerge->boundary[0].prev;
		for(i=0; i<2; ++i)
		{
			pBoundary = pMerge->boundary[i].next;
			pBoundary2 = pMerge->boundary[1-i].prev;
			if(pBoundary->twin->cluster == pBoundary2->twin->cluster)
			{
				pBoundary->merge->length += pBoundary2->merge->length;
				GI_LIST_REMOVE(pCluster->boundaries, pBoundary2);
				if(pBoundary->twin->cluster)
				{
					GI_LIST_REMOVE(pBoundary2->twin->cluster
						->boundaries, pBoundary2->twin);
				}
				GIHeap_remove(&qFringe, pBoundary2->merge);
				GI_FREE_SINGLE(pBoundary2->merge, sizeof(GIClusterMerge));
			}
		}
		GI_FREE_SINGLE(pMerge, sizeof(GIClusterMerge));

		/* adjust costs of all affected merges */
		GIHash_clear(&hNeighbours, 0);
		GI_LIST_FOREACH(pCluster->boundaries, pBoundary)
			if(pBoundary->twin->cluster)
				GIHash_insert(&hNeighbours, pBoundary->twin->cluster, NULL);
		GI_LIST_NEXT(pCluster->boundaries, pBoundary)
		GI_LIST_FOREACH(pCluster->boundaries, pBoundary)
			pCluster2 = pBoundary->twin->cluster;
			dError = 0.0;
			GI_LIST_FOREACH(pCluster2->boundaries, pBoundary2)
				if((pBoundary2 != pBoundary->twin->prev || 
					pBoundary2->twin->cluster != pBoundary->next->twin->cluster) && 
					(pBoundary2 != pBoundary->twin->next || 
					pBoundary2->twin->cluster != pBoundary->prev->twin->cluster) && 
					GIHash_find(&hNeighbours, pBoundary2->twin->cluster))
				{
					dError = DBL_MAX;
					break;
				}
			GI_LIST_NEXT(pCluster2->boundaries, pBoundary2)
			if(dError == 0.0)
				dError = GIClusterMerge_error(
					pBoundary->merge, orientation, shape);
			GIHeap_enqueue(&qFringe, pBoundary->merge, dError);
		GI_LIST_NEXT(pCluster->boundaries, pBoundary)
	}
	printf("cluster merging end\n");

	/* clean up */
	for(i=1; i<=qFringe.count; ++i)
	{
		GI_FREE_SINGLE(qFringe.items[i].data, sizeof(GIClusterMerge));
	}
	return uiClusters;
}

/** \internal
 *  \brief Compute merging error.
 *  \param merge merge operation to compute error for
 *  \param orientation orientation bias weight
 *  \param shape shape bias weight
 *  \return merging error or \c DBL_MAX if invalid
 *  \ingroup mrm
 */
GIdouble GIClusterMerge_error(GIClusterMerge *merge, 
							  GIdouble orientation, GIdouble shape)
{
	static const GIdouble dFourPI = 4.0 * GI_PI;
	GIFaceCluster *pCluster0 = merge->boundary[0].cluster, 
		*pCluster1 = merge->boundary[1].cluster;
	GIdouble qdFit[10], qdDir[10], Z[9], eig[3], tmp[2];
	GIdouble normal[3], n[3], d;
	GIdouble dError, dArea = pCluster0->area + pCluster1->area, dInvC;
	GIuint uiFCount;
	if(!pCluster0 || !pCluster1)
		return DBL_MAX;

	/* compute average normal and fit quadric */
	uiFCount = pCluster0->fcount + pCluster1->fcount;
	GI_VEC3_SCALE(normal, pCluster0->normal, 
		(GIdouble)(pCluster0->fcount));
	GI_VEC3_ADD_SCALED(normal, normal, pCluster1->normal, 
		(GIdouble)(pCluster1->fcount));
	GI_QUADRIC_ADD(qdFit, pCluster0->P, pCluster1->P);

	/* compute best fitting plane and fit error */
	dInvC = 1.0 / qdFit[9];
	Z[0] = qdFit[0] - qdFit[6]*qdFit[6]*dInvC;
	Z[1] = qdFit[1] - qdFit[6]*qdFit[7]*dInvC;
	Z[2] = qdFit[2] - qdFit[6]*qdFit[8]*dInvC;
	Z[3] = Z[1];
	Z[4] = qdFit[3] - qdFit[7]*qdFit[7]*dInvC;
	Z[5] = qdFit[4] - qdFit[7]*qdFit[8]*dInvC;
	Z[6] = Z[2];
	Z[7] = Z[5];
	Z[8] = qdFit[5] - qdFit[8]*qdFit[8]*dInvC;
	GImat3d_tridiagonalize(Z, eig, tmp);
	if(GImat3d_ql(Z, eig, tmp, 64) <= 64)
	{
		GIuint i = (eig[1] < eig[0]) ? 1 : 0;
		if(eig[2] < eig[i])
			i = 2;
		GI_VEC3_COPY(n, Z+3*i);
		if(GI_VEC3_DOT(n, normal) < 0.0)
		{
			GI_VEC3_SCALE(n, n, -1.0);
		}
	}
	else
	{
		GI_VEC3_COPY(n, normal);
	}
	GIvec3d_normalize(n);
	d = -GI_VEC3_DOT(n, qdFit+6) / qdFit[9];
	dError = GI_QUADRIC_PLANE_ERROR(qdFit, n, d) / (GIdouble)(3*uiFCount);

	/* compute orientation bias */
	if(orientation != 0.0)
	{
		GI_QUADRIC_ADD(qdDir, pCluster0->R, pCluster1->R);
		dError += orientation * (GI_QUADRIC_VECTOR_ERROR(qdDir, n)/dArea);
	}

	/* compute shape bias */
	if(shape != 0.0)
	{
		GIdouble dPerimeter = pCluster0->perimeter + 
			pCluster1->perimeter - 2.0*merge->length;
		GIdouble dGamma0 = (pCluster0->perimeter*pCluster0->perimeter) / 
			(dFourPI*pCluster0->area), dGamma1 = (pCluster1->perimeter*
			pCluster1->perimeter) / (dFourPI*pCluster1->area), dGamma = 
			(dPerimeter*dPerimeter) / (dFourPI*dArea);
		dError += shape * ((dGamma-GI_MAX(dGamma0, dGamma1))/dGamma);
	}
	return dError;
}
