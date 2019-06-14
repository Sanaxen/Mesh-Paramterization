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
 *  \brief Declaration of structures and functions for multiresolution modeling.
 */

#ifndef __GI_MULTIRESOLUTION_H__
#define __GI_MULTIRESOLUTION_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_mesh.h"


/*************************************************************************/
/* Macros */

/** \internal
 *  \brief Quadric constructor.
 *  \ingroup mrm
 */
#define GI_QUADRIC_CONSTRUCT(q,v,d)		(q)[0]=(v)[0]*(v)[0]; (q)[1]=(v)[0]*(v)[1]; \
										(q)[2]=(v)[0]*(v)[2]; (q)[3]=(v)[1]*(v)[1]; \
										(q)[4]=(v)[1]*(v)[2]; (q)[5]=(v)[2]*(v)[2]; \
										(q)[6]=(v)[0]*d; (q)[7]=(v)[1]*d; \
										(q)[8]=(v)[2]*d; (q)[9]=d*d;

/** \internal
 *  \brief Quadric error of vector.
 *  \ingroup mrm
 */
#define GI_QUADRIC_VECTOR_ERROR(q,v)	((q)[0]*(v)[0]*(v)[0] + \
										2.0*(q)[1]*(v)[0]*(v)[1] + \
										2.0*(q)[2]*(v)[0]*(v)[2] + \
										(q)[3]*(v)[1]*(v)[1] + \
										2.0*(q)[4]*(v)[1]*(v)[2] + \
										(q)[5]*(v)[2]*(v)[2] + \
										2.0*(q)[6]*(v)[0] + \
										2.0*(q)[7]*(v)[1] + \
										2.0*(q)[8]*(v)[2] + (q)[9])

/** \internal
 *  \brief Quadric error of plane.
 *  \ingroup mrm
 */
#define GI_QUADRIC_PLANE_ERROR(q,n,d)	((q)[0]*(n)[0]*(n)[0] + \
										2.0*(q)[1]*(n)[0]*(n)[1] + \
										2.0*(q)[2]*(n)[0]*(n)[2] + \
										(q)[3]*(n)[1]*(n)[1] + \
										2.0*(q)[4]*(n)[1]*(n)[2] + \
										(q)[5]*(n)[2]*(n)[2] + \
										2.0*(q)[6]*(n)[0]*(d) + \
										2.0*(q)[7]*(n)[1]*(d) + \
										2.0*(q)[8]*(n)[2]*(d) + \
										(q)[9]*(d)*(d))

/** \internal
 *  \brief Sum of two quadrics.
 *  \ingroup mrm
 */
#define GI_QUADRIC_ADD(d,q,p)			(d)[0]=(q)[0]+(p)[0]; (d)[1]=(q)[1]+(p)[1]; \
										(d)[2]=(q)[2]+(p)[2]; (d)[3]=(q)[3]+(p)[3]; \
										(d)[4]=(q)[4]+(p)[4]; (d)[5]=(q)[5]+(p)[5]; \
										(d)[6]=(q)[6]+(p)[6]; (d)[7]=(q)[7]+(p)[7]; \
										(d)[8]=(q)[8]+(p)[8]; (d)[9]=(q)[9]+(p)[9];

/** \internal
 *  \brief Scale a quadric.
 *  \ingroup mrm
 */
#define GI_QUADRIC_SCALE(d,q,s)			(d)[0]=(q)[0]*s; (d)[1]=(q)[1]*s; \
										(d)[2]=(q)[2]*s; (d)[3]=(q)[3]*s; \
										(d)[4]=(q)[4]*s; (d)[5]=(q)[5]*s; \
										(d)[6]=(q)[6]*s; (d)[7]=(q)[7]*s; \
										(d)[8]=(q)[8]*s; (d)[9]=(q)[9]*s;


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Face cluster.
 *  \details This structure represents a cluster of faces.
 *  \ingroup mrm
 */
typedef struct _GIFaceCluster
{
	GIuint						fcount;			/**< Number of faces in cluster. */
	GIFace						*faces;			/**< Start of face sublist. */
	struct _GIClusterBoundary	*boundaries;	/**< List of cluster boundaries. */
	GIdouble					P[10];			/**< Planarity quadric. */
	GIdouble					R[10];			/**< Orientation quadric. */
	GIdouble					area;			/**< Summed area of all faces. */
	GIdouble					perimeter;		/**< Perimeter of cluster. */
	GIdouble					normal[3];		/**< Average face normal. */
	struct _GIFaceCluster		*next;			/**< Next cluster in list. */
	struct _GIFaceCluster		*prev;			/**< Previous cluster in list. */
} GIFaceCluster;

/** \internal
 *  \brief Cluster boundary.
 *  \details This structure represents a half boundary of a face cluster.
 *  \ingroup mrm
 */
typedef struct _GIClusterBoundary
{
	GIFaceCluster				*cluster;		/**< Cluster this boundary belongs to. */
	struct _GIClusterMerge		*merge;			/**< Merge-edge this boundary belongs to. */
	struct _GIClusterBoundary	*twin;			/**< Opposite boundary. */
	struct _GIClusterBoundary	*next;			/**< Next boundary in list. */
	struct _GIClusterBoundary	*prev;			/**< Previous boundary in list. */
} GIClusterBoundary;

/** \internal
 *  \brief Cluster merge operation.
 *  \details this class represents a face cluster boundary to be merged away.
 *  \ingroup mrm
 */
typedef struct _GIClusterMerge
{
	GIdouble			length;					/**< Length of boundary. */
	GIClusterBoundary	boundary[2];			/**< Half boundaries. */
} GIClusterMerge;


/*************************************************************************/
/* Functions */

/** \name Multiresolution methods
 *  \{
 */
GIuint GIMultiresolution_face_clustering(GIMesh *mesh, 
	GIFaceCluster *clusters, GIdouble orientation, GIdouble shape, 
	GIdouble max_error, GIuint min_clusters, GIuint max_clusters);
/** \} */

/** \name Face cluster methods
 *  \{
 */
void GIFaceCluster_construct(GIFaceCluster *cluster, GIFace *face);
/** \} */

/** \name Cluster merge methods
 *  \{
 */
GIdouble GIClusterMerge_error(GIClusterMerge *merge, 
	GIdouble orientation, GIdouble shape);
/** \} */


#endif
