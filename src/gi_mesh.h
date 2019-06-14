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
 *  \brief Declaration of structures and functions for mesh handling.
 */

#ifndef __GI_MESH_H__
#define __GI_MESH_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_container.h"
#include "gi_math.h"

#define GI_ATTRIB_COUNT			16

#define GI_SEMANTIC_BASE		GI_POSITION_ATTRIB
#define GI_SEMANTIC_END			GI_PARAM_STRETCH_ATTRIB
#define GI_SEMANTIC_COUNT		(GI_SEMANTIC_END-GI_SEMANTIC_BASE+1)

#define GI_STRETCH_BASE			GI_MAX_GEOMETRIC_STRETCH
#define GI_STRETCH_END			GI_COMBINED_STRETCH
#define GI_STRETCH_COUNT		(GI_STRETCH_END-GI_STRETCH_BASE+1)

#define GI_VERTEX_EXACT_BIT		0x01
#define GI_VERTEX_CORNER_BIT	0x02


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Triangle mesh.
 *  \details This structure represents a triangular mesh.
 *  \ingroup mesh
 */
typedef struct _GIMesh
{
	GIuint				id;							/**< Mesh ID. */
	struct _GIContext	*context;					/**< Context the mesh belongs to. */
	GIuint				fcount;						/**< Number of faces. */
	GIuint				ecount;						/**< Number of edges. */
	GIuint				vcount;						/**< Number of vertices. */
	GIuint				pcount;						/**< Number of parameter coordinates. */
	GIuint				acount;						/**< Number of attributes. */
	struct _GIFace		*faces;						/**< Doubly linked list of faces. */
	struct _GIEdge		*edges;						/**< Doubly linked list of edges. */
	struct _GIVertex	*vertices;					/**< Doubly linked list of vertices. */
	struct _GIAttribute	*attributes;				/**< Doubly linked list of attributes. */
	GIuint				attrib_size;				/**< Size of attribute data. */
	GIint				aoffset[GI_ATTRIB_COUNT];	/**< Offsets of attributes in attribute data. */
	GIsizei				asize[GI_ATTRIB_COUNT];		/**< Number of components for attributes. */
	GIboolean			anorm[GI_ATTRIB_COUNT];		/**< Normalization flags of attributes. */
	GIenum				asemantic[GI_ATTRIB_COUNT];	/**< Semantics of attributes. */
	GIuint				semantic[GI_SEMANTIC_COUNT];/**< Attribute semantics. */
	GIuint				varray_size;				/**< Vertex array size of last retrieval. */
	GIuint				varray_attribs;				/**< Queried attributes of last retrieval. */
	struct _GIPatch		*varray_patch;				/**< Queried patch of last retrieval. */
	GIint				genus;						/**< Mesh genus. */
	GIdouble			aabb_min[3];				/**< Minimal coordinate values. */
	GIdouble			aabb_max[3];				/**< Maximal coordinate values. */
	GIdouble			radius;						/**< Radius of bounding sphere around origin. */
	GIdouble			mean_edge;					/**< Average edge length. */
	GIdouble			surface_area;				/**< Area of mesh in 3D. */
	GIdouble			param_area;					/**< Area of mesh in parameter space. */
	GIdouble			stretch[GI_STRETCH_COUNT];	/**< Stretch values for all metrics. */
	GIdouble			min_param_stretch;			/**< Minimum param stretch value. */
	GIdouble			max_param_stretch;			/**< Maximum param stretch value. */
	GIenum				param_metric;				/**< Current stretch metric for param stretch. */
	GIDynamicQueue		split_hedges;				/**< Stack of split half edges. */
	GIuint				patch_count;				/**< Number of patches. */
	GIuint				param_patches;				/**< Number of parameterized patches. */
	GIuint				resolution;					/**< Param resolution (if same for all patches) */
	GIuint				cut_splits;					/**< Number of splits before parameterization. */
	GIuint				pre_cut_splits;				/**< Number of splits before cut creation. */
	GIdouble			*old_coords;				/**< Original vertex coordinates. */
	struct _GIPatch		*patches;					/**< Patches. */
	struct _GIPatch		*active_patch;				/**< Currently selected patch. */
} GIMesh;

/** \internal
 *  \brief Mesh face.
 *  \details This structure represents a face in a triangular mesh.
 *  \ingroup mesh
 */
typedef struct _GIFace
{
	GIuint				id;						/**< Face ID. */
	struct _GIHalfEdge	*hedges;				/**< Half edge list. */
	struct _GIFace		*next;					/**< Next face in list. */
	struct _GIFace		*prev;					/**< Previous face in list. */
} GIFace;

/** \internal
 *  \brief Mesh half edge.
 *  \details This structure represents a half edge of a face.
 *  \ingroup mesh
 */
typedef struct _GIHalfEdge
{
	GIFace				*face;					/**< Face, this half edge belongs to or NULL for boundary. */
	struct _GIEdge		*edge;					/**< Edge, this half edge belongs to. */
	struct _GIVertex	*vstart;				/**< Vertex, this half edge starts at. */
	struct _GIAttribute	*astart;				/**< Attribute of starting vertex. */
	struct _GIParam		*pstart;				/**< Parameter coordinates of starting vertex. */
	struct _GIHalfEdge	*twin;					/**< Opposite half edge. */
	struct _GIHalfEdge	*next;					/**< Next half edge in face. */
	struct _GIHalfEdge	*prev;					/**< Previous half edge in face. */
} GIHalfEdge;

/** \internal
 *  \brief Mesh edge.
 *  \details This structure represents an edge in a triangular mesh.
 *  \ingroup mesh
 */
typedef struct _GIEdge
{
	GIuint			id;							/**< Edge ID. */
	GIdouble		length;						/**< Length of edge. */
	GIHalfEdge		hedge[2];					/**< Half edges. */
	struct _GIEdge	*next;						/**< Next edge in list. */
	struct _GIEdge	*prev;						/**< Previous edge in list. */
} GIEdge;

/** \internal
 *  \brief Mesh vertex.
 *  \details This structure represents a vertex of a triangular mesh.
 *  \ingroup mesh
 */
typedef struct _GIVertex
{
	GIuint				id;						/**< Vertex ID. */
	GIubyte				flags;					/**< Vertex properties. */
	GIubyte				cut_degree;				/**< Number of connected cut EDGES. */
	GIdouble			coords[3];				/**< Coordinates of this vertex. */
	GIHalfEdge			*hedge;					/**< Any half edge starting at this vertex. */
	struct _GIVertex	*next;					/**< Next vertex in list. */
	struct _GIVertex	*prev;					/**< Previous vertex in list. */
} GIVertex;

/** \internal
 *  \brief Vertex attribute.
 *  \details This structure represents an attribute associated 
 *  with the corner of a mesh face.
 *  \ingroup mesh
 */
typedef struct _GIAttribute
{
	GIuint				id;						/**< Attribute ID. */
	struct _GIAttribute	*next;					/**< Next attribute in list. */
	struct _GIAttribute *prev;					/**< Previous attribute in list. */
} GIAttribute;

/** \internal
 *  \brief Parameter coordinates.
 *  \details This structure represents a parameter coordinate 
 *  associated with the corner of a mesh face.
 *  \ingroup mesh
 */
typedef struct _GIParam
{
	GIuint				id;						/**< ID of parameter coordinate. */
	GIdouble			params[2];				/**< Parameter coordinates. */
	GIdouble			stretch;				/**< Per param stretch value. */
	GIVertex			*vertex;				/**< Vertex, this parameter coordinate belongs to. */
	GIHalfEdge			*cut_hedge;				/**< Cut half edge or NULL if parameter coord not on cut. */
	struct _GIParam		*next;					/**< Next parameter coordinate in list. */
	struct _GIParam		*prev;					/**< Previous parameter coordinate in list. */
} GIParam;


/*************************************************************************/
/* Functions */

/** \name Mesh methods
 *  \{
 */
void GIMesh_destruct(GIMesh *mesh);
void GIMesh_destroy_cut(GIMesh *mesh);
void GIMesh_revert_splits(GIMesh *mesh, GIint count);
GIint GIMesh_genus(GIMesh *mesh);
void GIMesh_compute_stretch(GIMesh *mesh, GIuint metric, GIboolean param_stretches);
/** \} */

/** \name Face methods
 *  \{
 */
GIHalfEdge* GIFace_halfedge_at(GIFace *face, GIuint i);
void GIFace_center(GIFace *face, GIdouble *center);
GIdouble GIFace_area(GIFace *face);
GIdouble GIFace_stretch(GIFace *face, GIenum metric, 
	GIvoid *args, GIdouble *area_2d);
/** \} */

/** \name Half edge methods
 *  \{
 */
GIuint GIHalfEdge_index(GIHalfEdge *hedge);
void GIHalfEdge_split(GIHalfEdge *hedge, struct _GIPatch *patch, 
	struct _GIPatch *twin_patch, GIdouble f, GIdouble *params);
void GIHalfEdge_half_split(GIHalfEdge *hedge, GIHalfEdge *hnew, 
	GIMesh *mesh, struct _GIPatch *patch, GIboolean hedge_first);
void GIHalfEdge_revert_half_split(GIHalfEdge *hedge, GIHalfEdge *hnew, 
	GIMesh *mesh, struct _GIPatch *patch, GIboolean hedge_first);
/** \} */

/** \name Attribute methods
 *  \{
 */
GIAttribute* GIAttribute_create_interpolated(GIAttribute *attrib1, 
	GIAttribute *attrib2, float f, GIMesh *mesh);

/** \name Param methods
 *  \{
 */
GIdouble GIParam_stretch(GIParam *param, GIenum metric, 
	GIvoid *args, GIHash *face_areas);
/** \} */


#endif
