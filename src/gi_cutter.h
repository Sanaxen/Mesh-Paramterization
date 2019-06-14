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
 *  \brief Declaration of structures and functions for mesh cutting.
 */

#ifndef __GI_CUTTER_H__
#define __GI_CUTTER_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_math.h"
#include "gi_mesh.h"
#include "gi_container.h"


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Cutter configuration.
 *  \ingroup cutting
 */
typedef struct _GICutter
{
	struct _GIContext	*context;					/**< Context this cutter belongs to */
	GIenum				cutter;						/**< Cutting algorithm to use. */
	GIboolean			straighten;					/**< Straighten cut paths. */
	GIuint				iterations;					/**< Number of subdivision iterations. */
	GIfloat				orientation_weight;			/**< Weight of orientation bias for face clustering. */
	GIfloat				shape_weight;				/**< Weight of shape bias for face clustering. */
} GICutter;

/** \internal
 *  \brief Mesh patch.
 *  \details This structure represents a patch as a subset of faces.
 *  \ingroup cutting
 */
typedef struct _GIPatch
{
	GIuint				id;							/**< Patch ID. */
	GIMesh				*mesh;						/**< Mesh this patch belongs to. */
	GIuint				fcount;						/**< Number of faces in the patch. */
	GIuint				pcount;						/**< Number of parameter coordinates. */
	GIuint				hcount;						/**< Number of half edges on the cut. */
	GIuint				path_count;					/**< Number of cut paths. */
	GIuint				groups;						/**< Number of undirected paths. */
	GIFace				*faces;						/**< Start of face sublist. */
	GIParam				*params;					/**< List of params. */
	GIParam				*corners[4];				/**< Parameter coordinates of corners. */
	struct _GICutPath	*paths;						/**< List of cut paths. */
	GIdouble			hlength;					/**< Total length of the cut's half edges. */
	GIdouble			side_lengths[4];			/**< Half edge lengths of one side (if corners fixed). */
	GIuint				resolution;					/**< Resolution of border params in texels. */
	GIdouble			surface_area;				/**< Area of patch in 3D. */
	GIdouble			param_area;					/**< Area of patch in parameter space. */
	GIdouble			stretch[GI_STRETCH_COUNT];	/**< Stretch values for all metrics. */
	GIdouble			min_param_stretch;			/**< Minimum param stretch value. */
	GIdouble			max_param_stretch;			/**< Maximum param stretch value. */
	GIenum				param_metric;				/**< Current stretch metric for param stretch. */
	GIboolean			parameterized;				/**< Patch has valid parameterization. */
	GIboolean			fixed_corners;				/**< Patch corners are permanent. */
	GIDynamicQueue		split_paths;				/**< Stack of path splits. */
	struct _GIPatch		*next;						/**< Next patch (for convenience). */
} GIPatch;

/** \internal
 *  \brief Cut path.
 *  \details This structure represents a directed path on the closed cut connecting two cut nodes.
 *  \ingroup cutting
 */
typedef struct _GICutPath
{
	GIuint				id;							/**< ID of cut path. */
	GIPatch				*patch;						/**< Patch this cut path belongs to. */
	GIuint				group;						/**< Undirected patch ID. */
	GIdouble			elength;					/**< Total length of path. */
	GIint				glength;					/**< Length of path in texels. */
	GIParam				*pstart;					/**< Param, this path starts at. */
	struct _GICutPath	*twin;						/**< Opposite path. */
	struct _GICutPath	*next;						/**< Next path in list. */
	struct _GICutPath	*prev;						/**< Previous path in list. */
} GICutPath;


/*************************************************************************/
/* Functions */

/** \name Cutter methods
 *  \{
 */
void GICutter_construct(GICutter *cutter, struct _GIContext *context);
GIint GICutter_from_params(GICutter *cutter, GIMesh *mesh);
GIint GICutter_initial_gim(GICutter *cutter, GIMesh *mesh);
GIint GICutter_face_clustering(GICutter *cutter, GIMesh *mesh);
GIint GICutter_seamless_atlas(GICutter *cutter, GIMesh *mesh);
GIint GICutter_catmull_clark(GICutter *cutter, GIMesh *mesh);
GIboolean GICutter_create_params(GICutter *cutter, 
	GIMesh *mesh, GIint *patch_ids, GIubyte *cut_flags);
void GICutter_straighten_paths(GICutter *cutter, 
	GIMesh *mesh, GIVertex *root, GIubyte *cut_flags);
void GICutter_sort_faces(GICutter *cutter, 
	GIMesh *mesh, GIint *patch_ids, GIubyte *cut_flags);
void GICutter_compute_paths(GICutter *cutter, GIMesh *mesh);
/** \} */

/** \name Patch methods
 *  \{
 */
void GIPatch_destruct(GIPatch *patch);
void GIPatch_prevent_singularities(GIPatch *patch);
void GIPatch_renumerate_params(GIPatch *patch);
GIboolean GIPatch_find_corners(GIPatch *patch);
GIboolean GIPatch_valid_parameterization(GIPatch *patch);
GIboolean GIPatch_split_path(GIPatch *patch, GICutPath *path, GIParam *param);
void GIPatch_revert_splits(GIPatch *patch, GIint count);
GIFace* GIPatch_compute_stretch(GIPatch *patch, GIuint metric, 
	GIboolean param_stretches, GIboolean init_stretches);
/** \} */


#endif
