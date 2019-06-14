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
 *  \brief Declaration of structures and functions for context handling.
 */

#ifndef __GI_CONTEXT_H__
#define __GI_CONTEXT_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_container.h"
#include "gi_mesh.h"
#include "gi_image.h"
#include "gi_cutter.h"
#include "gi_parameterizer.h"
#include "gi_sampler.h"
#include "gi_gl.h"

//#define OPENGI_DEBUG_OUTPUT

#ifdef OPENGI_DEBUG_OUTPUT
	#define GIDebug(a)	a
#else
	#define GIDebug(a)
#endif


/*************************************************************************/
/* Macros */

#define GI_VERSION_MAJOR		2
#define GI_VERSION_MINOR		1
#define GI_VERSION_REVISION		1

#define GI_SUBSET_BASE			GI_EXACT_MAPPING_SUBSET
#define GI_SUBSET_END			GI_PARAM_CORNER_SUBSET
#define GI_SUBSET_COUNT			(GI_SUBSET_END-GI_SUBSET_BASE+1)


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief GI context.
 *  \details This structure manages all GI states.
 *  \ingroup context
 */
typedef struct _GIContext
{
	GIHash			mesh_hash;							/**< Hash of meshes. */
	GIHash			image_hash;							/**< Hash of images. */
	GIMesh			*mesh;								/**< Current bound mesh. */
//	GIImage			*image;								/**< Current bound image. */
//	GIImage			*attrib_image[GI_ATTRIB_COUNT];		/**< Attribute images. */
	GIuint			next_mid;							/**< ID of next created mesh. */
	GIuint			next_iid;							/**< ID of next created image. */
	GIboolean		use_threads;						/**< Use multithreading. */
	GIfloat			*attrib_pointer[GI_ATTRIB_COUNT];	/**< Attribute pointers. */
	GIsizei			attrib_size[GI_ATTRIB_COUNT];		/**< Size values for attribute pointers. */
	GIsizei			attrib_stride[GI_ATTRIB_COUNT];		/**< Stride values for attribute pointers. */
	GIenum			attrib_normalized[GI_ATTRIB_COUNT];	/**< Normalization flag fo attributes. */
	GIboolean		attrib_enabled[GI_ATTRIB_COUNT];	/**< Enabled flags for attribute pointers. */
	GIenum			attrib_semantic[GI_ATTRIB_COUNT];	/**< Semantics of attributes. */
	GIuint			semantic[GI_SEMANTIC_COUNT];		/**< Attribute semantics. */
	const GIuint	*subset[GI_SUBSET_COUNT];			/**< Vertex subsets. */
	GIsizei			subset_count[GI_SUBSET_COUNT];		/**< Numbers of elements in vertex subsets. */
	GIboolean		subset_sorted[GI_SUBSET_COUNT];		/**< Sorted flags for vertex subsets. */
	GIboolean		subset_enabled[GI_SUBSET_COUNT];	/**< Enabled flags for vertex subsets. */
	GIenum			error;								/**< Error code of last encountered error. */
	GIerrorcb		error_cb;							/**< Error callback function. */
	GIvoid			*edata;								/**< User data for error callback. */
	GICutter		cutter;								/**< Cutting state. */
	GIParameterizer	parameterizer;						/**< Parameterizer state. */
//	GISampler		sampler;							/**< Sampler state. */
	GIRenderer		renderer;							/**< Render state. */
	GIGLManager		*gl_manager;						/**< OpenGL manager. */
} GIContext;


/*************************************************************************/
/* Functions */

//#ifdef _DEBUG
	GIContext* GIContext_current();
//#else
//	extern GIContext *gi_CurrentContext;
//	#define GIContext_current()		gi_CurrentContext;
//#endif
void GIContext_error(GIContext *context, GIenum error);


#endif
