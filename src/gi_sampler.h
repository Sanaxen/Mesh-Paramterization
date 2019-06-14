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
 *  \brief Declaration of functions for geometry image sampling.
 */

#ifndef __GI_SAMPLER_H__
#define __GI_SAMPLER_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_cutter.h"
#include "gi_image.h"
#include "gi_thread.h"
#include "gi_gl.h"


/*************************************************************************/
/* Macros */

#define GI_MAX_GROUP			4


/*************************************************************************/
/* Typedefs */

/** \internal
 *  \brief Texture filtering function.
 *  \ingroup sampling
 */
typedef void (*GIfilterfunc)(struct _GITexture*, const GIfloat*, GIfloat*);


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Sampler state.
 *  \ingroup sampling
 */
typedef struct _GISampler
{
	struct _GIContext	*context;							/**< Context this sampler belongs to */
	GIenum				sampler;							/**< Sampling mode. */
	GIboolean			use_shader;							/**< Sample using GLSL if possible. */
	GIboolean			use_fbo;							/**< Sample using FBOs if possible. */
	GIuint				sampled_attribs;					/**< Sampled attributes from last call to giSample(). */
	GIenum				attrib_mode[GI_ATTRIB_COUNT];		/**< Sampling modes of attributes. */
	GIuint				attrib_texture[GI_ATTRIB_COUNT];	/**< GL texture to resample. */
	GIuint				texture_dim[GI_ATTRIB_COUNT];		/**< Texture dimension. */
	GIfloat				attrib_matrix[GI_ATTRIB_COUNT][16];	/**< Custom transformation matrix. */
	GIbitfield			identity_matrix;					/**< Identity flags for attribute matrices. */
	GIsamplercb			callback[GI_ATTRIB_COUNT];			/**< Callback functions. */
	GIvoid				*cdata[GI_ATTRIB_COUNT];			/**< Data for callback functions. */
} GISampler;

/** \internal
 *  \brief 2D Texture.
 *  \details This structure represents the RAM-image of an 
 *  OpenGL-texture used for software sampling.
 *  \ingroup sampling
 */
typedef struct _GITexture
{
	GIfloat			*data;					/**< Texture data. */
	GIfloat			*image[6];				/**< Image offstets (for cube map). */
	GIuint			dimension;				/**< Texture dimension. */
	GIsizei			width[6];				/**< Width of texture. */
	GIsizei			height[6];				/**< Height of texture. */
	GIsizei			depth[6];				/**< Depth of texture. */
	GIint			mag_filter;				/**< Magnification filter. */
	GIint			min_filter;				/**< Minification filter. */
	GIint			wrap_s;					/**< Wrapping mode for S-coordinate. */
	GIint			wrap_t;					/**< Wrapping mode for T-coordinate. */
	GIint			wrap_r;					/**< Wrapping mode for R-coordinate. */
	GIenum			filter_mode;			/**< Filter to use. */
	GIfilterfunc	filter;					/**< Filtering method. */
} GITexture;

/** \internal
 *  \brief Information about group of attributes to sample.
 *  \ingroup sampling
 */
typedef struct _GIAttribGroup
{
	GIuint					count;					/**< Number of attributes in group. */
	GIuint					attrib[GI_MAX_GROUP];	/**< Attributes of group. */
	struct _GIAttribGroup	*prev;					/**< Previous group in list. */
	struct _GIAttribGroup	*next;					/**< Next group in list. */
} GIAttribGroup;

/** \internal
 *  \brief Parameters for software rasterizer.
 *  \ingroup sampling
 */
typedef struct _GIRasterizerData
{
	GISampler	*sampler;					/**< Sampler state. */
	GIPatch		*patch;						/**< Patch to work on. */
	GIFace		*next_face;					/**< Next face to rasterize. */
	GIFace		*end_face;					/**< Face after last face to rasterize. */
	GIuint		num_attribs;				/**< Number of attributes to sample. */
	GIuint		attribs[GI_ATTRIB_COUNT];	/**< Attributes to sample. */
	GIdouble	*params;					/**< Scaled params. */
	GITexture	*textures;					/**< Texture array. */
	GIpfunc		*setpixel_fn;				/**< Array of pixel setting functions. */
	GIMutex		mutex;						/**< Mutex for next face. */
} GIRasterizerData;


/*************************************************************************/
/* Functions */

/** \name Sampler methods.
 *  \{
 */
void GISampler_construct(GISampler *sampler, struct _GIContext *context);
void GISampler_sample_software(GISampler *sampler, GIPatch *patch);
void GISampler_sample_opengl(GISampler *sampler, GIPatch *patch);
void GISampler_sample_gl_fixed_function(GISampler *sampler, 
	GIPatch *patch, GIGLManager *gl, GIuint attrib);
void GISampler_sample_gl_shader(GISampler *sampler, 
	GIPatch *patch, GIGLManager *gl, GIuint attrib);
void GISampler_rasterize_triangle(GIFace *face, GIRasterizerData *data, GIfloat *temp);
GIthreadret GITHREADENTRY GISampler_rasterize_thread(void *arg);
/** \} */

/** \name Texture methods.
 *  \{
 */
void GITexture_construct(GITexture *texture, GIuint gl_texture, GIuint dim);
void GITexture_destruct(GITexture *texture);
const GIfloat* GITexture_getpixel1D(GITexture *texture, GIint x);
const GIfloat* GITexture_getpixel2D(GITexture *texture, GIint x, GIint y);
const GIfloat* GITexture_getpixel3D(GITexture *texture, GIint x, GIint y, GIint z);
const GIfloat* GITexture_getpixelCube(GITexture *texture, GIuint face, GIint x, GIint y);
void GITexture_filter1D(GITexture *texture, const GIfloat *coords, GIfloat *out);
void GITexture_filter2D(GITexture *texture, const GIfloat *coords, GIfloat *out);
void GITexture_filter3D(GITexture *texture, const GIfloat *coords, GIfloat *out);
void GITexture_filterCube(GITexture *texture, const GIfloat *coords, GIfloat *out);
/** \} */


#endif
#endif