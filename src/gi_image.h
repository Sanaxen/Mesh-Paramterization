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
 *  \brief Declaration of structures and functions for image handling.
 */

#ifndef __GI_IMAGE_H__
#define __GI_IMAGE_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Attribute image.
 *  \details This structure represents an image that can be bound to a specific attribute type.
 *  \ingroup image
 */
typedef struct _GIImage
{
	GIuint			id;						/**< Image ID. */
	GIsizei			width;					/**< Width of the image. */
	GIsizei			height;					/**< Height of the image. */
	GIsizei			offset_x;				/**< X-offset of sub-image. */
	GIsizei			offset_y;				/**< Y-offset of sub-image. */
	GIsizei			sub_width;				/**< Width of sub-image. */
	GIsizei			sub_height;				/**< Height of sub-image. */
	GIsizei			comp;					/**< Number of components in [1,4]. */
	GIenum			type;					/**< Data type of image elements. */
	GIuint			size;					/**< Size of image data in bytes. */
	GIuint			pixel_size;				/**< Size of a single pixel. */
	GIenum			gl_format;				/**< OpenGL texture format. */
	GIenum			gl_internal;			/**< OpenGL texture internal format. */
	GIvoid			*data;					/**< Data of the image (if RAM storage). */
	GIvoid			*sub_data;				/**< Pointer to sub-image (if RAM storage). */
	GIuint			texture;				/**< OpenGL texture object (if texture storage). */
	GIuint			buffer;					/**< OpenGL buffer object (if buffer storage). */
} GIImage;


/*************************************************************************/
/* Typedefs */

/** \internal
 *  \brief Pixel setting function.
 *  \ingroup image
 */
typedef void (*GIpfunc)(GIImage*, GIuint, GIuint, const GIfloat*);


/*************************************************************************/
/* Functions */

/** \name Image methods
 *  \{
 */
void GIImage_set_data(GIsizei width, GIsizei height, GIsizei components, 
	GIenum type, GIvoid *data, GIuint texture, GIuint buffer);
void GIImage_setpixel_ubyte(GIImage *image, GIuint x, GIuint y, const GIfloat *vec);
void GIImage_setpixel_half(GIImage *image, GIuint x, GIuint y, const GIfloat *vec);
void GIImage_setpixel_float(GIImage *image, GIuint x, GIuint y, const GIfloat *vec);
void GIImage_setpixel_float4(GIImage *image, GIuint x, GIuint y, const GIfloat *vec);
GIvoid* GIImage_border_pixel(GIImage *image, GIuint p);
/** \} */


#endif
#endif
