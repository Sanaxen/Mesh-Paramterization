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
 *  \brief Implementation of structures and functions for image handling.
 */

#include "gi_image.h"
#include "gi_context.h"
#include "gi_memory.h"

#include <stdlib.h>
#include <string.h>
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


/** Create new image object.
 *  \return id of new image
 *  \ingroup image
 */
GIuint GIAPIENTRY giGenImage()
{
	GIContext *pContext = GIContext_current();

	/* create image and add to hash */
	GIImage *pImage = (GIImage*)GI_CALLOC_SINGLE(sizeof(GIImage));
	pImage->id = pContext->next_iid++;
	GIHash_insert(&pContext->image_hash, &pImage->id, pImage);
	return pImage->id;
}

/** Create new image objects.
 *  \param n number of images to create
 *  \param images array to store ids of newly created images
 *  \ingroup image
 */
void GIAPIENTRY giGenImages(GIsizei n, GIuint *images)
{
	GIint i;

	/* create images */
	for(i=0; i<n; ++i)
		images[i] = giGenImage();
}

/** Bind image to a specified attribute type.
 *  \param image id of image to bind
 *  \ingroup image
 */
void GIAPIENTRY giBindImage(GIuint image)
{
	GIContext *pContext = GIContext_current();
	GIImage *pImage;

	/* image allready bound */
	if(pContext->image && pContext->image->id == image)
		return;
	if(!image)
	{
		pContext->image = NULL;
		return;
	}

	/* search and bind image */
	pImage = (GIImage*)GIHash_find(&pContext->image_hash, &image);
	if(!pImage)
		GIContext_error(pContext, GI_INVALID_ID);
	else
		pContext->image = pImage;
}

/** Delete image and all associated data.
 *  \param image id of image to delete
 *  \ingroup image
 */
void GIAPIENTRY giDeleteImage(GIuint image)
{
	GIContext *pContext = GIContext_current();
	GIImage *pImage;
	GIuint a;

	/* find image */
	pImage = (GIImage*)GIHash_remove(&pContext->image_hash, &image);
	if(!pImage)
	{
		GIContext_error(pContext, GI_INVALID_ID);
		return;
	}

	/* remove and clean up */
	if(pImage == pContext->image)
			pContext->image = NULL;
	for(a=0; a<GI_ATTRIB_COUNT; ++a)
		if(pImage == pContext->attrib_image[a])
			pContext->attrib_image[a] = NULL;
	if(image == pContext->next_iid-1)
		--pContext->next_iid;
	GI_FREE_SINGLE(pImage, sizeof(GIImage));
}

/** Delete image objects.
 *  \param n number of images to delete
 *  \param images ids of images to delete
 *  \ingroup image
 */
void GIAPIENTRY giDeleteImages(GIsizei n, const GIuint *images)
{
	GIint i;

	/* delete images */
	for(i=0; i<n; ++i)
		giDeleteImage(images[i]);
}

/** Set image to work on external data.
 *  This function tells the image to store its data at an external memory address.
 *  The caller is responsible for allocating and releasing the neccessary data.
 *  A call to this function destroys the effect of a previous call to giImage***Data().
 *  \param width width of image
 *  \param height height of image
 *  \param components number of components (3 or 4)
 *  \param type data type of image elements
 *  \param data address of appropriately sized image data in user-owned memory
 *  \ingroup image
 */
void GIAPIENTRY giImageExternalData(GIsizei width, GIsizei height, 
									GIsizei components, GIenum type, GIvoid *data)
{
	/* set data */
	GIImage_set_data(width, height, components, type, data, 0, 0);
}

/** Set image to work on OpenGL texture.
 *  This function tells the image to store its data in an OpenGL texture object.
 *  The caller is responsible for creating and deleting the texture object.
 *  A call to this function destroys the effect of a previous call to giImage***Data().
 *  \param width width of image
 *  \param height height of image
 *  \param components number of components (3 or 4)
 *  \param type data type of image elements
 *  \param texture name of appropriately sized user-owned OpenGL texture object
 *  \ingroup image
 */
void GIAPIENTRY giImageGLTextureData(GIsizei width, GIsizei height, 
									 GIsizei components, GIenum type, GIuint texture)
{
	/* set data */
	GIImage_set_data(width, height, components, type, 0, texture, 0);
}

/** Set image to work on OpenGL buffer.
 *  This function tells the image to store its data in an OpenGL buffer object.
 *  The caller is responsible for creating and deleting the buffer object.
 *  A call to this function destroys the effect of a previous call to giImage***Data().
 *  \param width width of image
 *  \param height height of image
 *  \param components number of components (3 or 4)
 *  \param type data type of image elements
 *  \param buffer name of appropriately sized user-owned OpenGL buffer object
 *  \ingroup image
 */
void GIAPIENTRY giImageGLBufferData(GIsizei width, GIsizei height, 
									GIsizei components, GIenum type, GIuint buffer)
{
	/* set data */
	GIImage_set_data(width, height, components, type, 0, 0, buffer);
}

/** Set sub-image to work on.
 *  This function selects a sub-rectangle of the currently bound 
 *  image on which every following image operation works.
 *  \param x x-offset
 *  \param y y-offset
 *  \param width width of sub-rectangle
 *  \param height height of sub-rectangle
 *  \ingroup image
 */
void GIAPIENTRY giSubImage(GIuint x, GIuint y, GIsizei width, GIsizei height)
{
	GIContext *pContext = GIContext_current();
	GIImage *pImage;

	/* error checking */
	pImage = pContext->image;
	if(!pImage)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}
	if(!width)
		width = pImage->width;
	if(!height)
		height = pImage->height;
	if(width < 2 || height < 2 || 
		x+width > pImage->width || y+height > pImage->height)
	{
		GIContext_error(pContext, GI_INVALID_VALUE);
		return;
	}

	/* set data */
	pImage->offset_x = x;
	pImage->offset_y = y;
	pImage->sub_width = width;
	pImage->sub_height = height;
	if(pImage->data)
		pImage->sub_data = (GIbyte*)pImage->data + pImage->pixel_size*
			(pImage->offset_y*pImage->width+pImage->offset_x);
}

/** Retrieve property of current attribute image.
 *  \param pname parameter to querry
 *  \param param address of variable to store value
 *  \ingroup image
 */
void GIAPIENTRY giGetImageiv(GIenum pname, GIint *param)
{
	GIContext *pContext = GIContext_current();
	GIImage *pImage;

	/* image bound */
	pImage = pContext->image;
	if(!pImage)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}

	/* select parameter and get value */
	switch(pname)
	{
	case GI_IMAGE_WIDTH:
		*param = pImage->width;
		break;
	case GI_IMAGE_HEIGHT:
		*param = pImage->height;
		break;
	case GI_IMAGE_COMPONENTS:
		*param = pImage->comp;
		break;
	case GI_IMAGE_TYPE:
		*param = pImage->type;
		break;
	case GI_GL_IMAGE_TEXTURE:
		*param = pImage->texture;
		break;
	case GI_GL_IMAGE_BUFFER:
		*param = pImage->buffer;
		break;
	case GI_IMAGE_STORAGE:
		if(pImage->data)
			*param = GI_EXTERNAL_DATA;
		else if(pImage->texture)
			*param = GI_GL_TEXTURE_DATA;
		else if(pImage->buffer)
			*param = GI_GL_BUFFER_DATA;
		else
			*param = GI_NO_IMAGE_DATA;
		break;
	case GI_SUBIMAGE_X:
		*param = pImage->offset_x;
		break;
	case GI_SUBIMAGE_Y:
		*param = pImage->offset_y;
		break;
	case GI_SUBIMAGE_WIDTH:
		*param = pImage->sub_width;
		break;
	case GI_SUBIMAGE_HEIGHT:
		*param = pImage->sub_height;
		break;
	case GI_SUBIMAGE:
		param[0] = pImage->offset_x;
		param[1] = pImage->offset_y;
		param[2] = pImage->sub_width;
		param[3] = pImage->sub_height;
		break;
	default:
		GIContext_error(pContext, GI_INVALID_ENUM);
	}
}

/** \internal
 *  \brief Set working data for image.
 *  \param width width of image
 *  \param height height of image
 *  \param components number of components (3 or 4)
 *  \param type data type of image elements
 *  \param data address of appropriately sized image data in user-owned memory
 *  \param texture name of appropriately sized user-owned OpenGL texture object
 *  \param buffer name of appropriately sized user-owned OpenGL buffer object
 *  \ingroup image
 */
void GIImage_set_data(GIsizei width, GIsizei height, GIsizei components, 
					  GIenum type, GIvoid *data, GIuint texture, GIuint buffer)
{
	GIContext *pContext = GIContext_current();
	GIImage *pImage;

	/* correct attribute type and image bound */
	pImage = pContext->image;
	if(!pImage)
	{
		GIContext_error(pContext, GI_INVALID_OPERATION);
		return;
	}

	/* reset image and check for errors */
	memset((char*)pImage+sizeof(GIuint), 0, sizeof(GIImage)-sizeof(GIuint));
	if(type != GI_UNSIGNED_BYTE && type != GI_FLOAT && type != GI_HALF_FLOAT)
		GIContext_error(pContext, GI_INVALID_ENUM);
	else if(width < 2 || height < 2 || components < 1 || components > 4)
		GIContext_error(pContext, GI_INVALID_VALUE);
	else if((data || texture || buffer) && width && height)
	{
		/* set data */
		pImage->data = pImage->sub_data = data;
		pImage->texture = texture;
		pImage->buffer = buffer;
		pImage->sub_width = pImage->width = width;
		pImage->sub_height = pImage->height = height;
		pImage->offset_x = pImage->offset_y = 0;
		pImage->comp = components;
		pImage->type = type;
		pImage->pixel_size = components;
		if(type == GI_HALF_FLOAT)
			pImage->pixel_size *= sizeof(GIhalf);
		else if(type == GI_FLOAT)
			pImage->pixel_size *= sizeof(GIfloat);
		pImage->size = width * height * pImage->pixel_size;
		switch(components)
		{
		case 1:
			pImage->gl_format = GL_RED;
			if(type == GI_FLOAT)
				pImage->gl_internal = GL_R32F;
			else if(type == GI_HALF_FLOAT)
				pImage->gl_internal = GL_R16F;
			else
				pImage->gl_internal = GL_R8;
			break;
		case 2:
			pImage->gl_format = GL_RG;
			if(type == GI_FLOAT)
				pImage->gl_internal = GL_RG32F;
			else if(type == GI_HALF_FLOAT)
				pImage->gl_internal = GL_RG16F;
			else
				pImage->gl_internal = GL_RG8;
			break;
		case 3:
			pImage->gl_format = GL_RGB;
			if(type == GI_FLOAT)
				pImage->gl_internal = GL_RGB32F;
			else if(type == GI_HALF_FLOAT)
				pImage->gl_internal = GL_RGB16F;
			else
				pImage->gl_internal = GL_RGB8;
			break;
		case 4:
			pImage->gl_format = GL_RGBA;
			if(type == GI_FLOAT)
				pImage->gl_internal = GL_RGBA32F;
			else if(type == GI_HALF_FLOAT)
				pImage->gl_internal = GL_RGBA16F;
			else
				pImage->gl_internal = GL_RGBA8;
		}
	}
}

/** \internal
 *  \brief Set pixel in unsigned byte image.
 *  \param image image to modify
 *  \param x x coordinate to write to
 *  \param y y coordinate to write to
 *  \param vec data to write
 *  \ingroup image
 */
void GIImage_setpixel_ubyte(GIImage *image, GIuint x, 
							GIuint y, const GIfloat *vec)
{
	GIubyte *pData = (GIubyte*)image->sub_data + 
		(y*image->width+x)*image->comp;

	/* set pixel */
#if OPENGI_SSE >= 2
	__m128 XMM0 = _mm_load_ps(vec);
	__m128 XMM1 = _mm_set1_ps(255.0f);
	__m128i XMM8;
	XMM0 = _mm_mul_ps(XMM0, XMM1);
	XMM8 = _mm_cvtps_epi32(XMM0);
	XMM8 = _mm_packs_epi32(XMM8, XMM8);
	XMM8 = _mm_packus_epi16(XMM8, XMM8);
	if(image->comp == 4)
		*(GIuint*)pData = _mm_cvtsi128_si32(XMM8);
	else if(image->comp == 3)
	{
		GIuint i = _mm_cvtsi128_si32(XMM8);
		*(GIushort*)pData = i;
		*((GIubyte*)pData+2) = *((GIubyte*)&i+2);
	}
	else if(image->comp == 2)
		*(GIushort*)pData = _mm_cvtsi128_si32(XMM8);
	else
		*pData = _mm_cvtsi128_si32(XMM8);
#else
	GIsizei i;
	for(i=0; i<image->comp; ++i)
	{
		if(vec[i] > 1.0f)
			pData[i] = 255;
		else if(vec[i] < 0.0f)
			pData[i] = 0;
		else
			pData[i] = (GIubyte)(255.0f*vec[i]+0.5f);
	}
#endif
}

/** \internal
 *  \brief Set pixel in half precision floating point image.
 *  \param image image to modify
 *  \param x x coordinate to write to
 *  \param y y coordinate to write to
 *  \param vec data to write
 *  \ingroup image
 */
void GIImage_setpixel_half(GIImage *image, GIuint x, 
						   GIuint y, const GIfloat *vec)
{
	GIhalf *pData = (GIhalf*)image->sub_data + 
		(y*image->width+x)*image->comp;

	/* set pixel */
#if OPENGI_SSE >= 2
	static const GIint infN = 0x7F800000;
	static const GIint maxN = 0x477FE000;
	static const GIint minN = 0x38800000;
	static const GIint signN = 0x80000000;
	static const GIint nanN = 0x7F802000;
	static const GIint maxC = 0x00023BFF;
	static const GIint mulN = 0x52000000;
	static const GIint subC = 0x000003FF;
	static const GIint maxD = 0x0001C000;
	static const GIint minD = 0x0001C000;
	__m128i XMM2, XMM3, XMM4;
	__m128i XMM0 = _mm_load_si128((const __m128i*)vec);
	__m128i XMM1 = _mm_set1_epi32(signN);
	__m128 XMMS = _mm_set1_ps(*(GIfloat*)&mulN);
	XMM1 = _mm_and_si128(XMM1, XMM0);
	XMM0 = _mm_xor_si128(XMM0, XMM1);
	XMM1 = _mm_srli_epi32(XMM1, 16);
	XMMS = _mm_mul_ps(XMMS, *(__m128*)&XMM0);
	XMM2 = _mm_cvttps_epi32(XMMS);
	XMM2 = _mm_xor_si128(XMM2, XMM0);
	XMM3 = _mm_set1_epi32(minN);
	XMM3 = _mm_cmpgt_epi32(XMM3, XMM0);
	XMM2 = _mm_and_si128(XMM2, XMM3);
	XMM0 = _mm_xor_si128(XMM0, XMM2);
	XMM2 = _mm_set1_epi32(maxN);
	XMM2 = _mm_cmpgt_epi32(XMM0, XMM2);
	XMM3 = _mm_set1_epi32(infN);
	XMM4 = _mm_cmpgt_epi32(XMM3, XMM0);
	XMM2 = _mm_and_si128(XMM2, XMM4);
	XMM4 = _mm_xor_si128(XMM0, XMM3);
	XMM2 = _mm_and_si128(XMM2, XMM4);
	XMM0 = _mm_xor_si128(XMM0, XMM2);
	XMM3 = _mm_cmpgt_epi32(XMM0, XMM3);
	XMM4 = _mm_set1_epi32(nanN);
	XMM2 = _mm_cmpgt_epi32(XMM4, XMM0);
	XMM2 = _mm_and_si128(XMM2, XMM3);
	XMM4 = _mm_xor_si128(XMM4, XMM0);
	XMM2 = _mm_and_si128(XMM2, XMM4);
	XMM0 = _mm_xor_si128(XMM0, XMM2);
	XMM0 = _mm_srli_epi32(XMM0, 13);
	XMM2 = _mm_set1_epi32(maxC);
	XMM2 = _mm_cmpgt_epi32(XMM0, XMM2);
	XMM3 = _mm_set1_epi32(maxD);
	XMM3 = _mm_sub_epi32(XMM0, XMM3);
	XMM3 = _mm_xor_si128(XMM3, XMM0);
	XMM2 = _mm_and_si128(XMM2, XMM3);
	XMM0 = _mm_xor_si128(XMM0, XMM2);
	XMM2 = _mm_set1_epi32(subC);
	XMM2 = _mm_cmpgt_epi32(XMM0, XMM2);
	XMM3 = _mm_set1_epi32(minD);
	XMM3 = _mm_sub_epi32(XMM0, XMM3);
	XMM3 = _mm_xor_si128(XMM3, XMM0);
	XMM2 = _mm_and_si128(XMM2, XMM3);
	XMM0 = _mm_xor_si128(XMM0, XMM2);
	XMM0 = _mm_or_si128(XMM0, XMM1);
	if(image->comp == 1)
		*pData = _mm_cvtsi128_si32(XMM0);
	else
	{
		XMM0 = _mm_shufflelo_epi16(XMM0, _MM_SHUFFLE(3, 1, 2, 0));
		*(GIuint*)pData = _mm_cvtsi128_si32(XMM0);
		if(image->comp > 2)
		{
			XMM0 = _mm_shufflehi_epi16(XMM0, _MM_SHUFFLE(3, 1, 2, 0));
			XMM0 = _mm_shuffle_epi32(XMM0, _MM_SHUFFLE(2, 2, 2, 2));
			if(image->comp == 3)
				*(pData+2) = _mm_cvtsi128_si32(XMM0);
			else
				*((GIuint*)pData+1) = _mm_cvtsi128_si32(XMM0);
		}
	}
#else
	GIsizei i;
	for(i=0; i<image->comp; ++i)
		pData[i] = GIhalf_from_float(vec[i]);
#endif
}

/** \internal
 *  \brief Set pixel in single precision floating point image.
 *  \param image image to modify
 *  \param x x coordinate to write to
 *  \param y y coordinate to write to
 *  \param vec data to write
 *  \ingroup image
 */
void GIImage_setpixel_float(GIImage *image, GIuint x, 
							GIuint y, const GIfloat *vec)
{
	GIfloat *pData = (GIfloat*)image->sub_data + 
		(y*image->width+x)*image->comp;

	/* set pixel */
#if OPENGI_SSE >= 1
	if(image->comp == 4)
	{
		__m128 XMM0 = _mm_load_ps(vec);
		_mm_storeu_ps(pData, XMM0);
	}
	else
#endif
	{
		GIsizei i;
		for(i=0; i<image->comp; ++i)
			pData[i] = vec[i];
	}
}

/** \internal
 *  \brief Set pixel in 4-component single precision floating point image.
 *  \param image image to modify
 *  \param x x coordinate to write to
 *  \param y y coordinate to write to
 *  \param vec data to write
 *  \ingroup image
 */
void GIImage_setpixel_float4(GIImage *image, GIuint x, 
							 GIuint y, const GIfloat *vec)
{
	GIfloat *pData = (GIfloat*)image->sub_data + 
		((y*image->width+x)<<2);

	/* set pixel */
#if OPENGI_SSE >= 1
	__m128 XMM0 = _mm_load_ps(vec);
	_mm_store_ps(pData, XMM0);
#else
	pData[0] = vec[0];
	pData[1] = vec[1];
	pData[2] = vec[2];
	pData[3] = vec[3];
#endif
}

/** \internal
 *  \brief Get pixel on image border.
 *  \param image image to work on
 *  \param p position on border (ccw from lower left)
 *  \return address of border pixel
 *  \ingroup image
 */
GIvoid* GIImage_border_pixel(GIImage *image, GIuint p)
{
	static const GIint corners[8] = { 0, 0, 1, 0, 1, 1, 0, 1 };
	static const GIint dirs[8] = { 1, 0, 0, 1, -1, 0, 0, -1 };
	GIuint uiWidth1 = image->sub_width-1, 
		uiHeight1 = image->sub_height-1, uiSide = 0;

	/* find side and position on side */
	if(p >= uiWidth1)
	{
		p -= uiWidth1;
		uiSide += 2;
		if(p >= uiHeight1)
		{
			p -= uiHeight1;
			uiSide += 2;
			if(p >= uiWidth1)
			{
				p -= uiWidth1;
				uiSide += 2;
			}
		}
	}

	/* get pixel */
	return (GIbyte*)image->data + image->pixel_size*((corners[uiSide+1]*
		uiHeight1 + p*dirs[uiSide+1]+image->offset_y)*image->width+
		corners[uiSide]*uiWidth1 + p*dirs[uiSide] + image->offset_x);
}
#endif
