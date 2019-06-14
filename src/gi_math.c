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
 *  \brief Implementation of simple mathematic helper functions.
 */

#include "gi_math.h"

#include <string.h>

#define INFN		0x7F800000
#define MAXN		0x477FE000
#define MINN		0x38800000
#define SIGNN		0x80000000
#define INFC		0x0003FC00
#define NANN		0x7F802000
#define MAXC		0x00023BFF
#define MINC		0x0001C400
#define SIGNC		0x00008000
#define MULN		0x52000000
#define MULC		0x33800000
#define SUBC		0x000003FF
#define NORC		0x00000400
#define MAXD		0x0001C000
#define MIND		0x0001C000


/** \internal
 *  \brief Vertices of platonic solids.
 *  \ingroup math
 */
const GIfloat g_platonic_v[50][3] = 
{
	{ 0.0f, 0.0f, 1.0f }, { 0.942809f, 0, -0.333333f }, 
	{ -0.471405f, 0.816497f, -0.333333 }, { -0.471405f, -0.816497f, -0.333333 }, 

	{ -0.57735f, -0.57735f, -0.57735f }, { 0.57735f, -0.57735f, -0.57735f }, 
	{ 0.57735f, 0.57735f, -0.57735f }, { -0.57735f, 0.57735f, -0.57735f }, 
	{ -0.57735f, -0.57735f, 0.57735f }, { 0.57735f, -0.57735f, 0.57735f }, 
	{ 0.57735f, 0.57735f, 0.57735f }, { -0.57735f, 0.57735f, 0.57735f }, 

	{ 1.0f, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, 
	{ 0.0f, 1.0f, 0.0f }, { 0.0f, -1.0f, 0.0f }, 
	{ 0.0f, 0.0f, 1.0f }, { 0.0f, 0.0f, -1.0f }, 

	{ 0.57735f, 0.57735f, 0.57735f }, { 0.57735f, 0.57735f, -0.57735f }, 
	{ 0.57735f, -0.57735f, 0.57735f }, { 0.57735f, -0.57735f, -0.57735f }, 
	{ -0.57735f, 0.57735f, 0.57735f }, { -0.57735f, 0.57735f, -0.57735f }, 
	{ -0.57735f, -0.57735f, 0.57735f }, { -0.57735f, -0.57735f, -0.57735f }, 
	{ 0.356822f, 0.934172, 0.0f }, { -0.356822f, 0.934172, 0.0f }, 
	{ 0.356822f, -0.934172, 0.0f }, { -0.356822f, -0.934172, 0.0f }, 
	{ 0.934172, 0.0f, 0.356822f }, { 0.934172, 0.0f, -0.356822f }, 
	{ -0.934172, 0.0f, 0.356822f }, { -0.934172, 0.0f, -0.356822f }, 
	{ 0.0f, 0.356822f, 0.934172 }, { 0.0f, -0.356822f, 0.934172 }, 
	{ 0.0f, 0.356822f, -0.934172 }, { 0.0f, -0.356822f, -0.934172 }, 

	{ 0.85065f, 0.525731f, 0.0f }, { -0.85065f, 0.525731f, 0.0f }, 
	{ 0.85065f, -0.525731f, 0.0f }, { -0.85065f, -0.525731f, 0.0f }, 
	{ 0.525731f, 0.0f, 0.85065f }, { 0.525731f, 0.0f, -0.85065f }, 
	{ -0.525731f, 0.0f, 0.85065f }, { -0.525731f, 0.0f, -0.85065f }, 
	{ 0.0f, 0.85065f, 0.525731f }, { 0.0f, -0.85065f, 0.525731f }, 
	{ 0.0f, 0.85065f, -0.525731f }, { 0.0f, -0.85065f, -0.525731f }, 
};


/** \internal
 *  \brief Convert float to half.
 *  \details http://stackoverflow.com/questions/1659440/32-bit-to-16-bit-floating-point-conversion
 *  \param value single precision floating point value
 *  \return half precision floating point value
 *  \ingroup math
 */
GIhalf GIhalf_from_float(GIfloat value)
{
	GIint v = *(GIint*)&value, s = MULN;
	GIuint sign = v & SIGNN;
	v ^= sign;
	sign >>= 16;
	s = (*(GIfloat*)&s) * (*(GIfloat*)&v);
	v ^= (s^v) & -(MINN>v);
	v ^= (INFN^v) & -((INFN>v)&(v>MAXN));
	v ^= (NANN^v) & -((NANN>v)&(v>INFN));
	*(GIuint*)&v >>= 13;
	v ^= ((v-MAXD)^v) & -(v>MAXC);
	v ^= ((v-MIND)^v) & -(v>SUBC);
	return *(GIuint*)&v | sign;
}

/** \internal
 *  \brief Convert half to float.
 *  \details http://stackoverflow.com/questions/1659440/32-bit-to-16-bit-floating-point-conversion
 *  \param value half precision floating point value
 *  \return single precision floating point value
 *  \ingroup math
 */
GIfloat GIhalf_to_float(GIhalf value)
{
	GIint v = value, s = MULC;
	GIint sign = v & SIGNC, mask;
	v ^= sign;
	sign <<= 16;
	v ^= ((v+MIND)^v) & -(v>SUBC);
	v ^= ((v+MAXD)^v) & -(v>MAXC);
	*(GIfloat*)&s *= v;
	mask = -(NORC>v);
	v <<= 13;
	v ^= (s^v) & mask;
	v |= sign;
	return *(GIfloat*)&v;
}

/** \internal
 *  \brief Set 2D vector.
 *  \param d vector to set
 *  \param x x coordinate to set
 *  \param y y coordinate to set
 *  \ingroup math
 */
void GIvec2f_set(GIvec2f d, GIfloat x, GIfloat y)
{
	GI_VEC2_SET(d, x, y);
}

/** \internal
 *  \brief Copy 2D vector.
 *  \param d vector take copy
 *  \param s vector to copy from
 *  \ingroup math
 */
void GIvec2f_copy(GIvec2f d, const GIvec2f s)
{
	GI_VEC2_COPY(d, s);
}

/** \internal
 *  \brief compare 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \retval GI_TRUE if equal
 *  \retval GI_FALSE if not equal
 *  \ingroup math
 */
GIboolean GIvec2f_equal(const GIvec2f v, const GIvec2f w)
{
	return GI_VEC2_EQUAL(v, w);
}

/** \internal
 *  \brief Sum of 2D vectors.
 *  \param d vector to store sum \a v + \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2f_add(GIvec2f d, const GIvec2f v, const GIvec2f w)
{
	GI_VEC2_ADD(d, v, w);
}

/** \internal
 *  \brief Difference of 2D vectors.
 *  \param d vector to store difference \a v - \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2f_sub(GIvec2f d, const GIvec2f v, const GIvec2f w)
{
	GI_VEC2_SUB(d, v, w);
}

/** \internal
 *  \brief Scale 2D vector.
 *  \param d vector to store scaled vector
 *  \param v vector to scale
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec2f_scale(GIvec2f d, const GIvec2f v, GIfloat s)
{
	GI_VEC2_SCALE(d, v, s);
}

/** \internal
 *  \brief Add scaled 2D vector to other vector.
 *  \param d vector to store sum \a v + \a s*w
 *  \param v first vector
 *  \param w second vector (to be scaled)
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec2f_add_scaled(GIvec2f d, const GIvec2f v, const GIvec2f w, GIfloat s)
{
	GI_VEC2_ADD_SCALED(d, v, w, s);
}

/** \internal
 *  \brief Negate 2D vector.
 *  \param d vector to store negative vector
 *  \param s vector to negate
 *  \ingroup math
 */
void GIvec2f_negate(GIvec2f d, const GIvec2f s)
{
	GI_VEC2_NEGATE(d, s);
}

/** \internal
 *  \brief Round 2D vector.
 *  \param d vector to store rounded vector
 *  \param s vector to round
 *  \ingroup math
 */
void GIvec2f_round(GIvec2f d, const GIvec2f s)
{
	GI_VEC2_ROUND(d, s);
}

/** \internal
 *  \brief Component-wise minimum of two 2D vectors.
 *  \param d vector to store minimum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2f_min(GIvec2f d, const GIvec2f v, const GIvec2f w)
{
	GI_VEC2_MIN(d, v, w);
}

/** \internal
 *  \brief Component-wise maximum of two 2D vectors.
 *  \param d vector to store maximum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2f_max(GIvec2f d, const GIvec2f v, const GIvec2f w)
{
	GI_VEC2_MAX(d, v, w);
}

/** \internal
 *  \brief Length of 2D vector.
 *  \param v vector
 *  \return eucildean norm of vector sqrt(<\a v,\a v>)
 *  \ingroup math
 */
GIfloat GIvec2f_length(const GIvec2f v)
{
	return GI_VEC2_LENGTH(v);
}

/** \internal
 *  \brief Squared length of 2D vector.
 *  \param v vector
 *  \return squared euclidean norm of vector <\a v,\a v>
 *  \ingroup math
 */
GIfloat GIvec2f_length_sqr(const GIvec2f v)
{
	return GI_VEC2_LENGTH_SQR(v);
}

/** \internal
 *  \brief Distance between 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return euclidean distance between vectors
 *  \ingroup math
 */
GIfloat GIvec2f_dist(const GIvec2f v, const GIvec2f w)
{
	GIfloat vec[2];
	GI_VEC2_SUB(vec, w, v);
	return GI_VEC2_LENGTH(vec);
}

/** \internal
 *  \brief Squared distance between 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return squared euclidean distance between vectors
 *  \ingroup math
 */
GIfloat GIvec2f_dist_sqr(const GIvec2f v, const GIvec2f w)
{
	GIfloat vec[2];
	GI_VEC2_SUB(vec, w, v);
	return GI_VEC2_LENGTH_SQR(vec);
}

/** \internal
 *  \brief Normalize 2D vector.
 *  \param v vector to normalize
 *  \ingroup math
 */
void GIvec2f_normalize(GIvec2f v)
{
	GIfloat fInvNorm = 1.0f / GI_VEC2_LENGTH(v);
	v[0] *= fInvNorm; v[1] *= fInvNorm;
}

/** \internal
 *  \brief Dot product of 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return dot product <\a v,\a w>
 *  \ingroup math
 */
GIfloat GIvec2f_dot(const GIvec2f v, const GIvec2f w)
{
	return GI_VEC2_DOT(v, w);
}

/** \internal
 *  \brief Cross product equivalent of 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return dot product det(\a v,\a w)
 *  \ingroup math
 */
GIfloat GIvec2f_det(const GIvec2f v, const GIvec2f w)
{
	return GI_VEC2_DET(v, w);
}

/** \internal
 *  \brief Set 2D vector.
 *  \param d vector to set
 *  \param x x coordinate to set
 *  \param y y coordinate to set
 *  \ingroup math
 */
void GIvec2d_set(GIvec2d d, GIdouble x, GIdouble y)
{
	GI_VEC2_SET(d, x, y);
}

/** \internal
 *  \brief Copy 2D vector.
 *  \param d vector take copy
 *  \param s vector to copy from
 *  \ingroup math
 */
void GIvec2d_copy(GIvec2d d, const GIvec2d s)
{
	GI_VEC2_COPY(d, s);
}

/** \internal
 *  \brief compare 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \retval GI_TRUE if equal
 *  \retval GI_FALSE if not equal
 *  \ingroup math
 */
GIboolean GIvec2d_equal(const GIvec2d v, const GIvec2d w)
{
	return GI_VEC2_EQUAL(v, w);
}

/** \internal
 *  \brief Sum of 2D vectors.
 *  \param d vector to store sum \a v + \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2d_add(GIvec2d d, const GIvec2d v, const GIvec2d w)
{
	GI_VEC2_ADD(d, v, w);
}

/** \internal
 *  \brief Difference of 2D vectors.
 *  \param d vector to store difference \a v - \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2d_sub(GIvec2d d, const GIvec2d v, const GIvec2d w)
{
	GI_VEC2_SUB(d, v, w);
}

/** \internal
 *  \brief Scale 2D vector.
 *  \param d vector to store scaled vector
 *  \param v vector to scale
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec2d_scale(GIvec2d d, const GIvec2d v, GIdouble s)
{
	GI_VEC2_SCALE(d, v, s);
}

/** \internal
 *  \brief Add scaled 2D vector to other vector.
 *  \param d vector to store sum \a v + \a s*w
 *  \param v first vector
 *  \param w second vector (to be scaled)
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec2d_add_scaled(GIvec2d d, const GIvec2d v, const GIvec2d w, GIdouble s)
{
	GI_VEC2_ADD_SCALED(d, v, w, s);
}

/** \internal
 *  \brief Negate 2D vector.
 *  \param d vector to store negative vector
 *  \param s vector to negate
 *  \ingroup math
 */
void GIvec2d_negate(GIvec2d d, const GIvec2d s)
{
	GI_VEC2_NEGATE(d, s);
}

/** \internal
 *  \brief Round 2D vector.
 *  \param d vector to store rounded vector
 *  \param s vector to round
 *  \ingroup math
 */
void GIvec2d_round(GIvec2d d, const GIvec2d s)
{
	GI_VEC2_ROUND(d, s);
}

/** \internal
 *  \brief Component-wise minimum of two 2D vectors.
 *  \param d vector to store minimum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2d_min(GIvec2d d, const GIvec2d v, const GIvec2d w)
{
	GI_VEC2_MIN(d, v, w);
}

/** \internal
 *  \brief Component-wise maximum of two 2D vectors.
 *  \param d vector to store maximum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec2d_max(GIvec2d d, const GIvec2d v, const GIvec2d w)
{
	GI_VEC2_MAX(d, v, w);
}

/** \internal
 *  \brief Length of 2D vector.
 *  \param v vector
 *  \return eucildean norm of vector sqrt(<\a v,\a v>)
 *  \ingroup math
 */
GIdouble GIvec2d_length(const GIvec2d v)
{
	return GI_VEC2_LENGTH(v);
}

/** \internal
 *  \brief Squared length of 2D vector.
 *  \param v vector
 *  \return squared euclidean norm of vector <\a v,\a v>
 *  \ingroup math
 */
GIdouble GIvec2d_length_sqr(const GIvec2d v)
{
	return GI_VEC2_LENGTH_SQR(v);
}

/** \internal
 *  \brief Distance between 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return Euclidean distance between vectors
 *  \ingroup math
 */
GIdouble GIvec2d_dist(const GIvec2d v, const GIvec2d w)
{
	GIdouble vec[2];
	GI_VEC2_SUB(vec, w, v);
	return GI_VEC2_LENGTH(vec);
}

/** \internal
 *  \brief Squared distance between 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return squared euclidean distance between vectors
 *  \ingroup math
 */
GIdouble GIvec2d_dist_sqr(const GIvec2d v, const GIvec2d w)
{
	GIdouble vec[2];
	GI_VEC2_SUB(vec, w, v);
	return GI_VEC2_LENGTH_SQR(vec);
}

/** \internal
 *  \brief Normalize 2D vector.
 *  \param v vector to normalize
 *  \ingroup math
 */
void GIvec2d_normalize(GIvec2d v)
{
	GIdouble dInvNorm = 1.0f / GI_VEC2_LENGTH(v);
	v[0] *= dInvNorm; v[1] *= dInvNorm;
}

/** \internal
 *  \brief Dot product of 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return dot product <\a v,\a w>
 *  \ingroup math
 */
GIdouble GIvec2d_dot(const GIvec2d v, const GIvec2d w)
{
	return GI_VEC2_DOT(v, w);
}

/** \internal
 *  \brief Cross product equivalent of 2D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return dot product det(\a v,\a w)
 *  \ingroup math
 */
GIdouble GIvec2d_det(const GIvec2d v, const GIvec2d w)
{
	return GI_VEC2_DET(v, w);
}

/** \internal
 *  \brief Set 3D vector.
 *  \param d vector to set
 *  \param x x coordinate to set
 *  \param y y coordinate to set
 *  \param z z coordinate to set
 *  \ingroup math
 */
void GIvec3f_set(GIvec3f d, GIfloat x, GIfloat y, GIfloat z)
{
	GI_VEC3_SET(d, x, y, z);
}

/** \internal
 *  \brief Copy 3D vector.
 *  \param d vector take copy
 *  \param s vector to copy from
 *  \ingroup math
 */
void GIvec3f_copy(GIvec3f d, const GIvec3f s)
{
	GI_VEC3_COPY(d, s);
}

/** \internal
 *  \brief compare 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \retval GI_TRUE if equal
 *  \retval GI_FALSE if not equal
 *  \ingroup math
 */
GIboolean GIvec3f_equal(const GIvec3f v, const GIvec3f w)
{
	return GI_VEC3_EQUAL(v, w);
}

/** \internal
 *  \brief Sum of 3D vectors.
 *  \param d vector to store sum \a v + \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3f_add(GIvec3f d, const GIvec3f v, const GIvec3f w)
{
	GI_VEC3_ADD(d, v, w);
}

/** \internal
 *  \brief Difference of 3D vectors.
 *  \param d vector to store difference \a v - \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3f_sub(GIvec3f d, const GIvec3f v, const GIvec3f w)
{
	GI_VEC3_SUB(d, v, w);
}

/** \internal
 *  \brief Scale 3D vector.
 *  \param d vector to store scaled vector
 *  \param v vector to scale
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec3f_scale(GIvec3f d, const GIvec3f v, GIfloat s)
{
	GI_VEC3_SCALE(d, v, s);
}

/** \internal
 *  \brief Add scaled 3D vector to other vector.
 *  \param d vector to store sum \a v + \a s*w
 *  \param v first vector
 *  \param w second vector (to be scaled)
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec3f_add_scaled(GIvec3f d, const GIvec3f v, const GIvec3f w, GIfloat s)
{
	GI_VEC3_ADD_SCALED(d, v, w, s);
}

/** \internal
 *  \brief Negate 3D vector.
 *  \param d vector to store negative vector
 *  \param s vector to negate
 *  \ingroup math
 */
void GIvec3f_negate(GIvec3f d, const GIvec3f s)
{
	GI_VEC3_NEGATE(d, s);
}

/** \internal
 *  \brief Round 3D vector.
 *  \param d vector to store rounded vector
 *  \param s vector to round
 *  \ingroup math
 */
void GIvec3f_round(GIvec3f d, const GIvec3f s)
{
	GI_VEC3_ROUND(d, s);
}

/** \internal
 *  \brief Component-wise minimum of two 3D vectors.
 *  \param d vector to store minimum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3f_min(GIvec3f d, const GIvec3f v, const GIvec3f w)
{
	GI_VEC3_MIN(d, v, w);
}

/** \internal
 *  \brief Component-wise maximum of two 3D vectors.
 *  \param d vector to store maximum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3f_max(GIvec3f d, const GIvec3f v, const GIvec3f w)
{
	GI_VEC3_MAX(d, v, w);
}

/** \internal
 *  \brief Length of 3D vector.
 *  \param v vector
 *  \return eucildean norm of vector sqrt(<\a v,\a v>)
 *  \ingroup math
 */
GIfloat GIvec3f_length(const GIvec3f v)
{
	return GI_VEC3_LENGTH(v);
}

/** \internal
 *  \brief Squared length of 3D vector.
 *  \param v vector
 *  \return squared euclidean norm of vector <\a v,\a v>
 *  \ingroup math
 */
GIfloat GIvec3f_length_sqr(const GIvec3f v)
{
	return GI_VEC3_LENGTH_SQR(v);
}

/** \internal
 *  \brief Distance between 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return Euclidean distance between vectors
 *  \ingroup math
 */
GIfloat GIvec3f_dist(const GIvec3f v, const GIvec3f w)
{
	GIfloat vec[3];
	GI_VEC3_SUB(vec, w, v);
	return GI_VEC3_LENGTH(vec);
}

/** \internal
 *  \brief Squared distance between 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return squared euclidean distance between vectors
 *  \ingroup math
 */
GIfloat GIvec3f_dist_sqr(const GIvec3f v, const GIvec3f w)
{
	GIfloat vec[3];
	GI_VEC3_SUB(vec, w, v);
	return GI_VEC3_LENGTH_SQR(vec);
}

/** \internal
 *  \brief Normalize 3D vector.
 *  \param v vector to normalize
 *  \ingroup math
 */
void GIvec3f_normalize(GIvec3f v)
{
	GIfloat fInvNorm = 1.0f / GI_VEC3_LENGTH(v);
	v[0] *= fInvNorm; v[1] *= fInvNorm; v[2] *= fInvNorm;
}

/** \internal
 *  \brief Dot product of 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return dot product <\a v,\a w>
 *  \ingroup math
 */
GIfloat GIvec3f_dot(const GIvec3f v, const GIvec3f w)
{
	return GI_VEC3_DOT(v, w);
}

/** \internal
 *  \brief Cross product of 3D vectors.
 *  \param d vector to store cross product \a v x \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3f_cross(GIvec3f d, const GIvec3f v, const GIvec3f w)
{
	GI_VEC3_CROSS(d, v, w);
}

/** \internal
 *  \brief Set 3D vector.
 *  \param d vector to set
 *  \param x x coordinate to set
 *  \param y y coordinate to set
 *  \param z z coordinate to set
 *  \ingroup math
 */
void GIvec3d_set(GIvec3d d, GIdouble x, GIdouble y, GIdouble z)
{
	GI_VEC3_SET(d, x, y, z);
}

/** \internal
 *  \brief Copy 3D vector.
 *  \param d vector take copy
 *  \param s vector to copy from
 *  \ingroup math
 */
void GIvec3d_copy(GIvec3d d, const GIvec3d s)
{
	GI_VEC3_COPY(d, s);
}

/** \internal
 *  \brief compare 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \retval GI_TRUE if equal
 *  \retval GI_FALSE if not equal
 *  \ingroup math
 */
GIboolean GIvec3d_equal(const GIvec3d v, const GIvec3d w)
{
	return GI_VEC3_EQUAL(v, w);
}

/** \internal
 *  \brief Sum of 3D vectors.
 *  \param d vector to store sum \a v + \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3d_add(GIvec3d d, const GIvec3d v, const GIvec3d w)
{
	GI_VEC3_ADD(d, v, w);
}

/** \internal
 *  \brief Difference of 3D vectors.
 *  \param d vector to store difference \a v - \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3d_sub(GIvec3d d, const GIvec3d v, const GIvec3d w)
{
	GI_VEC3_SUB(d, v, w);
}

/** \internal
 *  \brief Scale 3D vector.
 *  \param d vector to store scaled vector
 *  \param v vector to scale
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec3d_scale(GIvec3d d, const GIvec3d v, GIdouble s)
{
	GI_VEC3_SCALE(d, v, s);
}

/** \internal
 *  \brief Add scaled 3D vector to other vector.
 *  \param d vector to store sum \a v + \a s*w
 *  \param v first vector
 *  \param w second vector (to be scaled)
 *  \param s scalar value
 *  \ingroup math
 */
void GIvec3d_add_scaled(GIvec3d d, const GIvec3d v, const GIvec3d w, GIdouble s)
{
	GI_VEC3_ADD_SCALED(d, v, w, s);
}

/** \internal
 *  \brief Negate 3D vector.
 *  \param d vector to store negative vector
 *  \param s vector to negate
 *  \ingroup math
 */
void GIvec3d_negate(GIvec3d d, const GIvec3d s)
{
	GI_VEC3_NEGATE(d, s);
}

/** \internal
 *  \brief Round 3D vector.
 *  \param d vector to store rounded vector
 *  \param s vector to round
 *  \ingroup math
 */
void GIvec3d_round(GIvec3d d, const GIvec3d s)
{
	GI_VEC3_ROUND(d, s);
}

/** \internal
 *  \brief Component-wise minimum of two 3D vectors.
 *  \param d vector to store minimum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3d_min(GIvec3d d, const GIvec3d v, const GIvec3d w)
{
	GI_VEC3_MIN(d, v, w);
}

/** \internal
 *  \brief Component-wise maximum of two 3D vectors.
 *  \param d vector to store maximum vector
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3d_max(GIvec3d d, const GIvec3d v, const GIvec3d w)
{
	GI_VEC3_MAX(d, v, w);
}

/** \internal
 *  \brief Length of 3D vector.
 *  \param v vector
 *  \return eucildean norm of vector sqrt(<\a v,\a v>)
 *  \ingroup math
 */
GIdouble GIvec3d_length(const GIvec3d v)
{
	return GI_VEC3_LENGTH(v);
}

/** \internal
 *  \brief Squared length of 3D vector.
 *  \param v vector
 *  \return squared euclidean norm of vector <\a v,\a v>
 *  \ingroup math
 */
GIdouble GIvec3d_length_sqr(const GIvec3d v)
{
	return GI_VEC3_LENGTH_SQR(v);
}

/** \internal
 *  \brief Distance between 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return Euclidean distance between vectors
 *  \ingroup math
 */
GIdouble GIvec3d_dist(const GIvec3d v, const GIvec3d w)
{
	GIdouble vec[3];
	GI_VEC3_SUB(vec, w, v);
	return GI_VEC3_LENGTH(vec);
}

/** \internal
 *  \brief Squared distance between 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return squared euclidean distance between vectors
 *  \ingroup math
 */
GIdouble GIvec3d_dist_sqr(const GIvec3d v, const GIvec3d w)
{
	GIdouble vec[3];
	GI_VEC3_SUB(vec, w, v);
	return GI_VEC3_LENGTH_SQR(vec);
}

/** \internal
 *  \brief Normalize 3D vector.
 *  \param v vector to normalize
 *  \ingroup math
 */
void GIvec3d_normalize(GIvec3d v)
{
	GIdouble dInvNorm = 1.0f / GI_VEC3_LENGTH(v);
	v[0] *= dInvNorm; v[1] *= dInvNorm; v[2] *= dInvNorm;
}

/** \internal
 *  \brief Dot product of 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \return dot product <\a v,\a w>
 *  \ingroup math
 */
GIdouble GIvec3d_dot(const GIvec3d v, const GIvec3d w)
{
	return GI_VEC3_DOT(v, w);
}

/** \internal
 *  \brief Cross product of 3D vectors.
 *  \param d vector to store cross product \a v x \a w
 *  \param v first vector
 *  \param w second vector
 *  \ingroup math
 */
void GIvec3d_cross(GIvec3d d, const GIvec3d v, const GIvec3d w)
{
	GI_VEC3_CROSS(d, v, w);
}
