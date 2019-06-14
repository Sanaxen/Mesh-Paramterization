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
 *  \brief Declaration of simple mathematic helper functions.
 */

#ifndef __GI_MATH_H__
#define __GI_MATH_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include <stdlib.h>
#include <math.h>
#include <float.h>


/*************************************************************************/
/* Variables */

/* Vertices of platonic solids */
extern const GIfloat g_platonic_v[50][3];

/* Indices of platonic solids in array */
#define GI_TETRAHEDRON			0
#define GI_HEXAHEDRON			4
#define GI_OCTAHEDRON			12
#define GI_DODECAHEDRON			18
#define GI_ICOSAHEDRON			38


/*************************************************************************/
/* Typedefs */

/** \internal
 *  \brief single precision vectors and matrices.
 *  \ingroup math
 */
typedef GIfloat *GIvec2f, *GIvec3f, *GImat3f, *GImat4f;

/** \internal
 *  \brief double precision vectors and matrices.
 *  \ingroup math
 */
typedef GIdouble *GIvec2d, *GIvec3d, *GImat3d, *GImat4d;


/*************************************************************************/
/* Macros */
//#define MY_EPS	(1.0e-6)

/** \internal
 *  \brief Constant Pi.
 *  \ingroup math
 */
#define GI_PI						3.1415926535897932

/** \internal
 *  \brief Constant 2 * Pi.
 *  \ingroup math
 */
#define GI_TWO_PI					6.2831853071795865

/** \internal
 *  \brief Constant Pi / 2.
 *  \ingroup math
 */
#define GI_HALF_PI					1.5707963267948966

/** \internal
 *  \brief Constant 1 / Pi.
 *  \ingroup math
 */
#define GI_INV_PI					0.3183098861837907

/** \internal
 *  \brief Constant 1 / (2*Pi).
 *  \ingroup math
 */
#define GI_INV_TWO_PI				0.1591549430918953

/** \internal
 *  \brief Constant 2 / Pi.
 *  \ingroup math
 */
#define GI_INV_HALF_PI				0.6366197723675813

/** \internal
 *  \brief Convert radians to degrees.
 *  \ingroup math
 */
#define GI_RAD_TO_DEG(a)			((a)*57.295779513082321)

/** \internal
 *  \brief Convert degrees to radians.
 *  \ingroup math
 */
#define GI_DEG_TO_RAD(a)			((a)*0.0174532925199432)

/** \internal
 *  \brief Minimum of two values.
 *  \ingroup math
 */
#define GI_MIN(a,b)					(((a)<(b)) ? (a) : (b))

/** \internal
 *  \brief Maximum of two values.
 *  \ingroup math
 */
#define GI_MAX(a,b)					(((a)<(b)) ? (b) : (a))

/** \internal
 *  \brief Sign of floating point number.
 *  \ingroup math
 */
#define GI_SIGN(a)					(((a)>=0.0) ? 1.0 : -1.0)

/** \internal
 *  \brief Round floating point number.
 *  \ingroup math
 */
#define GI_ROUND(a)					(((a)-floor(a))<0.5 ? floor(a) : ceil(a))

/** \internal
 *  \brief Clamp number.
 *  \ingroup math
 */
#define GI_CLAMP(x,a,b)				GI_MAX(a, GI_MIN(x, b))

/** \internal
 *  \brief Swap arbitrary variables.
 *  \ingroup math
 */
#define GI_SWAP(a,b,t)				(t)=(a); (a)=(b); (b)=(t)

/** \internal
 *  \brief Check if power of two.
 *  \ingroup math
 */
#define GI_POWER_OF_2(x)			(!((x)&((x)-1)))

/** \internal
 *  \brief Check if not a power of two.
 *  \ingroup math
 */
#define GI_NON_POWER_OF_2(x)		((x)&((x)-1))

/** \internal
 *  \brief Round to next power of two.
 *  \ingroup math
 */
#define GI_TO_UPPER_POWER_OF_2(x)	--(x); (x)|=(x)>>1; (x)|=(x)>>2; (x)|=(x)>>4; \
									(x)|=(x)>>8; (x)|=(x)>>16; ++(x)

/** \internal
 *  \brief Set 2D vector.
 *  \ingroup math
 */
#define GI_VEC2_SET(d,x,y)			(d)[0]=x; (d)[1]=y;

/** \internal
 *  \brief Copy 2D vector.
 *  \ingroup math
 */
#define GI_VEC2_COPY(d,s)			(d)[0]=(s)[0]; (d)[1]=(s)[1];

/** \internal
 *  \brief Compare 2D vectors.
 *  \ingroup math
 */
#ifndef MY_EPS
#define GI_VEC2_EQUAL(v,w)			((v)[0]==(w)[0] && (v)[1]==(w)[1])
#else
#define GI_VEC2_EQUAL(v,w)			(fabs((v)[0]-(w)[0])<MY_EPS && fabs((v)[1]-(w)[1])<MY_EPS)
#endif

/** \internal
 *  \brief Sum of 2D vectors.
 *  \ingroup math
 */
#define GI_VEC2_ADD(d,v,w)			(d)[0]=(v)[0]+(w)[0]; (d)[1]=(v)[1]+(w)[1];

/** \internal
 *  \brief Difference of 2D vectors.
 *  \ingroup math
 */
#define GI_VEC2_SUB(d,v,w)			(d)[0]=(v)[0]-(w)[0]; (d)[1]=(v)[1]-(w)[1];

/** \internal
 *  \brief Scale 2D vector.
 *  \ingroup math
 */
#define GI_VEC2_SCALE(d,v,s)		(d)[0]=(v)[0]*(s); (d)[1]=(v)[1]*(s)

/** \internal
 *  \brief Add scaled 2D vector to other vector.
 *  \ingroup math
 */
#define GI_VEC2_ADD_SCALED(d,v,w,s)	(d)[0]=(v)[0]+(w)[0]*(s); (d)[1]=(v)[1]+(w)[1]*(s);

/** \internal
 *  \brief Negate 2D vector.
 *  \ingroup math
 */
#define GI_VEC2_NEGATE(d, s)		(d)[0]=-(s)[0]; (d)[1]=-(s)[1]

/** \internal
 *  \brief Round 2D vector.
 *  \ingroup math
 */
#define GI_VEC2_ROUND(d, s)			(d)[0]=GI_ROUND((s)[0]); (d)[1]=GI_ROUND((s)[1])

/** \internal
 *  \brief Minimum of two 2D vectors.
 *  \ingroup math
 */
#define GI_VEC2_MIN(d,v,w)			(d)[0]=GI_MIN((v)[0],(w)[0]); (d)[1]=GI_MIN((v)[1],(w)[1])

/** \internal
 *  \brief Maximum of two 2D vectors.
 *  \ingroup math
 */
#define GI_VEC2_MAX(d,v,w)			(d)[0]=GI_MAX((v)[0],(w)[0]); (d)[1]=GI_MAX((v)[1],(w)[1])

/** \internal
 *  \brief Dot product of 2D vectors.
 *  \ingroup math
 */
#define GI_VEC2_DOT(v,w)			((v)[0]*(w)[0]+(v)[1]*(w)[1])

/** \internal
 *  \brief Length of 2D vector.
 *  \ingroup math
 */
#define GI_VEC2_LENGTH_SQR(v)		GI_VEC2_DOT(v,v)

/** \internal
 *  \brief Squared length of 2D vector.
 *  \ingroup math
 */
#define GI_VEC2_LENGTH(v)			sqrt(GI_VEC2_DOT(v,v))

/** \internal
 *  \brief Cross product equivalent of 2D vectors
 *  \ingroup math
 */
#define GI_VEC2_DET(v,w)			((v)[0]*(w)[1]-(v)[1]*(w)[0])

/** \internal
 *  \brief Set 3D vector.
 *  \ingroup math
 */
#define GI_VEC3_SET(d,x,y,z)		(d)[0]=x; (d)[1]=y; (d)[2]=z

/** \internal
 *  \brief Copy 3D vector.
 *  \ingroup math
 */
#define GI_VEC3_COPY(d,s)			(d)[0]=(s)[0]; (d)[1]=(s)[1]; (d)[2]=(s)[2]

/** \internal
 *  \brief Compare 3D vectors.
 *  \ingroup math
 */
#ifndef MY_EPS
#define GI_VEC3_EQUAL(v,w)			((v)[0]==(w)[0] && (v)[1]==(w)[1] && (v)[2]==(w)[2])
#else
#define GI_VEC3_EQUAL(v,w)			(fabs((v)[0]-(w)[0])<MY_EPS && fabs((v)[1]-(w)[1])<MY_EPS && fabs((v)[2]-(w)[2])<MY_EPS)
#endif

/** \internal
 *  \brief Sum of 3D vectors.
 *  \ingroup math
 */
#define GI_VEC3_ADD(d,v,w)			(d)[0]=(v)[0]+(w)[0]; (d)[1]=(v)[1]+(w)[1]; (d)[2]=(v)[2]+(w)[2]

/** \internal
 *  \brief Difference of 3D vectors.
 *  \ingroup math
 */
#define GI_VEC3_SUB(d,v,w)			(d)[0]=(v)[0]-(w)[0]; (d)[1]=(v)[1]-(w)[1]; (d)[2]=(v)[2]-(w)[2]

/** \internal
 *  \brief Scale 3D vector.
 *  \ingroup math
 */
#define GI_VEC3_SCALE(d,v,s)		(d)[0]=(v)[0]*(s); (d)[1]=(v)[1]*(s); (d)[2]=(v)[2]*(s)

/** \internal
 *  \brief Add scaled 3D vector to other vector.
 *  \ingroup math
 */
#define GI_VEC3_ADD_SCALED(d,v,w,s)	(d)[0]=(v)[0]+(w)[0]*(s); (d)[1]=(v)[1]+(w)[1]*(s); (d)[2]=(v)[2]+(w)[2]*(s);

/** \internal
 *  \brief Negate 3D vector.
 *  \ingroup math
 */
#define GI_VEC3_NEGATE(d, s)		(d)[0]=-(s)[0]; (d)[1]=-(s)[1]; (d)[2]=-(s)[2];

/** \internal
 *  \brief Round 3D vector.
 *  \ingroup math
 */
#define GI_VEC3_ROUND(d, s)			(d)[0]=GI_ROUND((s)[0]); (d)[1]=GI_ROUND((s)[1]); (d)[2]=GI_ROUND((s)[2])

/** \internal
 *  \brief Minimum of two 3D vectors.
 *  \ingroup math
 */
#define GI_VEC3_MIN(d,v,w)			(d)[0]=GI_MIN((v)[0],(w)[0]); \
									(d)[1]=GI_MIN((v)[1],(w)[1]); \
									(d)[2]=GI_MIN((v)[2],(w)[2])

/** \internal
 *  \brief Maximum of two 3D vectors.
 *  \ingroup math
 */
#define GI_VEC3_MAX(d,v,w)			(d)[0]=GI_MAX((v)[0],(w)[0]); \
									(d)[1]=GI_MAX((v)[1],(w)[1]); \
									(d)[2]=GI_MAX((v)[2],(w)[2])

/** \internal
 *  \brief Dot product of 3D vectors.
 *  \ingroup math
 */
#define GI_VEC3_DOT(v,w)			((v)[0]*(w)[0]+(v)[1]*(w)[1]+(v)[2]*(w)[2])

/** \internal
 *  \brief Length of 3D vector.
 *  \ingroup math
 */
#define GI_VEC3_LENGTH_SQR(v)		GI_VEC3_DOT(v,v)

/** \internal
 *  \brief Squared length of 3D vector.
 *  \ingroup math
 */
#define GI_VEC3_LENGTH(v)			sqrt(GI_VEC3_DOT(v,v))

/** \internal
 *  \brief Cross product of 3D vectors.
 *  \ingroup math
 */
#define GI_VEC3_CROSS(d,v,w)		(d)[0]=(v)[1]*(w)[2]-(v)[2]*(w)[1]; \
									(d)[1]=(v)[2]*(w)[0]-(v)[0]*(w)[2]; \
									(d)[2]=(v)[0]*(w)[1]-(v)[1]*(w)[0]


/*************************************************************************/
/* Functions */

/** \name Half precision floating point methods
 *  \{
 */
GIhalf GIhalf_from_float(GIfloat value);
GIfloat GIhalf_to_float(GIhalf value);
/** \} */

/** \name Single precision 2D vector methods
 *  \{
 */
void GIvec2f_set(GIvec2f d, GIfloat x, GIfloat y);
void GIvec2f_copy(GIvec2f d, const GIvec2f s);
GIboolean GIvec2f_equal(const GIvec2f v, const GIvec2f w);
void GIvec2f_add(GIvec2f d, const GIvec2f v, const GIvec2f w);
void GIvec2f_sub(GIvec2f d, const GIvec2f v, const GIvec2f w);
void GIvec2f_scale(GIvec2f d, const GIvec2f v, GIfloat s);
void GIvec2f_add_scaled(GIvec2f d, const GIvec2f v, const GIvec2f w, GIfloat s);
void GIvec2f_negate(GIvec2f d, const GIvec2f s);
void GIvec2f_round(GIvec2f d, const GIvec2f s);
void GIvec2f_min(GIvec2f d, const GIvec2f v, const GIvec2f w);
void GIvec2f_max(GIvec2f d, const GIvec2f v, const GIvec2f w);
GIfloat GIvec2f_length(const GIvec2f v);
GIfloat GIvec2f_length_sqr(const GIvec2f v);
GIfloat GIvec2f_dist(const GIvec2f v, const GIvec2f w);
GIfloat GIvec2f_dist_sqr(const GIvec2f v, const GIvec2f w);
void GIvec2f_normalize(GIvec2f v);
GIfloat GIvec2f_dot(const GIvec2f v, const GIvec2f w);
GIfloat GIvec2f_det(const GIvec2f v, const GIvec2f w);
/** \} */

/** \name Double precision 2D vector methods
 *  \{
 */
void GIvec2d_set(GIvec2d d, GIdouble x, GIdouble y);
void GIvec2d_copy(GIvec2d d, const GIvec2d s);
GIboolean GIvec2d_equal(const GIvec2d v, const GIvec2d w);
void GIvec2d_add(GIvec2d d, const GIvec2d v, const GIvec2d w);
void GIvec2d_sub(GIvec2d d, const GIvec2d v, const GIvec2d w);
void GIvec2d_scale(GIvec2d d, const GIvec2d v, GIdouble s);
void GIvec2d_add_scaled(GIvec2d d, const GIvec2d v, const GIvec2d w, GIdouble s);
void GIvec2d_negate(GIvec2d d, const GIvec2d s);
void GIvec2d_round(GIvec2d d, const GIvec2d s);
void GIvec2d_min(GIvec2d d, const GIvec2d v, const GIvec2d w);
void GIvec2d_max(GIvec2d d, const GIvec2d v, const GIvec2d w);
GIdouble GIvec2d_length(const GIvec2d v);
GIdouble GIvec2d_length_sqr(const GIvec2d v);
GIdouble GIvec2d_dist(const GIvec2d v, const GIvec2d w);
GIdouble GIvec2d_dist_sqr(const GIvec2d v, const GIvec2d w);
void GIvec2d_normalize(GIvec2d v);
GIdouble GIvec2d_dot(const GIvec2d v, const GIvec2d w);
GIdouble GIvec2d_det(const GIvec2d v, const GIvec2d w);
/** \} */

/** \name Single precision 3D vector methods
 *  \{
 */
void GIvec3f_set(GIvec3f d, GIfloat x, GIfloat y, GIfloat z);
void GIvec3f_copy(GIvec3f d, const GIvec3f s);
GIboolean GIvec3f_equal(const GIvec3f v, const GIvec3f w);
void GIvec3f_add(GIvec3f d, const GIvec3f v, const GIvec3f w);
void GIvec3f_sub(GIvec3f d, const GIvec3f v, const GIvec3f w);
void GIvec3f_scale(GIvec3f d, const GIvec3f v, GIfloat s);
void GIvec3f_add_scaled(GIvec3f d, const GIvec3f v, const GIvec3f w, GIfloat s);
void GIvec3f_negate(GIvec3f d, const GIvec3f s);
void GIvec3f_round(GIvec3f d, const GIvec3f s);
void GIvec3f_min(GIvec3f d, const GIvec3f v, const GIvec3f w);
void GIvec3f_max(GIvec3f d, const GIvec3f v, const GIvec3f w);
GIfloat GIvec3f_length(const GIvec3f v);
GIfloat GIvec3f_length_sqr(const GIvec3f v);
GIfloat GIvec3f_dist(const GIvec3f v, const GIvec3f w);
GIfloat GIvec3f_dist_sqr(const GIvec3f v, const GIvec3f w);
void GIvec3f_normalize(GIvec3f v);
GIfloat GIvec3f_dot(const GIvec3f v, const GIvec3f w);
void GIvec3f_cross(GIvec3f d, const GIvec3f v, const GIvec3f w);
/** \{ */

/** \name Double precision 3D vector methods
 *  \{
 */
void GIvec3d_set(GIvec3d d, GIdouble x, GIdouble y, GIdouble z);
void GIvec3d_copy(GIvec3d d, const GIvec3d s);
GIboolean GIvec3d_equal(const GIvec3d v, const GIvec3d w);
void GIvec3d_add(GIvec3d d, const GIvec3d v, const GIvec3d w);
void GIvec3d_sub(GIvec3d d, const GIvec3d v, const GIvec3d w);
void GIvec3d_scale(GIvec3d d, const GIvec3d v, GIdouble s);
void GIvec3d_add_scaled(GIvec3d d, const GIvec3d v, const GIvec3d w, GIdouble s);
void GIvec3d_negate(GIvec3d d, const GIvec3d s);
void GIvec3d_round(GIvec3d d, const GIvec3d s);
void GIvec3d_min(GIvec3d d, const GIvec3d v, const GIvec3d w);
void GIvec3d_max(GIvec3d d, const GIvec3d v, const GIvec3d w);
GIdouble GIvec3d_length(const GIvec3d v);
GIdouble GIvec3d_length_sqr(const GIvec3d v);
GIdouble GIvec3d_dist(const GIvec3d v, const GIvec3d w);
GIdouble GIvec3d_dist_sqr(const GIvec3d v, const GIvec3d w);
void GIvec3d_normalize(GIvec3d v);
GIdouble GIvec3d_dot(const GIvec3d v, const GIvec3d w);
void GIvec3d_cross(GIvec3d d, const GIvec3d v, const GIvec3d w);
/** \{ */


#endif
