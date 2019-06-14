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
 *  \brief Implementation of common BLAS functions.
 */

#include "gi_blas.h"

#include <math.h>
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


/** \internal
 *  \brief Compute Givens rotation.
 *  \param a parameter a, r on output
 *  \param b parameter b, z on output
 *  \param c address to store c cosine of rotation at
 *  \param s address to store s sine of rotation at
 *  \ingroup numerics
 */
void drotg(MYdouble *a, MYdouble *b, MYdouble *c, MYdouble *s)
{
	MYdouble absa = fabs(*a), absb = fabs(*b), scale = absa + absb;
	if(scale == 0.0)
	{
		*c = 1.0;
		*s = *a = *b = 0.0;
	}
	else
	{
		MYdouble r, as = *a / scale, bs = *b / scale;
		r = scale*sqrt(as*as+bs*bs);
		if(((absa>absb) ? *a : *b) < 0.0)
			r = -r;
		*c = *a / r;
		*s = *b / r;
		*a = r;
		if(absa > absb)
			*b = *s;
		else if(*c != 0.0)
			*b = 1.0 / *c;
		else
			*b = 1.0;
	}
}

/** \internal
 *  \brief Rotate vectors by Givens rotation.
 *  \param n number of elements
 *  \param x first vector
 *  \param incx stride for vector x
 *  \param y second vector
 *  \param incy stride for vector y
 *  \param c cosine of rotation
 *  \param s sine of rotation
 *  \ingroup numerics
 */
void drot(int n, MYdouble *x, int incx, MYdouble *y, int incy, MYdouble c, MYdouble s)
{
	int i = 0;
	if(n <= 0)
		return;
	else if(incx == 1 && incy == 1)
	{
#if OPENGI_SSE >= 2
		__m128d XMM0 = _mm_set1_pd(c);
		__m128d XMM1 = _mm_set1_pd(s);
		for(; i<n; i+=2)
		{
			__m128d XMM2 = _mm_load_pd(x+i);
			__m128d XMM3 = _mm_load_pd(y+i);
			__m128d XMM4 = XMM2, XMM5 = XMM3;
			XMM2 = _mm_mul_pd(XMM2, XMM0);
			XMM3 = _mm_mul_pd(XMM3, XMM1);
			XMM2 = _mm_add_pd(XMM2, XMM3);
			XMM4 = _mm_mul_pd(XMM4, XMM1);
			XMM5 = _mm_mul_pd(XMM5, XMM0);
			XMM5 = _mm_sub_pd(XMM5, XMM4);
			_mm_store_pd(x+i, XMM2);
			_mm_store_pd(y+i, XMM5);
		}
#else
		int m = n & ~3;
		for(; i<m; i+=4)
		{
			register MYdouble temp = c*x[i] + s*y[i];
			y[i] = c*y[i] - s*x[i];
			x[i] = temp;
			temp = c*x[i+1] + s*y[i+1];
			y[i+1] = c*y[i+1] - s*x[i+1];
			x[i+1] = temp;
			temp = c*x[i+2] + s*y[i+2];
			y[i+2] = c*y[i+2] - s*x[i+2];
			x[i+2] = temp;
			temp = c*x[i+3] + s*y[i+3];
			y[i+3] = c*y[i+3] - s*x[i+3];
			x[i+3] = temp;
		}
		for(; i<n; ++i)
		{
			register MYdouble temp = c*x[i] + s*y[i];
			y[i] = c*y[i] - s*x[i];
			x[i] = temp;
		}
#endif
	}
	else
	{
		int ix = (incx < 0 ? -(n-1)*incx : 0);
		int iy = (incy < 0 ? -(n-1)*incy : 0);
		for(; i<n; ++i,ix+=incx,iy+=incy)
		{
			register MYdouble temp = c*x[ix] + s*y[iy];
			y[iy] = c*y[iy] - s*x[ix];
			x[ix] = temp;
		}
	}
}

/** \internal
 *  \brief Swap vectors.
 *  \param n number of elements
 *  \param x first vector
 *  \param incx stride for vector x
 *  \param y second vector
 *  \param incy stride for vector y
 *  \ingroup numerics
 */
void dswap(int n, MYdouble *x, int incx, MYdouble *y, int incy)
{
	int i = 0;
	if(n <= 0)
		return;
	else if(incx == 1 && incy == 1)
	{
#if OPENGI_SSE >= 2
		int m = n & ~7;
		for(; i<m; i+=8)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			__m128d XMM1 = _mm_load_pd(x+i+2);
			__m128d XMM2 = _mm_load_pd(x+i+4);
			__m128d XMM3 = _mm_load_pd(x+i+6);
			__m128d XMM4 = _mm_load_pd(y+i);
			__m128d XMM5 = _mm_load_pd(y+i+2);
			__m128d XMM6 = _mm_load_pd(y+i+4);
			__m128d XMM7 = _mm_load_pd(y+i+6);
			_mm_store_pd(x+i, XMM4);
			_mm_store_pd(x+i+2, XMM5);
			_mm_store_pd(x+i+4, XMM6);
			_mm_store_pd(x+i+6, XMM7);
			_mm_store_pd(y+i, XMM0);
			_mm_store_pd(y+i+2, XMM1);
			_mm_store_pd(y+i+4, XMM2);
			_mm_store_pd(y+i+6, XMM3);
		}
		for(; i<n; i+=2)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			__m128d XMM1 = _mm_load_pd(y+i);
			_mm_store_pd(x+i, XMM1);
			_mm_store_pd(y+i, XMM0);
		}
#else
		int m = n & ~3;
		for(; i<m; i+=4)
		{
			register MYdouble temp = x[i];
			x[i] = y[i];
			y[i] = temp;
			temp = x[i+1];
			x[i+1] = y[i+1];
			y[i+1] = temp;
			temp = x[i+2];
			x[i+2] = y[i+2];
			y[i+2] = temp;
			temp = x[i+3];
			x[i+3] = y[i+3];
			y[i+3] = temp;
		}
		for(; i<n; ++i)
		{
			register MYdouble temp = x[i];
			x[i] = y[i];
			y[i] = temp;
		}
#endif
	}
	else
	{
		int ix = (incx < 0 ? -(n-1)*incx : 0);
		int iy = (incy < 0 ? -(n-1)*incy : 0);
		for(; i<n; ++i,ix+=incx,iy+=incy)
		{
			register MYdouble temp = x[ix];
			x[ix] = y[iy];
			y[iy] = temp;
		}
	}
}

/** \internal
 *  \brief Scale vector.
 *  \param n number of elements
 *  \param a scale factor
 *  \param x vector to scale
 *  \param incx stride for vector x
 *  \ingroup numerics
 */
void dscal(int n, MYdouble a, MYdouble *x, int incx)
{
	int i = 0;
	if(n <= 0 || incx <= 0)
		return;
	else if(incx == 1)
	{
		int m = n & ~7;
#if OPENGI_SSE >= 2
		__m128d XMM7 = _mm_set1_pd(a);
		for(; i<m; i+=8)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			__m128d XMM1 = _mm_load_pd(x+i+2);
			__m128d XMM2 = _mm_load_pd(x+i+4);
			__m128d XMM3 = _mm_load_pd(x+i+6);
			XMM0 = _mm_mul_pd(XMM0, XMM7);
			XMM1 = _mm_mul_pd(XMM1, XMM7);
			XMM2 = _mm_mul_pd(XMM2, XMM7);
			XMM3 = _mm_mul_pd(XMM3, XMM7);
			_mm_store_pd(x+i, XMM0);
			_mm_store_pd(x+i+2, XMM1);
			_mm_store_pd(x+i+4, XMM2);
			_mm_store_pd(x+i+6, XMM3);
		}
		for(; i<n; i+=2)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			XMM0 = _mm_mul_pd(XMM0, XMM7);
			_mm_store_pd(x+i, XMM0);
		}
#else
		for(; i<m; i+=8)
		{
			x[i] *= a;
			x[i+1] *= a;
			x[i+2] *= a;
			x[i+3] *= a;
			x[i+4] *= a;
			x[i+5] *= a;
			x[i+6] *= a;
			x[i+7] *= a;
		}
		for(; i<n; ++i)
			x[i] *= a;
#endif
	}
	else
	{
		n *= incx;
		for(; i<n; i+=incx)
			x[i] *= a;
	}
}

/** \internal
 *  \brief Copy vector.
 *  \param n number of elements
 *  \param x source vector
 *  \param incx stride for vector x
 *  \param y destination vector
 *  \param incy stride for vector y
 *  \ingroup numerics
 */
void dcopy(int n, const MYdouble *x, int incx, MYdouble *y, int incy)
{
	int i = 0;
	if(n <= 0)
		return;
	else if(incx == 1 && incy == 1)
	{
#if OPENGI_SSE >= 2
		int m = n & ~15;
		for(; i<m; i+=16)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			__m128d XMM1 = _mm_load_pd(x+i+2);
			__m128d XMM2 = _mm_load_pd(x+i+4);
			__m128d XMM3 = _mm_load_pd(x+i+6);
			__m128d XMM4 = _mm_load_pd(x+i+8);
			__m128d XMM5 = _mm_load_pd(x+i+10);
			__m128d XMM6 = _mm_load_pd(x+i+12);
			__m128d XMM7 = _mm_load_pd(x+i+14);
			_mm_store_pd(y+i, XMM0);
			_mm_store_pd(y+i+2, XMM1);
			_mm_store_pd(y+i+4, XMM2);
			_mm_store_pd(y+i+6, XMM3);
			_mm_store_pd(y+i+8, XMM4);
			_mm_store_pd(y+i+10, XMM5);
			_mm_store_pd(y+i+12, XMM6);
			_mm_store_pd(y+i+14, XMM7);
		}
		for(; i<n; i+=2)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			_mm_store_pd(y+i, XMM0);
		}
#else
		int m = n & ~7;
		for(; i<m; i+=8)
		{
			y[i] = x[i];
			y[i+1] = x[i+1];
			y[i+2] = x[i+2];
			y[i+3] = x[i+3];
			y[i+4] = x[i+4];
			y[i+5] = x[i+5];
			y[i+6] = x[i+6];
			y[i+7] = x[i+7];
		}
		for(; i<n; ++i)
			y[i] = x[i];
#endif
	}
	else
	{
		int ix = (incx < 0 ? -(n-1)*incx : 0);
		int iy = (incy < 0 ? -(n-1)*incy : 0);

		for(; i<n; ++i,ix+=incx,iy+=incy)
			y[iy] = x[ix];
	}
}

/** \internal
 *  \brief Add scaled vector (y = a*x + y).
 *  \param n number of elements
 *  \param a scale factor
 *  \param x first vector
 *  \param incx stride for vector x
 *  \param y second vector
 *  \param incy stride for vector y
 *  \ingroup numerics
 */
void daxpy(int n, MYdouble a, const MYdouble *x, int incx, MYdouble *y, int incy)
{
	int i = 0;
	if(n <= 0 || a == 0.0)
		return;
	else if(incx == 1 && incy == 1)
	{
		int m = n & ~3;
#if OPENGI_SSE >= 2
		__m128d XMM7 = _mm_set1_pd(a);
		for(; i<m; i+=4)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			__m128d XMM1 = _mm_load_pd(x+i+2);
			__m128d XMM2 = _mm_load_pd(y+i);
			__m128d XMM3 = _mm_load_pd(y+i+2);
			XMM0 = _mm_mul_pd(XMM0, XMM7);
			XMM1 = _mm_mul_pd(XMM1, XMM7);
			XMM2 = _mm_add_pd(XMM2, XMM0);
			XMM3 = _mm_add_pd(XMM3, XMM1);
			_mm_store_pd(y+i, XMM2);
			_mm_store_pd(y+i+2, XMM3);
		}
		for(; i<n; i+=2)
		{
			__m128d XMM0 = _mm_load_pd(x+i);
			__m128d XMM1 = _mm_load_pd(y+i);
			XMM0 = _mm_mul_pd(XMM0, XMM7);
			XMM1 = _mm_add_pd(XMM1, XMM0);
			_mm_store_pd(y+i, XMM1);
		}
#else
		for(; i<m; i+=4)
		{
			y[i] += a * x[i];
			y[i+1] += a * x[i+1];
			y[i+2] += a * x[i+2];
			y[i+3] += a * x[i+3];
		}
		for(; i<n; ++i)
			y[i] += a * x[i];
#endif
	}
	else
	{
		int ix = (incx < 0 ? -(n-1)*incx : 0);
		int iy = (incy < 0 ? -(n-1)*incy : 0);
		for(; i<n; ++i,ix+=incx,iy+=incy)
			y[iy] += a * x[ix];
	}
}

/** \internal
 *  \brief Dot product.
 *  \param n number of elements
 *  \param x first vector
 *  \param incx stride for vector x
 *  \param y second vector
 *  \param incy stride for vector y
 *  \return dot product <\a x,\a y>
 *  \ingroup numerics
 */
MYdouble ddot(int n, const MYdouble *x, int incx, const MYdouble *y, int incy)
{
	int i = 0;
	if(n <= 0)
		return 0.0;
	else if(incx == 1 && incy == 1)
	{
		int m = n & ~3;
#if OPENGI_SSE >= 2
		MYdouble temp;
		__m128d XMM0 = _mm_setzero_pd();
#if OPENGI_SSE < 3
		__m128d XMM1;
#endif
		for(; i<m; i+=4)
		{
			__m128d XMM2 = _mm_load_pd(x+i);
			__m128d XMM3 = _mm_load_pd(x+i+2);
			__m128d XMM4 = _mm_load_pd(y+i);
			__m128d XMM5 = _mm_load_pd(y+i+2);
			XMM2 = _mm_mul_pd(XMM2, XMM4);
			XMM3 = _mm_mul_pd(XMM3, XMM5);
			XMM0 = _mm_add_pd(XMM0, XMM2);
			XMM0 = _mm_add_pd(XMM0, XMM3);
		}
		for(; i<n; i+=2)
		{
			__m128d XMM2 = _mm_load_pd(x+i);
			__m128d XMM3 = _mm_load_pd(y+i);
			XMM2 = _mm_mul_pd(XMM2, XMM3);
			XMM0 = _mm_add_pd(XMM0, XMM2);
		}
#if OPENGI_SSE >= 3
		XMM0 = _mm_hadd_pd(XMM0, XMM0);
#else
		XMM1 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(0, 1));
		XMM0 = _mm_add_pd(XMM0, XMM1);
#endif
		_mm_store_sd(&temp, XMM0);
#else
		register MYdouble temp = 0.0;
		for(; i<m; i+=4)
			temp += x[i]*y[i] + x[i+1]*y[i+1] + 
				x[i+2]*y[i+2] + x[i+3]*y[i+3];
		for(; i<n; ++i)
			temp += x[i] * y[i];
#endif
		return temp;
	}
	else
	{
		register MYdouble temp = 0.0;
		int ix = (incx < 0 ? -(n-1)*incx : 0);
		int iy = (incy < 0 ? -(n-1)*incy : 0);
		for(; i<n; ++i,ix+=incx,iy+=incy)
			temp += x[ix] * y[iy];
		return temp;
	}
}

/** \internal
 *  \brief Euclidean norm.
 *  \param n number of elements
 *  \param x vector to compute norm of
 *  \param incx stride for vector x
 *  \return euclidean norm sqrt(<\a x,\a x>)
 *  \ingroup numerics
 */
MYdouble dnrm2(int n, const MYdouble *x, int incx)
{
	int i = 0;
	if(n <= 0)
		return 0.0;
	else if(incx == 1)
	{
		int m = n & ~7;
#if OPENGI_SSE >= 2
		MYdouble temp;
		__m128d XMM0 = _mm_setzero_pd();
#if OPENGI_SSE < 3
		__m128d XMM1;
#endif
		for(; i<m; i+=8)
		{
			__m128d XMM2 = _mm_load_pd(x+i);
			__m128d XMM3 = _mm_load_pd(x+i+2);
			__m128d XMM4 = _mm_load_pd(x+i+4);
			__m128d XMM5 = _mm_load_pd(x+i+6);
			XMM2 = _mm_mul_pd(XMM2, XMM2);
			XMM3 = _mm_mul_pd(XMM3, XMM3);
			XMM4 = _mm_mul_pd(XMM4, XMM4);
			XMM5 = _mm_mul_pd(XMM5, XMM5);
			XMM0 = _mm_add_pd(XMM0, XMM2);
			XMM0 = _mm_add_pd(XMM0, XMM3);
			XMM0 = _mm_add_pd(XMM0, XMM4);
			XMM0 = _mm_add_pd(XMM0, XMM5);
		}
		for(; i<n; i+=2)
		{
			__m128d XMM2 = _mm_load_pd(x+i);
			XMM2 = _mm_mul_pd(XMM2, XMM2);
			XMM0 = _mm_add_pd(XMM0, XMM2);
		}
#if OPENGI_SSE >= 3
		XMM0 = _mm_hadd_pd(XMM0, XMM0);
#else
		XMM1 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(0, 1));
		XMM0 = _mm_add_pd(XMM0, XMM1);
#endif
		XMM0 = _mm_sqrt_sd(XMM0, XMM0);
		_mm_store_sd(&temp, XMM0);
		return temp;
#else
		register MYdouble temp = 0.0;
		for(; i<m; i+=8)
			temp += x[i]*x[i] + x[i+1]*x[i+1] + 
				x[i+2]*x[i+2] + x[i+3]*x[i+3] + 
				x[i+4]*x[i+4] + x[i+5]*x[i+5] + 
				x[i+6]*x[i+6] + x[i+7]*x[i+7];
		for(; i<n; ++i)
			temp += x[i] * x[i];
		return sqrt(temp);
#endif
	}
	else
	{
		register MYdouble temp = 0.0;
		int ix = (incx < 0 ? -(n-1)*incx : 0);
		for(; i<n; ++i,ix+=incx)
			temp += x[ix] * x[ix];
		return temp;
	}
	if(n <= 0 || incx <= 0)
		return 0.0;
	else if(n == 1)
		return fabs(x[0]);
	else
	{
		MYdouble scale = 0.0, ssq = 1.0, absxi, quot;
		int ix = 0;
		if(incx != 1)
			n = (n-1) * incx;
		for(; ix<n; ix+=incx)
		{
			if(x[ix] != 0.0)
			{
				absxi = fabs(x[ix]);
				if(scale < absxi)
				{
					quot = scale / absxi;
					ssq = 1.0 + ssq*quot*quot;
					scale = absxi;
				}
				else
				{
					quot = absxi / scale;
					ssq += quot * quot;
				}
			}
		}
		return scale * sqrt(ssq);
	}
}

/** \internal
 *  \brief Sum absolute values.
 *  \param n number of elements
 *  \param x vector to sum
 *  \param incx stride for vector x
 *  \return absolute sum of vector elements
 *  \ingroup numerics
 */
MYdouble dasum(int n, const MYdouble *x, int incx)
{
	int i = 0;
	if(n <= 0 || incx <= 0)
		return 0.0;
	if(incx == 1)
	{
		int m = n & ~7;
#if OPENGI_SSE >= 2
		MYdouble temp;
		__m128d XMM0 = _mm_setzero_pd();
		__m128d XMM1 = _mm_set1_pd(-0.0);
		for(; i<m; i+=8)
		{
			__m128d XMM2 = _mm_load_pd(x+i);
			__m128d XMM3 = _mm_load_pd(x+i+2);
			__m128d XMM4 = _mm_load_pd(x+i+4);
			__m128d XMM5 = _mm_load_pd(x+i+6);
			XMM2 = _mm_andnot_pd(XMM1, XMM2);
			XMM3 = _mm_andnot_pd(XMM1, XMM3);
			XMM4 = _mm_andnot_pd(XMM1, XMM4);
			XMM5 = _mm_andnot_pd(XMM1, XMM5);
			XMM0 = _mm_add_pd(XMM0, XMM2);
			XMM0 = _mm_add_pd(XMM0, XMM3);
			XMM0 = _mm_add_pd(XMM0, XMM4);
			XMM0 = _mm_add_pd(XMM0, XMM5);
		}
		for(; i<n; i+=2)
		{
			__m128d XMM2 = _mm_load_pd(x+i);
			XMM2 = _mm_andnot_pd(XMM1, XMM2);
			XMM0 = _mm_add_pd(XMM0, XMM2);
		}
#if OPENGI_SSE >= 3
		XMM0 = _mm_hadd_pd(XMM0, XMM0);
#else
		XMM1 = _mm_shuffle_pd(XMM0, XMM0, _MM_SHUFFLE2(0, 1));
		XMM0 = _mm_add_pd(XMM0, XMM1);
#endif
		_mm_store_sd(&temp, XMM0);
#else
		register MYdouble temp = 0.0;
		for(; i<m; i+=8)
			temp += fabs(x[i]) + fabs(x[i+1]) + fabs(x[i+2]) + fabs(x[i+3]) + 
				fabs(x[i+4]) + fabs(x[i+5]) + fabs(x[i+6]) + fabs(x[i+7]);
		for(; i<n; ++i)
			temp += fabs(x[i]);
#endif
		return temp;
	}
	else
	{
		register MYdouble temp = 0.0;
		n *= incx;
		for(; i<n; i+=incx)
			temp += fabs(x[i]);
		return temp;
	}
}

/** Compute largest absolute value.
 *  \param n number of elements
 *  \param x vector to work on
 *  \param incx stride for vector x
 *  \return index of element with largest absolute value
 *  \ingroup numerics
 */
int idamax(int n, const MYdouble *x, int incx)
{
	register MYdouble fmax = fabs(x[0]);
	int i = 1, imax = 0;
	if(n <= 0 || incx <= 0)
		return -1;
	else if(n == 1)
		return 0;
	if(incx == 1)
	{
		for(; i<n; ++i)
		{
			if(fabs(x[i]) > fmax)
			{
				imax = i;
				fmax = fabs(x[i]);
			}
		}
	}
	else
	{
		int ix = incx;
		for(; i<n; ++i,ix+=incx)
		{
			if(fabs(x[ix]) > fmax)
			{
				imax = i;
				fmax = fabs(x[ix]);
			}
		}
	}
	return imax;
}

/** \internal
 *  \brief General matrix vector multiplication (y = alpha*A(')*x + beta*y).
 *  \param trans 'T' for transpose, 'N' else
 *  \param m number of rows of matrix
 *  \param n number of columns of matrix
 *  \param alpha scalar
 *  \param A general m-by-n-matrix in column-major format
 *  \param lda first dimension of matrix
 *  \param x vector to multiply matrix with
 *  \param incx stride for vector x
 *  \param beta scalar
 *  \param y accumulator
 *  \param incy stride for vector y
 *  \ingroup numerics
 */
void dgemv(char trans, int m, int n, MYdouble alpha, const MYdouble *A, int lda, 
		   const MYdouble *x, int incx, MYdouble beta, MYdouble *y, int incy)
{
	int i, ix, iy, j, jx, jy, ij, kx, ky, lenx, leny;
	if(m <= 0 || n <= 0 || lda < m || !incx || !incy || 
	   (alpha == 0.0 && beta == 1.0) || 
	   (trans != 'N' && trans != 'T' && trans != 'C'))
		return;

	if(trans == 'N')
	{
		lenx = n;
		leny = m;
	}
	else
	{
		lenx = m;
		leny = n;
	}
	ky = (incy>0) ? 0 : (-(leny-1)*incy);

	if(beta != 1.0)
	{
		if(incy == 1)
		{
			if(beta == 0.0)
				for(i=0; i<leny; ++i)
					y[i] = 0.0;
			else
				for(i=0; i<leny; ++i)
					y[i] *= beta;
		}
		else
		{
			if(beta == 0.0)
				for(i=0,iy=ky; i<leny; ++i,iy+=incy)
					y[iy] = 0.0;
			else
				for(i=0,iy=ky; i<leny; ++i,iy+=incy)
					y[iy] *= beta;
		}
	}
	if(alpha == 0.0)
		return;

	kx = (incx>0) ? 0 : (-(lenx-1)*incx);
	if(trans == 'N')
	{
		if(incy == 1)
		{
			for(j=0,jx=kx; j<n; ++j,jx+=incx)
			{
				if(x[jx] != 0.0)
				{
					register MYdouble temp = alpha * x[jx];
					for(i=0,ij=lda*j; i<m; ++i,++ij)
						y[i] += temp * A[ij];
				}
			}
		}
		else
		{
			for(j=0,jx=kx; j<n; ++j,jx+=incx)
			{
				if(x[jx] != 0.0)
				{
					register MYdouble temp = alpha * x[jx];
					for(i=0,iy=ky,ij=lda*j; i<m; ++i,iy+=incy,++ij)
						y[iy] += temp * A[ij];
				}
			}
		}
	}
	else
	{
		if(incx == 1)
		{
			for(j=0,jy=ky; j<n; ++j,jy+=incy)
			{
				register MYdouble temp = 0.0;
				for(i=0,ij=lda*j; i<m; ++i,++ij)
					temp += A[ij] * x[i];
				y[j] += alpha * temp;
			}
		}
		else
		{
			for(j=0,jy=ky; j<n; ++j,jy+=incy)
			{
				register MYdouble temp = 0.0;
				for(i=0,ix=kx,ij=lda*j; i<m; ++i,ix+=incx,++ij)
					temp += A[ij] * x[ix];
				y[jy] += alpha * temp;
			}
		}
	}
}

/** \internal
 *  \brief Symmetric matrix vector multiplication (y = alpha*A(')*x + beta*y).
 *  \param uplo 'U' for upper triangle, 'L' for lower
 *  \param n number of rows/columns of matrix
 *  \param alpha scalar
 *  \param A symmetric n-by-n-matrix in column-major format
 *  \param lda first dimension of matrix
 *  \param x vector to multiply matrix with
 *  \param incx stride for vector x
 *  \param beta scalar
 *  \param y accumulator
 *  \param incy stride for vector y
 *  \ingroup numerics
 */
void dsymv(char uplo, int n, MYdouble alpha, const MYdouble *A, int lda, 
		   const MYdouble *x, int incx, MYdouble beta, MYdouble *y, int incy)
{
	int i, ix, iy, j, jx, jy, ij, kx, ky;
	if(n <= 0 || lda < n || !incx || !incy || 
	   (alpha == 0.0 && beta == 1.0) || (uplo != 'U' && uplo != 'L'))
		return;
	ky = (incy>0) ? 0 : (-(n-1)*incy);

	if(beta != 1.0)
	{
		if(incy == 1)
		{
			if(beta == 0.0)
				for(i=0; i<n; ++i)
					y[i] = 0.0;
			else
				for(i=0; i<n; ++i)
					y[i] *= beta;
		}
		else
		{
			if(beta == 0.0)
				for(i=0,iy=ky; i<n; ++i,iy+=incy)
					y[iy] = 0.0;
			else
				for(i=0,iy=ky; i<n; ++i,iy+=incy)
					y[iy] *= beta;
		}
	}
	if(alpha == 0.0)
		return;

	kx = (incx>0) ? 0 : (-(n-1)*incx);
	if(uplo == 'U')
	{
		if(incx == 1 && incy == 1)
		{
			for(j=0; j<n; ++j)
			{
				register MYdouble temp1 = alpha * x[j];
				register MYdouble temp2 = 0.0;
				for(i=0,ij=lda*j; i<j; ++i,++ij)
				{
					y[i] += temp1 * A[ij];
					temp2 += A[ij] * x[i];
				}
				y[j] += temp1*A[ij] + alpha*temp2;
			}
		}
		else
		{
			for(j=0,jx=kx,jy=ky; j<n; ++j,jx+=incx,jy+=incy)
			{
				register MYdouble temp1 = alpha * x[jx];
				register MYdouble temp2 = 0.0;
				for(i=0,ix=kx,iy=ky,ij=lda*j; i<j; ++i,ix+=incx,iy+=incy,++ij)
				{
					y[iy] += temp1 * A[ij];
					temp2 += A[ij] * x[ix];
				}
				y[jy] += temp1*A[ij] + alpha*temp2;
			}
		}
	}
	else
	{
		if(incx == 1 && incy == 1)
		{
			for(j=0; j<n; ++j)
			{
				register MYdouble temp1 = alpha * x[j];
				register MYdouble temp2 = 0.0;
				y[j] += temp1 * A[ij=lda*j+j];
				for(i=j+1,++ij; i<n; ++i,++ij)
				{
					y[i] += temp1 * A[ij];
					temp2 += A[ij] * x[i];
				}
				y[j] += alpha * temp2;
			}
		}
		else
		{
			for(j=0,jx=kx,jy=ky; j<n; ++j,jx+=incx,jy+=incy)
			{
				register MYdouble temp1 = alpha * x[jx];
				register MYdouble temp2 = 0.0;
				y[jy] += temp1 * A[ij=lda*j+j];
				for(i=j+1,ix=kx,iy=ky,++ij; i<n; ++i,ix+=incx,iy+=incy,++ij)
				{
					y[iy] += temp1 * A[ij];
					temp2 += A[ij] * x[ix];
				}
				y[jy] += alpha * temp2;
			}
		}
	}
}

/** \internal
 *  \brief Solve packed triangular system.
 *  \param uplo 'U' for upper triangle, 'L' for lower
 *  \param trans 'T' for transpose, 'N' else
 *  \param diag 'U' for unit diagonal 'N' else
 *  \param n number of unknowns
 *  \param A packed triangular matrix
 *  \param x right hand side and solution
 *  \param incx stride for vector x
 *  \ingroup numerics
 */
void dtpsv(char uplo, char trans, char diag, 
		   int n, const MYdouble *A, MYdouble *x, int incx)
{
	int i, ix, j, jx, k, kk, kx, nounit = (diag == 'N');
	if(n <= 0 || !incx || (uplo != 'U' && uplo != 'L') || 
	   (trans != 'N' && trans != 'T' && trans != 'C') || 
	   (diag != 'U' && !nounit))
		return;
	if(incx <= 0)
		kx = -(n-1) * incx;
	else if(incx != 1)
		kx = 0;
	if(trans == 'N')
	{
		register MYdouble temp;
		if(uplo == 'U')
		{
			kk = (n*(n+1))/2 - 1;
			if(incx == 1)
			{
				for(j=n-1; j>=0; --j)
				{
					if(x[j] != 0.0)
					{
						if(nounit)
							x[j] /= A[kk];
						temp = x[j];
						for(i=j-1,k=kk-1; i>=0; --i,--k)
							x[i] -= A[k] * temp;
					}
					kk -= j + 1;
				}
			}
			else
			{
				for(j=n-1,jx=kx+(n-1)*incx; j>=0; --j,jx-=incx)
				{
					if(x[jx] != 0.0)
					{
						if(nounit)
							x[jx] /= A[kk];
						temp = x[jx];
						for(k=kk-1,ix=jx-incx; k>=kk-j; --k,ix-=incx)
							x[ix] -= A[k] * temp;
					}
					kk -= j + 1;
				}
			}
		}
		else
		{
			kk = 0;
			if(incx == 1)
			{
				for(j=0; j<n; ++j)
				{
					if(x[j] != 0.0)
					{
						if(nounit)
							x[j] /= A[kk];
						temp = x[j];
						for(i=j+1,k=kk+1; i<n; ++i,++k)
							x[i] -= A[k] * temp;
					}
					kk += n - j;
				}
			}
			else
			{
				for(j=0,jx=kx; j<n; ++j,jx+=incx)
				{
					if(x[jx] != 0.0)
					{
						if(nounit)
							x[jx] /= A[kk];
						temp = x[jx];
						for(k=kk+1,ix=jx+incx; k<kk+n-j; ++k,ix+=incx)
							x[ix] -= A[k] * temp;
					}
					kk += n - j;
				}
			}
		}
	}
	else
	{
		if(uplo == 'U')
		{
			kk = 0;
			if(incx == 1)
			{
				for(j=0; j<n; ++j)
				{
					register MYdouble temp = x[j];
					for(i=0,k=kk; i<j; ++i,++k)
						temp -= A[k] * x[i];
					if(nounit)
						temp /= A[kk+j];
					x[j] = temp;
					kk += j + 1;
				}
			}
			else
			{
				for(j=0,jx=kx; j<n; ++j,jx+=incx)
				{
					register MYdouble temp = x[jx];
					for(k=kk,ix=kx; k<kk+j; ++k,ix+=incx)
						temp -= A[k] * x[ix];
					if(nounit)
						temp /= A[kk+j];
					x[jx] = temp;
					kk += j + 1;
				}
			}
		}
		else
		{
			kk = (n*(n+1))/2 - 1;
			if(incx == 1)
			{
				for(j=n-1; j>=0; --j)
				{
					register MYdouble temp = x[j];
					for(i=n-1,k=kk; i>j; --i,--k)
						temp -= A[k] * x[i];
					if(nounit)
						temp /= A[kk-n+j+1];
					x[j] = temp;
					kk -= n - j;
				}
			}
			else
			{
				kx += (n-1) * incx;
				for(j=n-1,jx=kx; j>=0; --j,jx-=incx)
				{
					register MYdouble temp = x[jx];
					for(k=kk,ix=kx; k>kk-n+j+1; --k,ix-=incx)
						temp -= A[k] * x[ix];
					if(nounit)
						temp /= A[kk-n+j+1];
					x[jx] = temp;
					kk -= n - j;
				}
			}
		}
	}
}

/** \internal
 *  \brief General rank 1 operation (A += alpha*x*y').
 *  \param m number of rows of matrix
 *  \param n number of columns of matrix
 *  \param alpha scalar
 *  \param x first vector
 *  \param incx stride for vector x
 *  \param y second vector
 *  \param incy stride for vector y
 *  \param A general m-by-n-matrix in column-major format
 *  \param lda first dimension of matrix
 *  \ingroup numerics
 */
void dger(int m, int n, MYdouble alpha, const MYdouble *x, int incx, 
		  const MYdouble *y, int incy, MYdouble *A, int lda)
{
	int i, ix, j, jy, ij, kx;
	if(m <= 0 || n <= 0 || lda < m || !incx || !incy || alpha == 0.0)
		return;
	if(incx == 1)
	{
		for(j=0,jy = (incy>0) ? 0 : (-(n-1)*incy); j<n; ++j,jy+=incy)
		{
			if(y[jy] != 0.0)
			{
				register MYdouble temp = alpha * y[jy];
				for(i=0,ij=lda*j; i<m; ++i,++ij)
					A[ij] += x[i] * temp;
			}
		}
	}
	else
	{
		kx = (incx>0) ? 0 : (-(m-1)*incx);
		for(j=0,jy = (incy>0) ? 0 : (-(n-1)*incy); j<n; ++j,jy+=incy)
		{
			if(y[jy] != 0.0)
			{
				register MYdouble temp = alpha * y[jy];
				for(i=0,ix=kx,ij=lda*j; i<m; ++i,ix+=incx,++ij)
					A[ij] += x[ix] * temp;
			}
		}
	}
}

/** \internal
 *  \brief Symmetric rank 1 operation (A += alpha*x*x').
 *  \param uplo 'U' for upper triangle, 'L' for lower
 *  \param n number of columns of matrix
 *  \param alpha scalar
 *  \param x first vector
 *  \param incx stride for vector x
 *  \param A symmetric n-by-n-matrix in column-major format
 *  \param lda first dimension of matrix
 *  \ingroup numerics
 */
void dsyr(char uplo, int n, MYdouble alpha, 
		  const MYdouble *x, int incx, MYdouble *A, int lda)
{
	int i, ix, j, jx, ij, kx;
	if(n <= 0 || lda < n || !incx || 
		alpha == 0.0 || (uplo != 'U' && uplo != 'L'))
		return;
	kx = (incx>0) ? 0 : (-(n-1)*incx);
	if(uplo == 'U')
	{
		if(incx == 1)
		{
			for(j=0; j<n; ++j)
			{
				if(x[j] != 0.0)
				{
					register MYdouble temp = alpha * x[j];
					for(i=0,ij=lda*j; i<=j; ++i,++ij)
						A[ij] += x[i] * temp;
				}
			}
		}
		else
		{
			for(j=0,jx=kx; j<n; ++j,jx+=incx)
			{
				if(x[jx] != 0.0)
				{
					register MYdouble temp = alpha * x[jx];
					for(i=0,ix=kx,ij=lda*j; i<=j; ++i,ix+=incx,++ij)
						A[ij] += x[ix] * temp;
				}
			}
		}
	}
	else
	{
		if(incx == 1)
		{
			for(j=0; j<n; ++j)
			{
				if(x[j] != 0.0)
				{
					register MYdouble temp = alpha * x[j];
					for(i=j,ij=lda*j+j; i<n; ++i,++ij)
						A[ij] += x[i] * temp;
				}
			}
		}
		else
		{
			for(j=0,jx=kx; j<n; ++j,jx+=incx)
			{
				if(x[jx] != 0.0)
				{
					register MYdouble temp = alpha * x[jx];
					for(i=j,ix=jx,ij=lda*j+j; i<n; ++i,ix+=incx,++ij)
						A[ij] += x[ix] * temp;
				}
			}
		}
	}
}

/** \internal
 *  \brief Symmetric rank 2 operation (A += alpha*x*y' + alpha*y*x').
 *  \param uplo 'U' for upper triangle, 'L' for lower
 *  \param n number of columns of matrix
 *  \param alpha scalar
 *  \param x first vector
 *  \param incx stride for vector x
 *  \param y second vector
 *  \param incy stride for vector y
 *  \param A symmetric n-by-n-matrix in column-major format
 *  \param lda first dimension of matrix
 *  \ingroup numerics
 */
void dsyr2(char uplo, int n, MYdouble alpha, const MYdouble *x, int incx, 
		   const MYdouble *y, int incy, MYdouble *A, int lda)
{
	int i, ix, iy, j, jx, jy, ij, kx, ky;
	if(n <= 0 || lda < n || !incx || !incy || 
		alpha == 0.0 || (uplo != 'U' && uplo != 'L'))
		return;
	kx = (incx>0) ? 0 : (-(n-1)*incx);
	ky = (incy>0) ? 0 : (-(n-1)*incy);
	if(uplo == 'U')
	{
		if(incx == 1 && incy == 1)
		{
			for(j=0; j<n; ++j)
			{
				if(x[j] != 0.0 || y[j] != 0.0)
				{
					register MYdouble temp1 = alpha * y[j];
					register MYdouble temp2 = alpha * x[j];
					for(i=0,ij=lda*j; i<=j; ++i,++ij)
						A[ij] += x[i]*temp1 + y[i]*temp2;
				}
			}
		}
		else
		{
			for(j=0,jx=kx,jy=ky; j<n; ++j,jx+=incx,jy+=incy)
			{
				if(x[jx] != 0.0 || y[jy] != 0.0)
				{
					register MYdouble temp1 = alpha * y[jy];
					register MYdouble temp2 = alpha * x[jx];
					for(i=0,ix=kx,iy=ky,ij=lda*j; i<=j; ++i,ix+=incx,iy+=incy,++ij)
						A[ij] += x[ix]*temp1 + y[iy]*temp2;
				}
			}
		}
	}
	else
	{
		if(incx == 1 && incy == 1)
		{
			for(j=0; j<n; ++j)
			{
				if(x[j] != 0.0 || y[j] != 0.0)
				{
					register MYdouble temp1 = alpha * y[j];
					register MYdouble temp2 = alpha * x[j];
					for(i=j,ij=lda*j+j; i<n; ++i,++ij)
						A[ij] += x[i]*temp1 + y[i]*temp2;
				}
			}
		}
		else
		{
			for(j=0,jx=kx,jy=ky; j<n; ++j,jx+=incx,jy+=incy)
			{
				if(x[jx] != 0.0 || y[jy] != 0.0)
				{
					register MYdouble temp1 = alpha * y[jy];
					register MYdouble temp2 = alpha * x[jx];
					for(i=j,ix=jx,iy=jy,ij=lda*j+j; i<n; ++i,ix+=incx,iy+=incy,++ij)
						A[ij] += x[ix]*temp1 + y[iy]*temp2;
				}
			}
		}
	}
}

void dgemm(char transa, char transb, int m, int n, int k, 
		   MYdouble alpha, const MYdouble *A, int lda, const MYdouble *B, 
		   int ldb, MYdouble beta, MYdouble *C, int ldc)
{
	int i, j, l, ij, lj, il, jl, li, ij0;
	if(m <= 0 || n <= 0 || k<0 || lda < ((transa=='N') ? m : k) || 
		ldb < ((transb=='N') ? k : n) || ldc < m || 
		((alpha==0.0 || k==0) && beta==1.0) || 
		(transa != 'N' && transa != 'T' && transa != 'C') || 
		(transb != 'N' && transb != 'T' && transb != 'C'))
		return;
	if(alpha == 0.0)
	{
		if(beta == 0.0)
			for(j=0; j<n; ++j)
				for(i=0,ij=ldc*j; i<m; ++i,++ij)
					C[ij] = 0.0;
		else
			for(j=0; j<n; ++j)
				for(i=0,ij=ldc*j; i<m; ++i,++ij)
					C[ij] *= beta;
	}
	else if(transb == 'N')
	{
		if(transa == 'N')
		{
			for(j=0,ij0=0; j<n; ++j,ij0+=ldc)
			{
				if(beta == 0.0)
					for(i=0,ij=ij0; i<m; ++i,++ij)
						C[ij] = 0.0;
				else if(beta != 1.0)
					for(i=0,ij=ij0; i<m; ++i,++ij)
						C[ij] *= beta;
				for(l=0,lj=ldb*j; l<k; ++l,++lj)
				{
					if(B[lj] != 0.0)
					{
						register MYdouble temp = alpha * B[lj];
						for(i=0,ij=ij0,il=lda*l; i<m; ++i,++ij,++il)
							C[ij] += temp * A[il];
					}
				}
			}
		}
		else
		{
			for(j=0; j<n; ++j)
			{
				for(i=0,ij=ldc*j; i<m; ++i,++ij)
				{
					register MYdouble temp = 0.0;
					for(l=0,li=lda*i,lj=ldb*j; l<k; ++l,++li,++lj)
						temp += A[li] * B[lj];
					if(beta == 0.0)
						C[ij] = alpha * temp;
					else
						C[ij] = alpha*temp + beta*C[ij];
				}
			}
		}
	}
	else
	{
		if(transa == 'N')
		{
			for(j=0,ij0=0; j<n; ++j,ij0+=ldc)
			{
				if(beta == 0.0)
					for(i=0,ij=ij0; i<m; ++i,++ij)
						C[ij] = 0.0;
				else if(beta != 1.0)
					for(i=0,ij=ij0; i<m; ++i,++ij)
						C[ij] *= beta;
				for(l=0,jl=j; l<k; ++l,jl+=ldb)
				{
					if(B[jl] != 0.0)
					{
						register MYdouble temp = alpha * B[jl];
						for(i=0,ij=ij0,il=lda*l; i<m; ++i,++ij,++il)
							C[ij] += temp * A[il];
					}
				}
			}
		}
		else
		{
			for(j=0; j<n; ++j)
			{
				for(i=0,ij=ldc*j; i<m; ++i,++ij)
				{
					register MYdouble temp = 0.0;
					for(l=0,li=lda*i,jl=j; l<k; ++l,++li,jl+=ldb)
						temp += A[li] * B[jl];
					if(beta == 0.0)
						C[ij] = alpha * temp;
					else
						C[ij] = alpha*temp + beta*C[ij];
				}
			}
		}
	}
}
