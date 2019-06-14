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
 *  \brief Declaration of common BLAS functions.
 */

#ifndef __GI_BLAS_H__
#define __GI_BLAS_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif

#include "my_types.h"
/*************************************************************************/
/* Functions */

/** \name Level 1 BLAS routines
 *  \{
 */
void drotg(MYdouble *a, MYdouble *b, MYdouble *c, MYdouble *s);
void drot(int n, MYdouble *x, int incx, MYdouble *y, int incy, MYdouble c, MYdouble s);
void dswap(int n, MYdouble *x, int incx, MYdouble *y, int incy);
void dscal(int n, MYdouble a, MYdouble *x, int incx);
void dcopy(int n, const MYdouble *x, int incx, MYdouble *y, int incy);
void daxpy(int n, MYdouble a, const MYdouble *x, int incx, MYdouble *y, int incy);
MYdouble ddot(int n, const MYdouble *x, int incx, const MYdouble *y, int incy);
MYdouble dnrm2(int n, const MYdouble *x, int incx);
MYdouble dasum(int n, const MYdouble *x, int incx);
int idamax(int n, const MYdouble *x, int incx);
/** \} */

/** \name Level 2 BLAS routines
 *  \{
 */
void dgemv(char trans, int m, int n, MYdouble alpha, const MYdouble *A, int lda, 
	const MYdouble *x, int incx, MYdouble beta, MYdouble *y, int incy);
void dsymv(char uplo, int n, MYdouble alpha, const MYdouble *A, int lda, 
	const MYdouble *x, int incx, MYdouble beta, MYdouble *y, int incy);
void dger(int m, int n, MYdouble alpha, const MYdouble *x, int incx, 
	const MYdouble *y, int incy, MYdouble *A, int lda);
void dsyr(char uplo, int n, MYdouble alpha, 
	const MYdouble *x, int incx, MYdouble *A, int lda);
void dsyr2(char uplo, int n, MYdouble alpha, const MYdouble *x, 
	int incx, const MYdouble *y, int incy, MYdouble *A, int lda);
void dtpsv(char uplo, char trans, char diag, int n, 
	const MYdouble *A, MYdouble *x, int incx);
/** \} */

/** \name Level 3 BLAS routines
 *  \{
 */
void dgemm(char transa, char transb, int m, int n, int k, 
	MYdouble alpha, const MYdouble *A, int lda, const MYdouble *B, 
	int ldb, MYdouble beta, MYdouble *C, int ldc);
/** \} */


#endif
