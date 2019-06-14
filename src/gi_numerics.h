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
 *  \brief Declaration of structures and functions for numerical linear algebra.
 */

#ifndef __GI_NUMERICS_H__
#define __GI_NUMERICS_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include <stdio.h>


/*************************************************************************/
/* Typedefs */

/* forward declaration */
//typedef struct _GISparseMatrix GISparseMatrix;

/** \internal
 *  \brief Matrix-vector multiplication function
 *  \ingroup numerics
 */
typedef void (*GImvfunc)(const struct _GISparseMatrix*, const GIdouble*, GIdouble*);

/** \internal
 *  \brief Matrix manipulation function
 *  \ingroup numerics
 */
typedef void (*GImfunc)(struct _GISparseMatrix*);

/** \internal
 *  \brief Iterative solving function
 *  \ingroup numerics
 */
typedef GIuint (*GIsolverfunc)(const struct _GISparseMatrix*, const GIdouble*, GIdouble*, GImvfunc, GImvfunc, GIdouble, GIuint);


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Element of sparse vector.
 *  \details This structure represents an element of a dynamic sparse vector.
 *  \ingroup numerics
 */
typedef struct _GIVectorElement
{
	GIuint		index;							/**< Index of item in vector. */
	GIdouble	value;							/**< Value of item. */
} GIVectorElement;

/** \internal
 *  \brief Sparse vector.
 *  \details This structure represents a dynamic sparse vector.
 *  \ingroup numerics
 */
typedef struct _GISparseVector
{
	GIuint			size;						/**< Number of non-zero elements. */
	GIuint			capacity;					/**< Size of vector. */
	GIVectorElement	*elements;					/**< Vector elements. */
} GISparseVector;

/** \internal
 *  \brief Square matrix.
 *  \details This structure represents the base class for square matrices.
 *  \ingroup numerics
 */
typedef struct _GISparseMatrix
{
	GIuint		n;								/**< Size of matrix. */
	GIboolean	symmetric;						/**< Symmetric matrix. */
	GIvoid		*data;							/**< Custom data (used by preconditioners). */
} GISparseMatrix;

/** \internal
 *  \brief Dynamic sparse matrix.
 *  \details This structure represents a dynamic sparse matrix in LIL format.
 *  \ingroup numerics
 */
typedef struct _GISparseMatrixLIL
{
	GIuint			n;							/**< Numbr of rows/columns. */
	GIboolean		symmetric;					/**< Symmetric matrix. */
	GIvoid			*data;						/**< Custom data (used by preconditioners) */
	GISparseVector	*rows;						/**< Rows of matrix. */
} GISparseMatrixLIL;

/** \internal
 *  \brief Compressed sparse matrix.
 *  \details This structure represents a sparse matrix in CSR format.
 *  \ingroup numerics
 */
typedef struct _GISparseMatrixCSR
{
	GIuint		n;								/**< Numbr of rows/columns. */
	GIboolean	symmetric;						/**< Symmetric matrix. */
	GIvoid		*data;							/**< Custom data (used by preconditioners). */
	GIuint		nnz;							/**< Number of non-zero elements. */
	GIdouble	*values;						/**< Non-zero elements. */
	GIuint		*idx;							/**< Column indices. */
	GIuint		*ptr;							/**< Start indices of rows. */
} GISparseMatrixCSR;

/** \internal
 *  \brief Block compressed sparse matrix.
 *  \details This structure represents a sparse matrix in 2x2-BCSR format.
 *  \ingroup numerics
 */
typedef struct _GISparseMatrixBCSR2
{
	GIuint		n;								/**< Numbr of rows/columns. */
	GIboolean	symmetric;						/**< Symmetric matrix. */
	GIvoid		*data;							/**< Custom data (used by preconditioners). */
	GIuint		nnz;							/**< Number of non-zero blocks. */
	GIdouble	*values;						/**< Non-zero elements. */
	GIuint		*idx;							/**< Column indices. */
	GIuint		*ptr;							/**< Start indices of rows. */
} GISparseMatrixBCSR2;


/*************************************************************************/
/* Functions */

/** \name Sparse vector methods
 *  \{
 */
void GISparseVector_construct(GISparseVector *vec);
void GISparseVector_grow(GISparseVector *vec);
void GISparseVector_set(GISparseVector *vec, GIuint index, GIdouble value);
void GISparseVector_add(GISparseVector *vec, GIuint index, GIdouble value);
void GISparseVector_append(GISparseVector *vec, GIuint index, GIdouble value);
void GISparseVector_to_zero(GISparseVector *vec);
void GISparseVector_clear(GISparseVector *vec);
/** \} */

/** \name Sparse matrix methods
 *  \{
 */
void GISparseMatrixLIL_construct(GISparseMatrixLIL *mat, GIuint n, GIboolean symmetric);
void GISparseMatrixLIL_destruct(GISparseMatrixLIL *mat);
void GISparseMatrixLIL_set(GISparseMatrixLIL *mat, GIuint i, GIuint j, GIdouble value);
void GISparseMatrixLIL_add(GISparseMatrixLIL *mat, GIuint i, GIuint j, GIdouble value);
void GISparseMatrixLIL_append(GISparseMatrixLIL *mat, GIuint i, GIuint j, GIdouble value);
void GISparseMatrixLIL_to_zero(GISparseMatrixLIL *mat);
void GISparseMatrixLIL_clear(GISparseMatrixLIL *mat);
GIuint GISparseMatrixLIL_nnz(const GISparseMatrixLIL *mat);
void GISparseMatrixLIL_print(const GISparseMatrixLIL *mat, FILE *file);
void GISparseMatrixLIL_ax(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
/** \} */

/** \name Compressed matrix methods
 *  \{
 */
void GISparseMatrixCSR_construct(GISparseMatrixCSR *mat, const GISparseMatrixLIL *src);
void GISparseMatrixCSR_destruct(GISparseMatrixCSR *mat);
void GISparseMatrixCSR_print(const GISparseMatrixCSR *mat, FILE *file);
void GISparseMatrixCSR_ax(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
/** \} */

/** \name Block compressed matrix methods
 *  \{
 */
void GISparseMatrixBCSR2_construct(GISparseMatrixBCSR2 *mat, const GISparseMatrixLIL *src);
void GISparseMatrixBCSR2_destruct(GISparseMatrixBCSR2 *mat);
void GISparseMatrixBCSR2_print(const GISparseMatrixBCSR2 *mat, FILE *file);
void GISparseMatrixBCSR2_ax(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
/** \} */

/** \name Preconditioner preparation
 *  \{
 */
void GISparseMatrixLIL_prepare_jacobi(GISparseMatrixLIL *mat);
void GISparseMatrixLIL_prepare_ssor(GISparseMatrixLIL *mat, GIdouble omega);
void GISparseMatrixLIL_prepare_ilu(GISparseMatrixLIL *mat);
void GISparseMatrixCSR_prepare_jacobi(GISparseMatrixCSR *mat);
void GISparseMatrixCSR_prepare_ssor(GISparseMatrixCSR *mat, GIdouble omega);
void GISparseMatrixCSR_prepare_ilu(GISparseMatrixCSR *mat);
void GISparseMatrixBCSR2_prepare_jacobi(GISparseMatrixBCSR2 *mat);
/** \} */

/** \name Preconditioners
 *  \{
 */
void GISparseMatrix_pc_jacobi(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
void GISparseMatrixLIL_pc_ssor(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
void GISparseMatrixLIL_pc_ilu(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
void GISparseMatrixCSR_pc_ssor(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
void GISparseMatrixCSR_pc_ilu(const GISparseMatrix *A, const GIdouble *x, GIdouble *y);
#define GISparseMatrixLIL_pc_jacobi		GISparseMatrix_pc_jacobi
#define GISparseMatrixCSR_pc_jacobi		GISparseMatrix_pc_jacobi
#define GISparseMatrixBCSR2_pc_jacobi	GISparseMatrix_pc_jacobi
/** \} */

/** \name Solvers
 *  \{
 */
GIuint GISolver_cg(const GISparseMatrix *A, const GIdouble *b, GIdouble *x, 
	GImvfunc ax, GImvfunc pc, GIdouble eps, GIuint max_iter);
GIuint GISolver_bicgstab(const GISparseMatrix *A, const GIdouble *b, GIdouble *x, 
	GImvfunc ax, GImvfunc pc, GIdouble eps, GIuint max_iter);
GIuint GISolver_gmres(const GISparseMatrix *A, const GIdouble *b, GIdouble *x, 
	GImvfunc ax, GImvfunc pc, GIdouble eps, GIuint max_iter);
/** \} */

/** \name Matrix methods
 *  \{
 */
GIboolean GImat3d_tridiagonalize(GIdouble *mat, GIdouble *diag, GIdouble *subdiag);
GIuint GImat3d_ql(GIdouble *mat, GIdouble *diag, GIdouble *subdiag, GIuint max_iter);
/** \} */


#endif
