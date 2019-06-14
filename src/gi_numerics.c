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
 *  \brief Implementation of structures and functions for numerical linear algebra.
 */

#include "gi_numerics.h"
#include "gi_blas.h"
#include "gi_memory.h"
#include "gi_math.h"

#include <math.h>
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


/** \internal
 *  \brief Apply IC/ILU preconditioner with compressed matrix.
 *  \param A system matrix
 *  \param values value array to use instead of matrix values
 *  \param x vector to multiply preconditioning matrix with
 *  \param y vector to store result
 *  \ingroup numerics
 */
static void incomplete_lu(const GISparseMatrixCSR *A, 
						  const GIdouble *values, 
						  const GIdouble *x, GIdouble *y)
{
	GIint i, ij, N = A->n;
	if(!values)
		values = A->values;

	/* forward-eliminate for lower triangle */
	for(i=0; i<N; ++i)
	{
		register GIdouble temp = 0.0;
		for(ij=A->ptr[i]; A->idx[ij]<i; ++ij)
			temp += values[ij] * y[A->idx[ij]];
		y[i] = (x[i]-temp) * values[ij];
	}

	/* backward-eliminate for upper triangle */
	if(A->symmetric)
	{
		for(i=N-1,ij=A->nnz-1; i>0; --i)
		{
			register GIdouble temp = y[i] *= values[ij--];
			for(; ij>=A->ptr[i]; --ij)
				y[A->idx[ij]] -= values[ij] * temp;
		}
		y[0] *= values[ij];
	}
	else
	{
		for(i=N-2; i>=0; --i)
		{
			register GIdouble temp = 0.0;
			for(ij=A->ptr[i+1]-1; A->idx[ij]>i; --ij)
				temp += values[ij] * y[A->idx[ij]];
			y[i] -= temp;
		}
	}
}

/** \internal
 *  \brief Sparse vector constructor.
 *  \param vec vector to construct
 *  \ingroup numerics
 */
void GISparseVector_construct(GISparseVector *vec)
{
	/* empty vector */
	vec->size = vec->capacity = 0;
	vec->elements = NULL;
}

/** \internal
 *  \brief Let sparse vector grow.
 *  \param vec vector to grow
 *  \ingroup numerics
 */
void GISparseVector_grow(GISparseVector *vec)
{
	if(vec->capacity)
	{
		/* double size */
		vec->capacity <<= 1;
		vec->elements = (GIVectorElement*)GI_REALLOC_ARRAY(vec->elements, 
			vec->capacity, sizeof(GIVectorElement));
	}
	else
	{
		/* init elements */
		vec->capacity = 4;
		vec->elements = (GIVectorElement*)GI_MALLOC_ARRAY(
			vec->capacity, sizeof(GIVectorElement));
	}
}

/** \internal
 *  \brief Set value of sparse vector element.
 *  \param vec vector to work on
 *  \param index index of element
 *  \param value new value of element
 *  \ingroup numerics
 */
void GISparseVector_set(GISparseVector *vec, GIuint index, GIdouble value)
{
	GIuint i;

	/* element allready in vector */
	for(i=0; i<vec->size && vec->elements[i].index<index; ++i) ;
	if(i < vec->size && vec->elements[i].index == index)
		vec->elements[i].value = value;
	else
		GISparseVector_append(vec, index, value);
}

/** \internal
 *  \brief Add value to sparse vector element.
 *  \param vec vector to work on
 *  \param index index of element
 *  \param value value to add
 *  \ingroup numerics
 */
void GISparseVector_add(GISparseVector *vec, GIuint index, GIdouble value)
{
	GIuint i;

	for(i=0; i<vec->size && vec->elements[i].index<index; ++i) ;
	if(i < vec->size && vec->elements[i].index == index)
		vec->elements[i].value += value;
	else
		GISparseVector_append(vec, index, value);
}

/** \internal
 *  \brief Append element to sparse vector.
 *  \param vec vector to work on
 *  \param index index of element
 *  \param value value of element
 *  \ingroup numerics
 */
void GISparseVector_append(GISparseVector *vec, GIuint index, GIdouble value)
{
	GIuint i;

	/* just add new element (grow if neccessary) */
	if(vec->size == vec->capacity)
		GISparseVector_grow(vec);
	for(i=vec->size; i>0 && index<vec->elements[i-1].index; --i)
		vec->elements[i] = vec->elements[i-1];
	vec->elements[i].index = index;
	vec->elements[i].value = value;
	++vec->size;
}

/** \internal
 *  \brief Set sparse vector to zero vector.
 *  \param vec vector to zero
 *  \ingroup numerics
 */
void GISparseVector_to_zero(GISparseVector *vec)
{
	/* just no elements */
	vec->size = 0;
}

/** \internal
 *  \brief Clear sparse vector.
 *  \param vec vector to clear
 *  \ingroup numerics
 */
void GISparseVector_clear(GISparseVector *vec)
{
	/* clear data */
	vec->size = vec->capacity = 0;
	if(vec->elements)
	{
		GI_FREE_ARRAY(vec->elements);
		vec->elements = NULL;
	}
}

/** \internal
 *  \brief Sparse matrix constructor.
 *  \param mat matrix to construct
 *  \param n number of rows/columns
 *  \param symmetric GI_TRUE if matrix is symmetric GI_FALSE else
 *  \ingroup numerics
 */
void GISparseMatrixLIL_construct(GISparseMatrixLIL *mat, 
								 GIuint n, GIboolean symmetric)
{
	/* init data */
	mat->n = n;
	mat->symmetric = symmetric;
	mat->data = NULL;

	/* create rows */
	mat->rows = (GISparseVector*)GI_CALLOC_ARRAY(n, sizeof(GISparseVector));
}

/** \internal
 *  \brief Sparse matrix destructor.
 *  \param mat matrix to destruct
 *  \ingroup numerics
 */
void GISparseMatrixLIL_destruct(GISparseMatrixLIL *mat)
{
	GIuint i;

	/* clear and delete rows */
	for(i=0; i<mat->n; ++i)
		GISparseVector_clear(mat->rows+i);
	GI_FREE_ARRAY(mat->rows);

	/* clear data */
	if(mat->data)
		GI_FREE_ARRAY(mat->data);
	memset(mat, 0, sizeof(GISparseMatrixLIL));
}

/** \internal
 *  \brief Set value of sparse matrix element.
 *  \param mat matrix to work on
 *  \param i column of element
 *  \param j row of element
 *  \param value new value of element
 *  \ingroup numerics
 */
void GISparseMatrixLIL_set(GISparseMatrixLIL *mat, GIuint i, 
						   GIuint j, GIdouble value)
{
	/* set row element */
	if(mat->symmetric && (j > i))
	{
		register GIuint temp;
		GI_SWAP(i, j, temp);
	}
	GISparseVector_set(mat->rows+i, j, value);
}

/** \internal
 *  \brief Add value to sparse matrix element.
 *  \param mat matrix to work on
 *  \param i column of element
 *  \param j row of element
 *  \param value value to add
 *  \ingroup numerics
 */
void GISparseMatrixLIL_add(GISparseMatrixLIL *mat, GIuint i, 
						   GIuint j, GIdouble value)
{
	/* add to row element */
	if(mat->symmetric && (j > i))
	{
		register GIuint temp;
		GI_SWAP(i, j, temp);
	}
	GISparseVector_add(mat->rows+i, j, value);
}

/** \internal
 *  \brief Append value to sparse matrix.
 *  \param mat matrix to work on
 *  \param i column of element
 *  \param j row of element
 *  \param value value to add
 *  \ingroup numerics
 */
void GISparseMatrixLIL_append(GISparseMatrixLIL *mat, GIuint i, GIuint j, GIdouble value)
{
	/* append element to row */
	if(mat->symmetric)
	{
		if(j > i)
		{
			register GIuint temp;
			GI_SWAP(i, j, temp);
		}
		GISparseVector_set(mat->rows+i, j, value);
	}
	else
		GISparseVector_append(mat->rows+i, j, value);
}

/** \internal
 *  \brief Set sparse matrix to zero matrix.
 *  \param mat matrix to clear
 *  \ingroup numerics
 */
void GISparseMatrixLIL_to_zero(GISparseMatrixLIL *mat)
{
	GIuint i, N = mat->n;

	/* zero row */
	for(i=0; i<N; ++i)
		GISparseVector_to_zero(mat->rows+i);
}

/** \internal
 *  \brief Clear sparse matrix (set to zero matrix and release memory).
 *  \param mat matrix to clear
 *  \ingroup numerics
 */
void GISparseMatrixLIL_clear(GISparseMatrixLIL *mat)
{
	GIuint i, N = mat->n;

	/* clear row */
	for(i=0; i<N; ++i)
		GISparseVector_clear(mat->rows+i);
}

/** \internal
 *  \brief Get number of non-zero elements in sparse matrix.
 *  \param mat matrix to query
 *  \return number of non-zero elements
 *  \ingroup numerics
 */
GIuint GISparseMatrixLIL_nnz(const GISparseMatrixLIL *mat)
{
	GIuint i, N = mat->n;
	register GIuint uiResult = 0;

	/* count row sizes */
	for(i=0; i<N; ++i)
		uiResult += mat->rows[i].size;
	return uiResult;
}

/** \internal
 *  \brief print sparse matrix to file.
 *  \param mat matrix to print
 *  \param file file to print to
 *  \ingroup numerics
 */
void GISparseMatrixLIL_print(const GISparseMatrixLIL *mat, FILE *file)
{
	GISparseVector *pVec;
	GIVectorElement *pElement;
	GIuint i, j, ij, N = mat->n;

	if(mat->symmetric)
	{
		/* convert to dense matrix and print */
		GIdouble *out = (GIdouble*)GI_CALLOC_ARRAY(N*N, sizeof(GIdouble));
		for(i=0; i<N; ++i)
		{
			pVec = mat->rows + i;
			for(ij=0; ij<pVec->size; ++ij)
			{
				pElement = pVec->elements + ij;
				out[i*N+pElement->index] = pElement->value;
				out[pElement->index*N+i] = pElement->value;
			}
		}
		for(i=0,ij=0; i<N; ++i)
		{
			for(j=0; j<N; ++j,++ij)
				fprintf(file, "%.17f ", out[ij]);
			fprintf(file, "\n");
		}
		GI_FREE_ARRAY(out);
	}
	else
	{
		/* print rows and fill gaps with zeros */
		for(i=0,j=0; i<N; ++i,j=0)
		{
			pVec = mat->rows + i;
			for(ij=0; ij<pVec->size; ++ij,++j)
			{
				pElement = pVec->elements + ij;
				for(; j<pElement->index; ++j)
					fprintf(file, "0 ");
				fprintf(file, "%.17f ", pElement->value);
			}
			for(; j<N; ++j)
				fprintf(file, "0 ");
			fprintf(file, "\n");
		}
	}
}

/** \internal
 *  \brief Multiply sparse matrix by vector
 *  \param A matrix
 *  \param x vector to multiply with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrixLIL_ax(const GISparseMatrix *A, 
						  const GIdouble *x, GIdouble *y)
{
	const GISparseMatrixLIL *mat = (const GISparseMatrixLIL*)A;
	const GISparseVector *pVec;
	const GIVectorElement *pElement;
	GIuint i, j, N = mat->n;

	if(mat->symmetric)
	{
		/* symmetric mutliplication */
		for(i=0; i<N; ++i)
		{
			register GIdouble temp = 0.0, xi = x[i];
			for(pElement=mat->rows[i].elements; pElement->index<i; ++pElement)
			{
				temp += pElement->value * x[pElement->index];
				y[pElement->index] += pElement->value * xi;
			}
			y[i] = temp + pElement->value*x[pElement->index];
		}
	}
	else
	{
		/* general mutliplication */
		for(i=0; i<N; ++i)
		{
			register GIdouble temp = 0.0;
			pVec = mat->rows + i;
			for(j=0,pElement=pVec->elements; j<pVec->size; ++j,++pElement)
				temp += pElement->value * x[pElement->index];
			y[i] = temp;
		}
	}
}

/** \internal
 *  \brief Compressed matrix constructor.
 *  \param mat matrix to construct
 *  \param src dynamic sparse matrix to construct from
 *  \ingroup numerics
 */
void GISparseMatrixCSR_construct(GISparseMatrixCSR *mat, 
								 const GISparseMatrixLIL *src)
{
	const GISparseVector *pVec;
	GIint i, j, c, N = src->n;

	/* create matrix */
	mat->n = N;
	mat->symmetric = src->symmetric;
	mat->data = NULL;
	mat->nnz = GISparseMatrixLIL_nnz(src);
	mat->values = (GIdouble*)GI_MALLOC_ARRAY(mat->nnz, sizeof(GIdouble));
	mat->idx = (GIuint*)GI_MALLOC_ARRAY(mat->nnz, sizeof(GIuint));
	mat->ptr = (GIuint*)GI_MALLOC_ARRAY(N+1, sizeof(GIuint));

	/* copy data */
	for(i=0,c=0; i<N; ++i)
	{
		pVec = src->rows + i;
		mat->ptr[i] = c;
		for(j=0; j<pVec->size; ++j,++c)
		{
			mat->values[c] = pVec->elements[j].value;
			mat->idx[c] = pVec->elements[j].index;
		}
	}
	mat->ptr[N] = c;
}

/** \internal
 *  \brief Compressed matrix destructor.
 *  \param mat matrix to destruct
 *  \ingroup numerics
 */
void GISparseMatrixCSR_destruct(GISparseMatrixCSR *mat)
{
	/* delete arrays and clear data */
	GI_FREE_ARRAY(mat->values);
	GI_FREE_ARRAY(mat->idx);
	GI_FREE_ARRAY(mat->ptr);
	if(mat->data)
		GI_FREE_ARRAY(mat->data);
	memset(mat, 0, sizeof(GISparseMatrixCSR));
}

/** \internal
 *  \brief print compressed matrix to file.
 *  \param mat matrix to print
 *  \param file file to print to
 *  \ingroup numerics
 */
void GISparseMatrixCSR_print(const GISparseMatrixCSR *mat, FILE *file)
{
	GIuint i, j, ij, N = mat->n;

	if(mat->symmetric)
	{
		/* convert to dense matrix and print */
		GIdouble *out = (GIdouble*)GI_CALLOC_ARRAY(N*N, sizeof(GIdouble));
		for(i=0,ij=0; i<N; ++i)
		{
			for(; ij<mat->ptr[i+1]; ++ij)
			{
				out[i*N+mat->idx[ij]] = mat->values[ij];
				out[mat->idx[ij]*N+i] = mat->values[ij];
			}
		}
		for(i=0,ij=0; i<N; ++i)
		{
			for(j=0; j<N; ++j,++ij)
				fprintf(file, "%.17f ", out[ij]);
			fprintf(file, "\n");
		}
		GI_FREE_ARRAY(out);
	}
	else
	{
		/* print rows and fill gaps with zeros */
		for(i=0,j=0,ij=0; i<N; ++i,j=0)
		{
			for(; ij<mat->ptr[i+1]; ++ij,++j)
			{
				for(; j<mat->idx[ij]; ++j)
					fprintf(file, "0 ");
				fprintf(file, "%.17f ", mat->values[ij]);
			}
			for(; j<N; ++j)
				fprintf(file, "0 ");
			fprintf(file, "\n");
		}
	}
}

/** \internal
 *  \brief Multiply compressed matrix by vector
 *  \param A matrix
 *  \param x vector to multiply with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrixCSR_ax(const GISparseMatrix *A, 
						  const GIdouble *x, GIdouble *y)
{
	const GISparseMatrixCSR *mat = (const GISparseMatrixCSR*)A;
	GIuint i, ij, N = mat->n;

	if(mat->symmetric)
	{
		/* symmetric mutliplication */
		for(i=0,ij=0; i<N; ++i,++ij)
		{
			register GIdouble temp = 0.0, xi = x[i];
			for(; ij<mat->ptr[i+1]-1; ++ij)
			{
				temp += mat->values[ij] * x[mat->idx[ij]];
				y[mat->idx[ij]] += mat->values[ij] * xi;
			}
			y[i] = temp + mat->values[ij]*x[mat->idx[ij]];
		}
	}
	else
	{
		/* general mutliplication */
		for(i=0,ij=0; i<N; ++i)
		{
			register GIdouble temp = 0.0;
			for(; ij<mat->ptr[i+1]; ++ij)
				temp += mat->values[ij] * x[mat->idx[ij]];
			y[i] = temp;
		}
	}
}

/** \internal
 *  \brief Block compressed matrix constructor.
 *  \param mat matrix to construct
 *  \param src dynamic sparse matrix to construct from
 *  \ingroup numerics
 */
void GISparseMatrixBCSR2_construct(GISparseMatrixBCSR2 *mat, 
								   const GISparseMatrixLIL *src)
{
	GIint i, i2, c = -1, N = (src->n+1) >> 1, NNZ = GISparseMatrixLIL_nnz(src);
	GIdouble *pValues;
	GIuint *pIdx;

	/* create matrix */
	mat->n = src->n;
	mat->symmetric = src->symmetric;
	mat->data = NULL;
	mat->nnz = 0;
	mat->ptr = (GIuint*)GI_MALLOC_ARRAY(N+1, sizeof(GIuint));
	pValues = (GIdouble*)GI_MALLOC_ARRAY(NNZ<<2, sizeof(GIdouble));
	pIdx = (GIuint*)GI_MALLOC_ARRAY(NNZ, sizeof(GIuint));

	/* copy data */
	for(i=0,i2=0; i<N; ++i,i2+=2)
	{
		const GISparseVector *pVec[2] = { src->rows+i2, 
			(i2+1<src->n) ? src->rows+i2+1 : NULL };
		GIuint ij[2] = { 0, 0 };
		GIint jCur = -1, j, r;
		mat->ptr[i] = c + 1;
		while(pVec[0] || pVec[1])
		{
			if(pVec[0] && pVec[1])
				r = (pVec[0]->elements[ij[0]].index<
					pVec[1]->elements[ij[1]].index) ? 0 : 1;
			else
				r = pVec[0] ? 0 : 1;
			j = pVec[r]->elements[ij[r]].index >> 1;
			if(j > jCur)
			{
				pIdx[++c] = jCur = j;
				pValues[c<<2] = pValues[(c<<2)+1] = 
					pValues[(c<<2)+2] = pValues[(c<<2)+3] = 0.0;
			}
			pValues[(c<<2)+(r<<1)+(pVec[r]->elements[ij[r]].index&1)] = 
				pVec[r]->elements[ij[r]].value;
			if(++ij[r] == pVec[r]->size)
				pVec[r] = NULL;
		}
	}
	mat->ptr[N] = mat->nnz = ++c;

	/* complete blocks for diagonal of symmetric matrix */
	if(mat->symmetric)
	{
		for(i=0; i<N; ++i)
		{
			GIuint ij = (mat->ptr[i+1]-1) << 2;
			pValues[ij+1] = pValues[ij+2];
		}
	}

	/* compress data */
	mat->values = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE((mat->nnz<<2)*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	memcpy(mat->values, pValues, (mat->nnz<<2)*sizeof(GIdouble));
	GI_FREE_ARRAY(pValues);
	mat->idx = (GIuint*)GI_MALLOC_ARRAY(mat->nnz, sizeof(GIuint));
	memcpy(mat->idx, pIdx, mat->nnz*sizeof(GIuint));
	GI_FREE_ARRAY(pIdx);
}

/** \internal
 *  \brief Block compressed matrix destructor.
 *  \param mat matrix to destruct
 *  \ingroup numerics
 */
void GISparseMatrixBCSR2_destruct(GISparseMatrixBCSR2 *mat)
{
	/* delete arrays and clear data */
	GI_FREE_ALIGNED(mat->values);
	GI_FREE_ARRAY(mat->idx);
	GI_FREE_ARRAY(mat->ptr);
	if(mat->data)
		GI_FREE_ARRAY(mat->data);
	memset(mat, 0, sizeof(GISparseMatrixBCSR2));
}

/** \internal
 *  \brief print block compressed matrix to file.
 *  \param mat matrix to print
 *  \param file file to print to
 *  \ingroup numerics
 */
void GISparseMatrixBCSR2_print(const GISparseMatrixBCSR2 *mat, FILE *file)
{
	GIuint i, j, ij, ij4, N = (mat->n+1) >> 1, M = N << 1;

	if(mat->symmetric)
	{
		/* convert to dense matrix and print */
		GIdouble *out = (GIdouble*)GI_CALLOC_ARRAY(M*M, sizeof(GIdouble));
		for(i=0,ij=0,ij4=0; i<N; ++i)
		{
			for(; ij<mat->ptr[i+1]; ++ij,ij4+=4)
			{
				GIuint I = i << 1, J = mat->idx[ij] << 1;
				out[I*M+J] = out[J*M+I] = mat->values[ij4];
				out[(I+1)*M+J] = out[J*M+I+1] = mat->values[ij4+1];
				out[I*M+J+1] = out[(J+1)*M+I] = mat->values[ij4+2];
				out[(I+1)*M+J+1] = out[(J+1)*M+I+1] = mat->values[ij4+3];
			}
		}
		for(i=0,ij=0; i<M; ++i)
		{
			for(j=0; j<M; ++j,++ij)
				fprintf(file, "%.17f ", out[ij]);
			fprintf(file, "\n");
		}
		GI_FREE_ARRAY(out);
	}
	else
	{
		/* print rows and fill gaps with zeros */
		for(i=0,j=0; i<N; ++i)
		{
			for(ij=mat->ptr[i],j=0; ij<mat->ptr[i+1]; ++ij,++j)
			{
				for(; j<mat->idx[ij]; ++j)
					fprintf(file, "0 0 ");
				fprintf(file, "%.17f %.17f ", 
					mat->values[ij<<2], mat->values[(ij<<2)+1]);
			}
			for(; j<N; ++j)
				fprintf(file, "0 0 ");
			fprintf(file, "\n");
			for(ij=mat->ptr[i],j=0; ij<mat->ptr[i+1]; ++ij,++j)
			{
				for(; j<mat->idx[ij]; ++j)
					fprintf(file, "0 0 ");
				fprintf(file, "%.17f %.17f ", 
					mat->values[(ij<<2)+2], mat->values[(ij<<2)+3]);
			}
			for(; j<N; ++j)
				fprintf(file, "0 0 ");
			fprintf(file, "\n");
		}
	}
}

/** \internal
 *  \brief Multiply block compressed matrix by vector
 *  \param A matrix
 *  \param x vector to multiply with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrixBCSR2_ax(const GISparseMatrix *A, const GIdouble *x, GIdouble *y)
{
	const GISparseMatrixCSR *mat = (const GISparseMatrixCSR*)A;
	GIuint i, j2, ij, ij4, N = (mat->n+1) >> 1;

	if(mat->symmetric)
	{
		/* symmetric mutliplication */
		for(i=0,ij=0,ij4=0; i<N; ++i,++ij,ij4+=4)
		{
#if OPENGI_SSE >= 2
			GIuint i2 = i << 1;
			__m128d XMM0 = _mm_setzero_pd(), XMM1;
			__m128d XMM2 = _mm_load1_pd(x+i2);
			__m128d XMM3 = _mm_load1_pd(x+i2+1);
			__m128d XMM4, XMM5, XMM6, XMM7;
			for(; ij<mat->ptr[i+1]-1; ++ij,ij4+=4)
			{
				j2 = mat->idx[ij] << 1;
				XMM1 = _mm_load_pd(y+j2);
				XMM4 = _mm_load_pd(mat->values+ij4);
				XMM5 = _mm_load_pd(mat->values+ij4+2);
				XMM6 = _mm_mul_pd(XMM4, XMM2);
				XMM7 = _mm_mul_pd(XMM5, XMM3);
				XMM1 = _mm_add_pd(XMM1, XMM6);
				XMM1 = _mm_add_pd(XMM1, XMM7);
				_mm_store_pd(y+j2, XMM1);
				XMM1 = _mm_load_pd(x+j2);
				XMM4 = _mm_mul_pd(XMM4, XMM1);
				XMM5 = _mm_mul_pd(XMM5, XMM1);
#if OPENGI_SSE >= 3
				XMM4 = _mm_hadd_pd(XMM4, XMM5);
				XMM0 = _mm_add_pd(XMM0, XMM4);
#else
				XMM6 = _mm_shuffle_pd(XMM4, XMM5, _MM_SHUFFLE2(0, 0));
				XMM7 = _mm_shuffle_pd(XMM4, XMM5, _MM_SHUFFLE2(1, 1));
				XMM0 = _mm_add_pd(XMM0, XMM6);
				XMM0 = _mm_add_pd(XMM0, XMM7);
#endif	/* SSE3 */
			}
			XMM1 = _mm_load_pd(x+(mat->idx[ij]<<1));
			XMM2 = _mm_load_pd(mat->values+ij4);
			XMM3 = _mm_load_pd(mat->values+ij4+2);
			XMM2 = _mm_mul_pd(XMM2, XMM1);
			XMM3 = _mm_mul_pd(XMM3, XMM1);
#if OPENGI_SSE >= 3
			XMM2 = _mm_hadd_pd(XMM2, XMM3);
			XMM0 = _mm_add_pd(XMM0, XMM2);
#else
			XMM4 = _mm_shuffle_pd(XMM2, XMM3, _MM_SHUFFLE2(0, 0));
			XMM5 = _mm_shuffle_pd(XMM2, XMM3, _MM_SHUFFLE2(1, 1));
			XMM0 = _mm_add_pd(XMM0, XMM4);
			XMM0 = _mm_add_pd(XMM0, XMM5);
#endif	/* SSE3 */
			_mm_store_pd(y+i2, XMM0);
#else
			register GIdouble temp0 = 0.0, temp1 = 0.0;
			GIuint i2 = i << 1;
			for(; ij<mat->ptr[i+1]-1; ++ij,ij4+=4)
			{
				j2 = mat->idx[ij] << 1;
				temp0 += mat->values[ij4]*x[j2] + mat->values[ij4+1]*x[j2+1];
				temp1 += mat->values[ij4+2]*x[j2] + mat->values[ij4+3]*x[j2+1];
				y[j2] += mat->values[ij4]*x[i2] + mat->values[ij4+2]*x[i2+1];
				y[j2+1] += mat->values[ij4+1]*x[i2] + mat->values[ij4+3]*x[i2+1];
			}
			j2 = mat->idx[ij] << 1;
			y[i2] = temp0 + mat->values[ij4]*x[j2] + mat->values[ij4+1]*x[j2+1];
			y[i2+1] = temp1 + mat->values[ij4+2]*x[j2] + mat->values[ij4+3]*x[j2+1];
#endif /* SSE2 */
		}
	}
	else
	{
		/* general mutliplication */
		for(i=0,ij=0,ij4=0; i<N; ++i,y+=2)
		{
#if OPENGI_SSE >= 2
			__m128d XMM0 = _mm_setzero_pd();
			__m128d XMM1, XMM2, XMM3;
#if OPENGI_SSE < 3
			__m128d XMM4, XMM5;
#endif
			for(; ij<mat->ptr[i+1]; ++ij,ij4+=4)
			{
				XMM1 = _mm_load_pd(x+(mat->idx[ij]<<1));
				XMM2 = _mm_load_pd(mat->values+ij4);
				XMM3 = _mm_load_pd(mat->values+ij4+2);
				XMM2 = _mm_mul_pd(XMM2, XMM1);
				XMM3 = _mm_mul_pd(XMM3, XMM1);
#if OPENGI_SSE >= 3
				XMM2 = _mm_hadd_pd(XMM2, XMM3);
				XMM0 = _mm_add_pd(XMM0, XMM2);
#else
				XMM4 = _mm_shuffle_pd(XMM2, XMM3, _MM_SHUFFLE2(0, 0));
				XMM5 = _mm_shuffle_pd(XMM2, XMM3, _MM_SHUFFLE2(1, 1));
				XMM0 = _mm_add_pd(XMM0, XMM4);
				XMM0 = _mm_add_pd(XMM0, XMM5);
#endif	/* SSE3 */
			}
			_mm_store_pd(y, XMM0);
#else
			register GIdouble temp0 = 0.0, temp1 = 0.0;
			for(; ij<mat->ptr[i+1]; ++ij,ij4+=4)
			{
				j2 = mat->idx[ij] << 1;
				temp0 += mat->values[ij4]*x[j2] + mat->values[ij4+1]*x[j2+1];
				temp1 += mat->values[ij4+2]*x[j2] + mat->values[ij4+3]*x[j2+1];
			}
			y[0] = temp0;
			y[1] = temp1;
#endif	/* SSE2 */
		}
	}
}

/** \internal
 *  \brief Prepare data for Jacobi preconditioner with sparse matrix.
 *  \param mat matrix to create data for
 *  \ingroup numerics
 */
void GISparseMatrixLIL_prepare_jacobi(GISparseMatrixLIL *mat)
{
	GIdouble *pData;
	GIuint i, N = mat->n, uiSize = N * sizeof(GIdouble);

	/* create data if neccessary */
	if(!mat->data || *((GIuint*)mat->data) != uiSize)
	{
		if(mat->data)
			GI_FREE_ARRAY(mat->data);
		mat->data = GI_MALLOC_ARRAY(uiSize+sizeof(GIuint), 1);
		*((GIuint*)mat->data) = uiSize;
	}
	pData = (GIdouble*)((GIuint*)mat->data+1);

	/* store inverted diagonal */
	if(mat->symmetric)
		for(i=0; i<N; ++i)
			pData[i] = 1.0 / mat->rows[i].elements[mat->rows[i].size-1].value;
	else
	{
		GIVectorElement *pElement;

		/* search for diagonal entry */
		for(i=0; i<N; ++i)
		{
			for(pElement=mat->rows[i].elements; pElement->index<i; ++pElement) ;
			pData[i] = 1.0 / pElement->value;
		}
	}
}

/** \internal
 *  \brief Prepare data for SSOR preconditioner with sparse matrix.
 *  \param mat matrix to create data for
 *  \param omega relaxation parameter in [0,2]
 *  \ingroup numerics
 */
void GISparseMatrixLIL_prepare_ssor(GISparseMatrixLIL *mat, GIdouble omega)
{
	GIdouble *pData;
	GIuint i, N = mat->n, uiSize = (N+1) * sizeof(GIdouble);

	/* create data if neccessary */
	if(!mat->data || *((GIuint*)mat->data) != uiSize)
	{
		if(mat->data)
			GI_FREE_ARRAY(mat->data);
		mat->data = GI_MALLOC_ARRAY(uiSize+sizeof(GIuint), 1);
		*((GIuint*)mat->data) = uiSize;
	}
	pData = (GIdouble*)((GIuint*)mat->data+1);

	/* save diagonal and omega */
	if(mat->symmetric)
		for(i=0; i<N; ++i)
			pData[i] = mat->rows[i].elements[mat->rows[i].size-1].value;
	else
	{
		GIVectorElement *pElement;

		/* search for diagonal entry */
		for(i=0; i<N; ++i)
		{
			for(pElement=mat->rows[i].elements; pElement->index<i; ++pElement) ;
			pData[i] = pElement->value;
		}
	}
	pData[N] = omega;
}

/** \internal
 *  \brief Prepare data for IC/ILU preconditioner with sparse matrix.
 *  \param mat matrix to create data for
 *  \ingroup numerics
 */
void GISparseMatrixLIL_prepare_ilu(GISparseMatrixLIL *mat)
{
	GISparseMatrixCSR *pData;
	GIVectorElement *pElement;
	GIuint uiSize = sizeof(GISparseMatrixCSR);
	GIuint e, i, j, k, ij, ik, jk, kj, N = mat->n;

	/* create data if neccessary */
	if(!mat->data || *((GIuint*)mat->data) != uiSize)
	{
		if(mat->data)
			GI_FREE_ARRAY(mat->data);
		mat->data = GI_MALLOC_ARRAY(uiSize+sizeof(GIuint), 1);
		*((GIuint*)mat->data) = uiSize;
	}
	pData = (GISparseMatrixCSR*)((GIuint*)mat->data+1);

	/* create compressed matrix */
	pData->n = N;
	pData->nnz = GISparseMatrixLIL_nnz(mat);
	pData->symmetric = mat->symmetric;
	pData->data = NULL;
	pData->values = (GIdouble*)GI_MALLOC_ARRAY(pData->nnz, sizeof(GIdouble));
	pData->idx = (GIuint*)GI_MALLOC_ARRAY(pData->nnz, sizeof(GIuint));
	pData->ptr = (GIuint*)GI_MALLOC_ARRAY(N+1, sizeof(GIuint));

	if(mat->symmetric)
	{
		/* compute incomplete Cholesky factorization */
		for(i=0,ij=0; i<N; ++i,++ij)
		{
			register GIdouble temp = 0.0;
			pData->ptr[i] = ij;
			for(pElement=mat->rows[i].elements; pElement->index<i; ++pElement,++ij)
			{
				register GIdouble temp2 = 0.0;
				pData->idx[ij] = j = pElement->index;
				for(ik=pData->ptr[i],jk=pData->ptr[j]; ik<ij; ++ik)
				{
					for(; pData->idx[jk]<pData->idx[ik]; ++jk) ;
					if(pData->idx[ik] == pData->idx[jk])
						temp2 += pData->values[ik] * pData->values[jk];
				}
				pData->values[ij] = (pElement->value-temp2) * 
					pData->values[pData->ptr[j+1]-1];
				temp += pData->values[ij] * pData->values[ij];
			}
			pData->values[ij] = 1.0 / sqrt(pElement->value-temp);
			pData->idx[ij] = i;
		}
		pData->ptr[N] = pData->nnz;
	}
	else
	{
		/* compute incomplete LU factorization */
		for(i=0,ij=0; i<N; ++i)
		{
			pData->ptr[i] = ij;
			for(e=0,pElement=mat->rows[i].elements; 
				e<mat->rows[i].size; ++e,++pElement,++ij)
			{
				register GIdouble temp = 0.0;
				pData->idx[ij] = j = pElement->index;
				if(i >= j)
				{
					for(ik=pData->ptr[i]; ik<ij; ++ik)
					{
						k = pData->idx[ik];
						for(kj=pData->ptr[k]; kj<pData->ptr[k+1]; ++kj)
						{
							if(pData->idx[kj] == j)
							{
								temp += pData->values[ik] * pData->values[kj];
								break;
							}
						}
					}
					pData->values[ij] = pElement->value - temp;
					if(i == j)
						pData->values[ij] = 1.0 / pData->values[ij];
				}
				else
				{
					for(ik=pData->ptr[i]; pData->idx[ik]<i; ++ik)
					{
						k = pData->idx[ik];
						for(kj=pData->ptr[k]; kj<pData->ptr[k+1]; ++kj)
						{
							if(pData->idx[kj] == j)
							{
								temp += pData->values[ik] * pData->values[kj];
								break;
							}
						}
					}
					pData->values[ij] = (pElement->value-temp) * pData->values[ik];
				}
			}
		}
		pData->ptr[N] = pData->nnz;
	}
}

/** \internal
 *  \brief Prepare data for Jacobi preconditioner with compressed matrix.
 *  \param mat matrix to create data for
 *  \ingroup numerics
 */
void GISparseMatrixCSR_prepare_jacobi(GISparseMatrixCSR *mat)
{
	GIdouble *pData;
	GIuint i, N = mat->n, uiSize = N * sizeof(GIdouble);

	/* create data if neccessary */
	if(!mat->data || *((GIuint*)mat->data) != uiSize)
	{
		if(mat->data)
			GI_FREE_ARRAY(mat->data);
		mat->data = GI_MALLOC_ARRAY(uiSize+sizeof(GIuint), 1);
		*((GIuint*)mat->data) = uiSize;
	}
	pData = (GIdouble*)((GIuint*)mat->data+1);

	/* save inverted diagonal */
	if(mat->symmetric)
	{
		for(i=0; i<N; ++i)
			pData[i] = 1.0 / mat->values[mat->ptr[i+1]-1];
	}
	else
	{
		GIuint ij;

		for(i=0; i<N; ++i)
		{
			for(ij=mat->ptr[i]; mat->idx[ij]<i; ++ij) ;
			pData[i] = 1.0 / mat->values[ij];
		}
	}
}

/** \internal
 *  \brief Prepare data for SSOR preconditioner with compressed matrix.
 *  \param mat matrix to create data for
 *  \param omega relaxation parameter in [0,2]
 *  \ingroup numerics
 */
void GISparseMatrixCSR_prepare_ssor(GISparseMatrixCSR *mat, GIdouble omega)
{
	GIdouble *pData;
	GIuint i, N = mat->n, uiSize = (N+1) * sizeof(GIdouble);

	/* create data if neccessary */
	if(!mat->data || *((GIuint*)mat->data) != uiSize)
	{
		if(mat->data)
			GI_FREE_ARRAY(mat->data);
		mat->data = GI_MALLOC_ARRAY(uiSize+sizeof(GIuint), 1);
		*((GIuint*)mat->data) = uiSize;
	}
	pData = (GIdouble*)((GIuint*)mat->data+1);

	/* save diagonal and omega */
	if(mat->symmetric)
	{
		for(i=0; i<N; ++i)
			pData[i] = mat->values[mat->ptr[i+1]-1];
	}
	else
	{
		GIuint ij;

		for(i=0; i<N; ++i)
		{
			for(ij=mat->ptr[i]; mat->idx[ij]<i; ++ij) ;
			pData[i] = mat->values[ij];
		}
	}
	pData[N] = omega;
}

/** \internal
 *  \brief Prepare data for IC/ILU preconditioner with compressed matrix.
 *  \param mat matrix to create data for
 *  \ingroup numerics
 */
void GISparseMatrixCSR_prepare_ilu(GISparseMatrixCSR *mat)
{
	GIdouble *pData;
	GIuint uiSize = mat->nnz * sizeof(GIdouble);
	GIuint i, j, k, ij, ik, jk, kj, N = mat->n;

	/* create data if neccessary */
	if(!mat->data || *((GIuint*)mat->data) != uiSize)
	{
		if(mat->data)
			GI_FREE_ARRAY(mat->data);
		mat->data = GI_MALLOC_ARRAY(uiSize+sizeof(GIuint), 1);
		*((GIuint*)mat->data) = uiSize;
	}
	pData = (GIdouble*)((GIuint*)mat->data+1);

	if(mat->symmetric)
	{
		/* compute incomplete Cholesky factorization */
		for(i=0,ij=0; i<N; ++i,++ij)
		{
			register GIdouble temp = 0.0;
			for(; mat->idx[ij]<i; ++ij)
			{
				register GIdouble temp2 = 0.0;
				j = mat->idx[ij];
				for(ik=mat->ptr[i],jk=mat->ptr[j]; ik<ij; ++ik)
				{
					for(; mat->idx[jk]<mat->idx[ik]; ++jk) ;
					if(mat->idx[ik] == mat->idx[jk])
						temp2 += pData[ik] * pData[jk];
				}
				pData[ij] = (mat->values[ij]-temp2) * pData[mat->ptr[j+1]-1];
				temp += pData[ij] * pData[ij];
			}
			pData[ij] = 1.0 / sqrt(mat->values[ij]-temp);
		}
	}
	else
	{
		/* compute incomplete LU factorization */
		for(i=0,ij=0; i<N; ++i)
		{
			for(; ij<mat->ptr[i+1]; ++ij)
			{
				register GIdouble temp = 0.0;
				j = mat->idx[ij];
				if(i >= j)
				{
					for(ik=mat->ptr[i]; ik<ij; ++ik)
					{
						k = mat->idx[ik];
						for(kj=mat->ptr[k]; kj<mat->ptr[k+1]; ++kj)
						{
							if(mat->idx[kj] == j)
							{
								temp += pData[ik] * pData[kj];
								break;
							}
						}
					}
					pData[ij] = mat->values[ij] - temp;
					if(i == j)
						pData[ij] = 1.0 / pData[ij];
				}
				else
				{
					for(ik=mat->ptr[i]; mat->idx[ik]<i; ++ik)
					{
						k = mat->idx[ik];
						for(kj=mat->ptr[k]; kj<mat->ptr[k+1]; ++kj)
						{
							if(mat->idx[kj] == j)
							{
								temp += pData[ik] * pData[kj];
								break;
							}
						}
					}
					pData[ij] = (mat->values[ij]-temp) * pData[ik];
				}
			}
		}
	}
}

/** \internal
 *  \brief Prepare data for Jacobi preconditioner with block compressed matrix.
 *  \param mat matrix to create data for
 *  \ingroup numerics
 */
void GISparseMatrixBCSR2_prepare_jacobi(GISparseMatrixBCSR2 *mat)
{
	GIdouble *pData;
	GIuint i, N = (mat->n+1) >> 1, uiSize = (N<<1) * sizeof(GIdouble);

	/* create data if neccessary */
	if(!mat->data || *((GIuint*)mat->data) != uiSize)
	{
		if(mat->data)
			GI_FREE_ARRAY(mat->data);
		mat->data = GI_MALLOC_ARRAY(uiSize+sizeof(GIuint), 1);
		*((GIuint*)mat->data) = uiSize;
	}
	pData = (GIdouble*)((GIuint*)mat->data+1);

	/* save inverted diagonal */
	if(mat->symmetric)
	{
		for(i=0; i<N; ++i,pData+=2)
		{
			pData[0] = 1.0 / mat->values[(mat->ptr[i+1]-1)<<2];
			pData[1] = 1.0 / mat->values[((mat->ptr[i+1]-1)<<2)+3];
		}
	}
	else
	{
		GIuint ij;

		for(i=0; i<N; ++i,pData+=2)
		{
			for(ij=mat->ptr[i]; mat->idx[ij]<i; ++ij) ;
			pData[0] = 1.0 / mat->values[ij<<2];
			pData[1] = 1.0 / mat->values[(ij<<2)+3];
		}
	}
}

/** \internal
 *  \brief Apply Jacobi preconditioner with generic matrix.
 *  \param A system matrix
 *  \param x vector to multiply preconditioning matrix with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrix_pc_jacobi(const GISparseMatrix *A, const GIdouble *x, GIdouble *y)
{
	GIuint i, N = A->n;
	const GIdouble *pInvDiag = (const GIdouble*)((const GIuint*)A->data+1);

	/* divide vector by diagonal of matrix */
	for(i=0; i<N; ++i)
		y[i] = x[i] * pInvDiag[i];
}

/** \internal
 *  \brief Apply SSOR preconditioner with sparse matrix.
 *  \param A system matrix
 *  \param x vector to multiply preconditioning matrix with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrixLIL_pc_ssor(const GISparseMatrix *A, const GIdouble *x, GIdouble *y)
{
	const GISparseMatrixLIL *mat = (const GISparseMatrixLIL*)A;
	const GIVectorElement *pElement;
	const GIdouble *pDiag = (const GIdouble*)((const GIuint*)A->data+1);
	GIint i, N = mat->n;
	GIdouble dOmega = pDiag[N];

	/* forward-eliminate for lower triangle */
	for(i=0; i<N; ++i)
	{
		register GIdouble temp = 0.0;
		for(pElement=mat->rows[i].elements; pElement->index<i; ++pElement)
			temp += pElement->value * y[pElement->index];
		y[i] = (x[i]-dOmega*temp) / pDiag[i];
	}

	/* multiply by diagonal */
	for(i=0; i<N; ++i)
		y[i] *= pDiag[i];

	/* backward-eliminate for upper triangle */
	if(mat->symmetric)
	{
		for(i=N-1; i>0; --i)
		{
			register GIdouble temp = dOmega * (y[i]/=pDiag[i]);
			for(pElement=mat->rows[i].elements; pElement->index<i; ++pElement)
				y[pElement->index] -= pElement->value * temp;
		}
		y[0] /= pDiag[0];
	}
	else
	{
		for(i=N-1; i>=0; --i)
		{
			register GIdouble temp = 0.0;
			for(pElement=mat->rows[i].elements+mat->rows[i].size-1; 
				pElement->index>i; --pElement)
				temp += pElement->value * y[pElement->index];
			y[i] = (y[i]-dOmega*temp) / pDiag[i];
		}
	}
}

/** \internal
 *  \brief Apply IC/ILU preconditioner with sparse matrix.
 *  \param A system matrix
 *  \param x vector to multiply preconditioning matrix with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrixLIL_pc_ilu(const GISparseMatrix *A, const GIdouble *x, GIdouble *y)
{
	incomplete_lu((const GISparseMatrixCSR*)((GIuint*)A->data+1), NULL, x, y);
}

/** \internal
 *  \brief Apply SSOR preconditioner with compressed matrix.
 *  \param A system matrix
 *  \param x vector to multiply preconditioning matrix with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrixCSR_pc_ssor(const GISparseMatrix *A, const GIdouble *x, GIdouble *y)
{
	const GISparseMatrixCSR *mat = (const GISparseMatrixCSR*)A;
	const GIdouble *pDiag = (const GIdouble*)((const GIuint*)A->data+1);
	GIint i, ij, N = mat->n;
	GIdouble dOmega = pDiag[N];

	/* forward-eliminate for lower triangle */
	for(i=0,ij=0; i<N; ++i)
	{
		register GIdouble temp = 0.0;
		for(ij=mat->ptr[i]; mat->idx[ij]<i; ++ij)
			temp += mat->values[ij] * y[mat->idx[ij]];
		y[i] = (x[i]-dOmega*temp) / mat->values[ij];
	}

	/* multiply by diagonal */
	for(i=0; i<N; ++i)
		y[i] *= pDiag[i];

	/* backward-eliminate for upper triangle */
	if(mat->symmetric)
	{
		for(i=N-1,ij=mat->nnz-2; i>0; --i,--ij)
		{
			register GIdouble temp = dOmega * (y[i]/=pDiag[i]);
			for(; ij>=mat->ptr[i]; --ij)
				y[mat->idx[ij]] -= mat->values[ij] * temp;
		}
		y[0] /= mat->values[0];
	}
	else
	{
		for(i=N-1; i>=0; --i)
		{
			register GIdouble temp = 0.0;
			for(ij=mat->ptr[i+1]-1; mat->idx[ij]>i; --ij)
				temp += mat->values[ij] * y[mat->idx[ij]];
			y[i] = (y[i]-dOmega*temp) / mat->values[ij];
		}
	}
}

/** \internal
 *  \brief Apply IC/ILU preconditioner with compressed matrix.
 *  \param A system matrix
 *  \param x vector to multiply preconditioning matrix with
 *  \param y vector to store result
 *  \ingroup numerics
 */
void GISparseMatrixCSR_pc_ilu(const GISparseMatrix *A, const GIdouble *x, GIdouble *y)
{
	incomplete_lu((const GISparseMatrixCSR*)A, 
		(const GIdouble*)((const GIuint*)A->data+1), x, y);
}

/** \internal
 *  \brief Solve equation system by conjugate gradient method.
 *  \param A system matrix
 *  \param b right hand side vector
 *  \param x vector of unknowns
 *  \param ax matrix-vector-multiplication function
 *  \param pc preconditioning function or NULL if no preconditioning
 *  \param eps error threshold
 *  \param max_iter maximum number of iterations
 *  \return number of used iterations
 *  \ingroup numerics
 */
GIuint GISolver_cg(const GISparseMatrix *A, const GIdouble *b, GIdouble *x, 
				   GImvfunc ax, GImvfunc pc, GIdouble eps, GIuint max_iter)
{
	GIuint N = A->n;
	GIdouble *r = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *q = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *v = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *w = (pc ? v : r);
	GIdouble alpha, beta, gamma, tol = eps * eps * ddot(N, b, 1, b, 1);
	GIuint i = 0;

	/* initialize */
	ax(A, x, r);
	daxpy(N, -1.0, b, 1, r, 1);
	if(pc)
		pc(A, r, w);
	dcopy(N, w, 1, q, 1);
	gamma = ddot(N, r, 1, w, 1);

	/* iterate */
	for(i=1; i<=max_iter; ++i)
	{
		ax(A, q, v);
		alpha = -gamma / ddot(N, v, 1, q, 1);
		daxpy(N, alpha, q, 1, x, 1);
		daxpy(N, alpha, v, 1, r, 1);
		if(pc)
			pc(A, r, w);
		beta = 1.0 / gamma;
		if(pc)
		{
			if(ddot(N, r, 1, r, 1) <= tol)
				break;
			gamma = ddot(N, r, 1, w, 1);
		}
		else
		{
			gamma = ddot(N, r, 1, r, 1);
			if(gamma <= tol)
				break;
		}
		beta *= gamma;
		dscal(N, beta, q, 1);
		daxpy(N, 1.0, w, 1, q, 1);
	}

	/* clean up */
	GI_FREE_ALIGNED(r);
	GI_FREE_ALIGNED(q);
	GI_FREE_ALIGNED(v);
	return i;
}

/** \internal
 *  \brief Solve equation system by stabilized biconjugate gradient method.
 *  \param A system matrix
 *  \param b right hand side vector
 *  \param x vector of unknowns
 *  \param ax matrix-vector-multiplication function
 *  \param pc preconditioning function or NULL if no preconditioning
 *  \param eps error threshold
 *  \param max_iter maximum number of iterations
 *  \return number of used iterations
 *  \ingroup numerics
 */
GIuint GISolver_bicgstab(const GISparseMatrix *A, const GIdouble *b, GIdouble *x, 
						 GImvfunc ax, GImvfunc pc, GIdouble eps, GIuint max_iter)
{
	GIuint N = A->n;
	GIdouble *r = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *r0 = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *q = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *v = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *t = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *s = r, *tP = t, *rP = r, *sP = r, *vP = v;
	GIdouble alpha, beta, gamma, omega, tol = eps * eps * ddot(N, b, 1, b, 1);
	GIuint i = 0;

	/* initialize */
	ax(A, x, r);
	daxpy(N, -1.0, b, 1, r, 1);
	if(pc)
	{
		vP = (GIdouble*)GI_MALLOC_ALIGNED(
			GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
		sP = rP = (GIdouble*)GI_MALLOC_ALIGNED(
			GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
		tP = v;
		pc(A, r, rP);
	}
	dcopy(N, rP, 1, q, 1);
	dcopy(N, rP, 1, r0, 1);
	gamma = ddot(N, rP, 1, r0, 1);

	/* iterate */
	for(i=1; i<=max_iter; ++i)
	{
		ax(A, q, v);
		if(pc)
			pc(A, v, vP);
		alpha = gamma / ddot(N, vP, 1, r0, 1);
		daxpy(N, -alpha, q, 1, x, 1);
		daxpy(N, -alpha, v, 1, s, 1);
		if(ddot(N, r, 1, r, 1) <= tol)
			break;
		if(pc)
			daxpy(N, -alpha, vP, 1, sP, 1);
		ax(A, sP, t);
		if(pc)
			pc(A, t, tP);
		omega = ddot(N, tP, 1, sP, 1) / ddot(N, tP, 1, tP, 1);
		daxpy(N, -omega, sP, 1, x, 1);
		daxpy(N, -omega, t, 1, r, 1);
		if(ddot(N, r, 1, r, 1) <= tol)
			break;
		if(pc)
			daxpy(N, -omega, tP, 1, rP, 1);
		beta = alpha / (omega*gamma);
		gamma = ddot(N, rP, 1, r0, 1);
		beta *= gamma;
		daxpy(N, -omega, vP, 1, q, 1);
		dscal(N, beta, q, 1);
		daxpy(N, 1.0, rP, 1, q, 1);
	}

	/* clean up */
	GI_FREE_ALIGNED(r);
	GI_FREE_ALIGNED(r0);
	GI_FREE_ALIGNED(q);
	GI_FREE_ALIGNED(v);
	GI_FREE_ALIGNED(t);
	if(pc)
	{
		GI_FREE_ALIGNED(vP);
		GI_FREE_ALIGNED(rP);
	}
	return i;
}

/** \internal
 *  \brief Solve equation system by generalized minimized residual method.
 *  \param A system matrix
 *  \param b right hand side vector
 *  \param x vector of unknowns
 *  \param ax matrix-vector-multiplication function
 *  \param pc preconditioning function or NULL if no preconditioning
 *  \param eps error threshold
 *  \param max_iter maximum number of iterations
 *  \return number of used iterations
 *  \ingroup numerics
 */
GIuint GISolver_gmres(const GISparseMatrix *A, const GIdouble *b, GIdouble *x, 
					  GImvfunc ax, GImvfunc pc, GIdouble eps, GIuint max_iter)
{
	GIuint N = A->n, M = max_iter & 0xFF;
#if OPENGI_SSE >= 2
	GIuint LDQ = (N+1) & (~1);
#else
	GIuint LDQ = N;
#endif
	GIdouble *Q = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(LDQ*(M+1)*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *H = (GIdouble*)GI_MALLOC_ALIGNED(
		((M*(M+1))/2)*sizeof(GIdouble), sizeof(GIdouble));
	GIdouble *r = (GIdouble*)GI_MALLOC_ALIGNED(
		GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
	GIdouble *c = (GIdouble*)GI_MALLOC_ALIGNED(M*sizeof(GIdouble), sizeof(GIdouble));
	GIdouble *s = (GIdouble*)GI_MALLOC_ALIGNED(M*sizeof(GIdouble), sizeof(GIdouble));
	GIdouble *y = (GIdouble*)GI_MALLOC_ALIGNED((M+1)*sizeof(GIdouble), sizeof(GIdouble));
	GIdouble *w = r;
	GIdouble beta, tol, tmp, hjj;
	GIuint i = 0, j, k, Hij;

	/* outer initialization */
	max_iter >>= 8;
	if(pc)
	{
		w = (GIdouble*)GI_MALLOC_ALIGNED(
			GI_SSE_SIZE(N*sizeof(GIdouble)), GI_SSE_ALIGN_DOUBLE);
		pc(A, b, w);
		tol = eps * eps * ddot(N, w, 1, w, 1);
	}
	else
		tol = eps * eps * ddot(N, b, 1, b, 1);

	/* outer iteration */
	do
	{
		/* inner initialization */
		ax(A, x, r);
		daxpy(N, -1.0, b, 1, r, 1);
		if(pc)
			pc(A, r, w);
		y[0] = beta = dnrm2(N, w, 1);
		dcopy(N, w, 1, Q, 1);
		dscal(N, 1.0/beta, Q, 1);
		Hij = 0;

		/* inner iteration */
		for(j=0; j<M; ++j)
		{
			GIdouble *qj = Q + j*LDQ, *qj1 = Q + (j+1)*LDQ;

			/* compute and orthogonalize q[j+1] */
			if(pc)
			{
				ax(A, qj, w);
				pc(A, w, qj1);
			}
			else
				ax(A, qj, qj1);
			dgemv('T', N, j+1, 1.0, Q, LDQ, qj1, 1, 0.0, H+Hij, 1);
			dgemv('N', N, j+1, -1.0, Q, LDQ, H+Hij, 1, 1.0, qj1, 1);
			beta = dnrm2(N, qj1, 1);

			/* rotate new H-column */
			for(k=0; k<j; ++k,++Hij)
			{
				tmp = c[k]*H[Hij] - s[k]*H[Hij+1];
				H[Hij+1] = s[k]*H[Hij] + c[k]*H[Hij+1];
				H[Hij] = tmp;
			}

			/* compute new rotation */
			hjj = H[Hij];
			tmp = sqrt(hjj*hjj+beta*beta);
			c[j] = hjj / tmp;
			s[j] = -beta / tmp;
			H[Hij++] = tmp;

			/* rotate right hand side and normalize q[j+1] if needed further */
			y[j+1] = s[j] * y[j];
			y[j] = c[j] * y[j];
			if(fabs(y[j+1]) <= tol)
			{
				++j;
				break;
			}
			dscal(N, 1.0/beta, qj1, 1);
		}

		/* backward-eliminate for y and compute x */
		dtpsv('U', 'N', 'N', j, H, y, 1);
		dgemv('N', N, j, -1.0, Q, LDQ, y, 1, 1.0, x, 1);
		++i;
	}while(fabs(y[j]) > tol && i < max_iter);

	/* clean up */
	GI_FREE_ALIGNED(Q);
	GI_FREE_ALIGNED(H);
	GI_FREE_ALIGNED(r);
	GI_FREE_ALIGNED(c);
	GI_FREE_ALIGNED(s);
	GI_FREE_ALIGNED(y);
	if(pc)
		GI_FREE_ALIGNED(w);
	return (i-1)*M + j;
}

/** \internal
 *  \brief Tridiagonalize symmetric 3x3-matrix.
 *  \param mat symmetric 3x3-matrix in column-major format, contains transformation on return
 *  \param diag vector to take diagonal of tridiagonal matrix
 *  \param subdiag vector to take subdiagonal of tridiagonal matrix
 *  \retval GI_TRUE if transformation applied
 *  \retval GI_FALSE if already tridiagonal
 *  \ingroup numerics
 */
GIboolean GImat3d_tridiagonalize(GIdouble *mat, GIdouble *diag, 
								 GIdouble *subdiag)
{
	static const GIdouble identity[9] = { 
		1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
	diag[0] = mat[0];
	if(fabs(mat[6]) > DBL_EPSILON)
	{
		/* tridiagonalize */
		GIdouble tmp;
		subdiag[0] = sqrt(mat[3]*mat[3]+mat[6]*mat[6]);
		tmp = 1.0 / subdiag[0];
		mat[3] *= tmp;
		mat[6] *= tmp;
		tmp = 2.0*mat[3]*mat[7] + mat[6]*(mat[8]-mat[4]);
		diag[1] = mat[4] + mat[6]*tmp;
		diag[2] = mat[8] - mat[6]*tmp;
		subdiag[1] = mat[7] - mat[3]*tmp;
		mat[4] = mat[3];
		mat[7] = mat[6];
		mat[5] = mat[6];
		mat[8] = -mat[3];
		mat[0] = 1.0;
		mat[3] = mat[6] = mat[1] = mat[2] = 0.0;
		return GI_TRUE;
	}
	else
	{
		/* already tridiagonal */
		diag[1] = mat[4];
		diag[2] = mat[8];
		subdiag[0] = mat[3];
		subdiag[1] = mat[7];
		memcpy(mat, identity, 9*sizeof(GIdouble));
		return GI_FALSE;
	}
}

/** \internal
 *  \brief Apply QL algorithm to symmetric 3x3-tridiagonal matrix.
 *  \param mat transformation so far, contains eigenbasis on return
 *  \param diag diagonal of tridiagonal matrix, contains eigenvalues on return
 *  \param subdiag subdiagonal of tridiagonal matrix, zero on return
 *  \param max_iter maximum number of iterations
 *  \return number of iterations
 *  \ingroup numerics
 */
GIuint GImat3d_ql(GIdouble *mat, GIdouble *diag, 
				  GIdouble *subdiag, GIuint max_iter)
{
	GIuint k, i;
	for(k=1; k<=max_iter; ++k)
	{
		/* block-diagonal, lower 2x2 */
        GIdouble sum, diff, discr, eig0, eig1, cs, sn, 
			tmp, tmp0, tmp1, ratio, root, a, b;
        sum = fabs(diag[0]) + fabs(diag[1]);
		if(fabs(subdiag[0])+sum == sum)
		{
			/* compute eigenvalues as roots of quadratic equation */
			sum = fabs(diag[1]) + fabs(diag[2]);
			if(fabs(subdiag[1])+sum == sum)
				return k - 1;
			sum = diag[1] + diag[2];
			diff = diag[1] - diag[2];
			discr = sqrt(diff*diff+4.0*subdiag[1]*subdiag[1]);
			eig0 = 0.5 * (sum-discr);
			eig1 = 0.5 * (sum+discr);

			/* compute Givens rotation */
			if(diff >= 0.0)
			{
				cs = subdiag[1];
				sn = diag[1] - eig0;
			}
			else
			{
				cs = diag[2] - eig0;
				sn = subdiag[1];
			}
			tmp = 1.0 / sqrt(cs*cs+sn*sn);
			cs *= tmp;
			sn *= tmp;

			/* postmultiply current orthogonal matrix with Givens rotation */
			for(i=0; i<3; ++i)
			{
				tmp = mat[6+i];
				mat[6+i] = sn*mat[3+i] + cs*tmp;
				mat[3+i] = cs*mat[3+i] - sn*tmp;
			}

			/* update tridiagonal matrix */
			diag[1] = eig0;
			diag[2] = eig1;
			subdiag[0] = subdiag[1] = 0.0;
			break;
		}

		/* block-diagonal, upper 2x2 */
		sum = fabs(diag[1]) + fabs(diag[2]);
		if(fabs(subdiag[1])+sum == sum)
		{
			/* compute eigenvalues as roots of quadratic equation */
			sum = diag[0] + diag[1];
			diff = diag[0] - diag[1];
			discr = sqrt(diff*diff+4.0*subdiag[0]*subdiag[0]);
			eig0 = 0.5 * (sum-discr);
			eig1 = 0.5 * (sum+discr);

			/* compute Givens rotation */
			if(diff >= 0.0)
			{
				cs = subdiag[0];
				sn = diag[0] - eig0;
			}
			else
			{
				cs = diag[1] - eig0;
				sn = subdiag[0];
			}
			tmp = 1.0 / sqrt(cs*cs+sn*sn);
			cs *= tmp;
			sn *= tmp;

			/* postmultiply current orthogonal matrix with Givens rotation */
			for(i=0; i<3; ++i)
			{
				tmp = mat[3+i];
				mat[3+i] = sn*mat[i] + cs*tmp;
				mat[i] = cs*mat[i] - sn*tmp;
			}

			/* update tridiagonal matrix */
			diag[0] = eig0;
			diag[1] = eig1;
			subdiag[0] = subdiag[1] = 0.0;
			break;
		}

		/* set up parameters for first pass of QL step */
		ratio = (diag[1]-diag[0]) / (2.0*subdiag[0]);
		root = sqrt(1.0+ratio*ratio);
		b = subdiag[1];
		a = diag[2] - diag[0] + subdiag[0]/
			((ratio>=0.0) ? (ratio+root) : (ratio-root));

		/* compute Givens rotation for first pass */
		if(fabs(b) >= fabs(a))
		{
			ratio = a / b;
			sn = 1.0 / sqrt(1.0+ratio*ratio);
			cs = ratio * sn;
		}
		else
		{
			ratio = b / a;
			cs = 1.0 / sqrt(1.0+ratio*ratio);
			sn = ratio * cs;
		}

		/* postmultiply current orthogonal matrix with Givens rotation */
		for(i=0; i<3; ++i)
		{
			tmp = mat[6+i];
			mat[6+i] = sn*mat[3+i] + cs*tmp;
			mat[3+i] = cs*mat[3+i] - sn*tmp;
		}

		/* set up parameters for second pass of QL step */
		tmp0 = (diag[1]-diag[2])*sn + 2.0*subdiag[1]*cs;
		tmp1 = cs * subdiag[0];
		b = sn * subdiag[0];
		a = cs*tmp0 - subdiag[1];
		tmp0 *= sn;

		/* compute Givens rotation for first pass */
		if(fabs(b) >= fabs(a))
		{
			ratio = a / b;
			root = sqrt(1.0+ratio*ratio);
			subdiag[1] = b * root;
			sn = 1.0 / root;
			cs = ratio * sn;
		}
		else
		{
			ratio = b / a;
			root = sqrt(1.0+ratio*ratio);
			subdiag[1] = a * root;
			cs = 1.0 / root;
			sn = ratio * cs;
		}

		/* postmultiply current orthogonal matrix with Givens rotation */
		for(i=0; i<3; ++i)
		{
			tmp = mat[3+i];
			mat[3+i] = sn*mat[i] + cs*tmp;
			mat[i] = cs*mat[i] - sn*tmp;
		}

		/* update tridiagonal matrix */
		tmp = diag[1] - tmp0;
		diag[2] += tmp0;
		tmp0 = (diag[0]-tmp)*sn + 2.0*tmp1*cs;
		subdiag[0] = cs*tmp0 - tmp1;
		tmp0 *= sn;
		diag[1] = tmp + tmp0;
		diag[0] -= tmp0;
	}
	return k;
}
