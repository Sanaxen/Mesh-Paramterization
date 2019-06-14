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
 *  \brief Declaration of structures and functions for parameterizing triangular meshes.
 */

#ifndef __GI_PARAMETERIZER_H__
#define __GI_PARAMETERIZER_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_thread.h"
#include "gi_mesh.h"
#include "gi_cutter.h"
#include "gi_numerics.h"

#define GI_CALLBACK_BASE		GI_PARAM_STARTED
#define GI_CALLBACK_END			GI_PARAM_FINISHED
#define GI_CALLBACK_COUNT		(GI_CALLBACK_END-GI_CALLBACK_BASE+1)


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Parameterization configuration.
 *  \ingroup parameterization
 */
typedef struct _GIParameterizer
{
	struct _GIContext	*context;						/**< Context this parameterizer belongs to */
	GIenum				parameterizer;					/**< Parameterization algorithm to use. */
	GIenum				initial_param;					/**< Initial parameterization for stretch minimizer. */
	GIenum				stretch_metric;					/**< Metric to use for stretch minimization. */
	GIfloat				conformal_weight;				/**< Conformal weight for intrinsic parameterization. */
	GIfloat				authalic_weight;				/**< Authalic weight for intrinsic parameterization. */
	GIfloat				stretch_weight;					/**< Eta parameter for stretch minimizer. */
	GIfloat				area_weight;					/**< Theta parameter for combined energy. */ 
	GIuint				source_attrib;					/**< Attribute to use as parameter coordinates. */
	GIuint				sampling_res;					/**< Desired minimal sampling resolution. */
	GIenum				solver;							/**< Solver for unsymmetric systems. */
	GIparamcb			callback[GI_CALLBACK_COUNT];	/**< Callback function. */
	GIvoid				*cdata[GI_CALLBACK_COUNT];		/**< User data for callback function. */
} GIParameterizer;

/** \internal
 *  \brief Data for linear equation system.
 *  \ingroup parameterization
 */
typedef struct _GILinearSystem
{
	GIParameterizer		*parameterizer;			/**< Parameterizer to use. */
	GIPatch				*patch;					/**< Patch to which system belongs. */
	GISparseMatrixCSR	*A;						/**< Matrix of coefficients. */
	GISparseMatrixLIL	*B;						/**< Separately stored coefficients of right hand side. */
	GIdouble			*bU;					/**< Right hand side for U coordinate. */
	GIdouble			*bV;					/**< Right hand side for V coordinate. */
	GIdouble			*u;						/**< Unknown vector for U coordinate. */
	GIdouble			*v;						/**< Unknown vector for V coordinate. */
} GILinearSystem;

/** \internal
 *  \brief Parameters for solver thread.
 *  \ingroup parameterization
 */
typedef struct _GISolverData
{
	GIsolverfunc			solver_func;			/**< Solving function. */
	const GISparseMatrix	*A;						/**< System matrix. */
	const GIdouble			*b;						/**< Right hand side vector. */
	GIdouble				*x;						/**< Vector of unknowns. */
	GImvfunc				ax;						/**< Matrix-vector-multiplication function. */
	GImvfunc				pc;						/**< Prconditioning function. */
	GIdouble				eps;					/**< Error threshold. */
	GIuint					max_iter;				/**< Maximum number of iterations. */
} GISolverData;


/*************************************************************************/
/* Functions */

/** \name Parameterizer methods
 *  \{
 */
void GIParameterizer_construct(GIParameterizer *par, struct _GIContext *context);
GIboolean GIParameterizer_arc_length_circle(GIParameterizer *par, GIPatch *patch);
GIboolean GIParameterizer_arc_length_square(GIParameterizer *par, GIPatch *patch);
GIboolean GIParameterizer_stretch_minimizing(GIParameterizer *par, GIPatch *patch);
GIboolean GIParameterizer_stretch_minimizing2(GIParameterizer *par, GIPatch *patch);
GIboolean GIParameterizer_gim(GIParameterizer *par, GIPatch *patch);
/** \} */

/** \name Linear system methods
 *  \{
 */
void GILinearSystem_construct(GILinearSystem *system, GIParameterizer *par, 
	GIPatch *patch, GIenum type, GIboolean force_non_symmetric, GIboolean store_rhs);
void GILinearSystem_destruct(GILinearSystem *system);
void GILinearSystem_unknowns_to_params(GILinearSystem *system);
GIboolean GILinearSystem_solve(GILinearSystem *system);
GIthreadret GITHREADENTRY GILinearSystem_solve_thread(GIvoid *arg);
/** \} */


#endif
