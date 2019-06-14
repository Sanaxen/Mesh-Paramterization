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
 *  \brief Declaration of types and functions for multithreading.
 */

#ifndef __GI_THREAD_H__
#define __GI_THREAD_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include <stdlib.h>

#ifndef OPENGI_MAX_THREADS
	#define OPENGI_MAX_THREADS		4
#endif

#if defined(_WIN32) || HAVE_PTHREAD_H
	#define OPENGI_NUM_THREADS		OPENGI_MAX_THREADS
#else
	#define OPENGI_NUM_THREADS		0
#endif

#if OPENGI_NUM_THREADS > 1
	#ifdef _WIN32
		#ifndef _WIN32_WINNT
			#define _WIN32_WINNT	WINVER
		#endif
		#include <windows.h>
		#define GITHREADENTRY		__stdcall
	#else
		#include <pthread.h>
		#define GITHREADENTRY
	#endif
#else
	#define GITHREADENTRY
	typedef GIuint GIthread;
	typedef GIuint GIMutex;
	typedef GIuint GIthreadret;
#endif


/*************************************************************************/
/* Typedefs */

#if OPENGI_NUM_THREADS > 1

/** \internal
 *  \brief Thread handle.
 *  \ingroup threads
 */
#ifdef _WIN32
	typedef HANDLE		GIthread;
#else
	typedef pthread_t	GIthread;
#endif

/** \internal
 *  \brief Mutual exclusion.
 *  \ingroup threads
 */
#ifdef _WIN32
	typedef CRITICAL_SECTION	GIMutex;
#else
	typedef pthread_mutex_t		GIMutex;
#endif

/** \internal
 *  \brief Condition variable.
 *  \ingroup threads
 */
#ifdef _WIN32
	typedef struct _GICondVar
	{
		GIuint waiters;						/**< Number of waiting threads. */
		CRITICAL_SECTION waiters_mutex;		/**< Mutex for waiters variable. */
		HANDLE events[2];					/**< Events for waking up. */
	} GICondVar;
#else
	typedef pthread_cond_t GICondVar;
#endif

/** \internal
 *  \brief Synchronization barrier.
 *  \ingroup threads
 */
typedef struct _GIBarrier
{
	GIuint		num_threads;				/**< Needed number of threads. */
	GIuint		counter;					/**< Current number of threads. */
	GIMutex		mutex;						/**< Mutual exclusion. */
	GICondVar	cond_var;					/**< Broadcast to signal completeness. */
} GIBarrier;

/** \internal
 *  \brief Thread function return type.
 *  \ingroup threads
 */
#ifdef _WIN32
	typedef unsigned int	GIthreadret;
#else
	typedef void*			GIthreadret;
#endif

/** \internal
 *  \brief Thread execution function.
 *  \ingroup threads
 */
typedef GIthreadret (GITHREADENTRY *GIthreadfunc)(void*);

extern volatile GIuint g_uiActiveThreads;


/*************************************************************************/
/* Functions */

/** \name Thread methods
 *  \{
 */
GIthread GIthread_create(GIthreadfunc fn, GIvoid* arg);
GIthreadret GIthread_join(GIthread thread);
/** \} */

/** \name Mutex methods
 *  \{
 */
void GIMutex_construct(GIMutex *mutex);
void GIMutex_destruct(GIMutex *mutex);
void GIMutex_lock(GIMutex *mutex);
void GIMutex_unlock(GIMutex *mutex);
GIboolean GIMutex_trylock(GIMutex *mutex);
/** \} */

/** \name Condition variable methods
 *  \{
 */
void GICondVar_construct(GICondVar *cond);
void GICondVar_destruct(GICondVar *cond);
void GICondVar_wait(GICondVar *cond, GIMutex *mutex);
GIboolean GICondVar_timedwait(GICondVar *cond, GIMutex *mutex, GIuint msecs);
void GICondVar_wake_one(GICondVar *cond);
void GICondVar_wake_all(GICondVar *cond);
/** \} */

/** \name Barrier methods
 *  \{
 */
void GIBarrier_construct(GIBarrier *barrier, GIuint num_threads);
void GIBarrier_destruct(GIBarrier *barrier);
void GIBarrier_enter(GIBarrier *barrier);
/** \} */

#endif	/* OPENGI_NUM_THREADS > 1 */


#endif
