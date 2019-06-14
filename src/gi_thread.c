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
 *  \brief Implementation of types and functions for multithreading.
 */

#include "gi_thread.h"
#include "gi_math.h"

#if OPENGI_NUM_THREADS > 1
	#ifdef _WIN32
		#include <process.h>
		#define GI_WAKE_ONE			0
		#define GI_WAKE_ALL			1
	#else
		#include <errno.h>
		#include <sys/time.h>
	#endif
#endif


#if OPENGI_NUM_THREADS > 1

/** \internal
 *  \brief Number of currently running threads.
 *  \ingroup threads
 */
volatile GIuint g_uiActiveThreads = 1;

/** \internal
 *  \brief Create and execute new thread.
 *  \param fn thread execution function
 *  \param arg arguments to thread function
 *  \return handle to newly created thread
 *  \ingroup threads
 */
GIthread GIthread_create(GIthreadfunc fn, GIvoid* arg)
{
	/* create thread */
	++g_uiActiveThreads;
#ifdef _WIN32
	return (GIthread)_beginthreadex(NULL, 0, fn, arg, 0, NULL);
#else
	GIthread thread;
	pthread_create(&thread, NULL, fn, arg);
	return thread;
#endif
}

/** \internal
 *  \brief Wiat for thread to exit.
 *  \param thread thread to wait for
 *  \return value returned by thread
 *  \ingroup threads
 */
GIthreadret GIthread_join(GIthread thread)
{
	GIthreadret retval;

	/* wait for thread */
#ifdef _WIN32
	WaitForSingleObject(thread, INFINITE);
	GetExitCodeThread(thread, &retval);
	CloseHandle(thread);
#else
	pthread_join(thread, &retval);
#endif

	/* return thread's return value */
	--g_uiActiveThreads;
	return retval;
}

/** \internal
 *  \brief Mutex constructor.
 *  \param mutex mutex to construct
 *  \ingroup threads
 */
void GIMutex_construct(GIMutex *mutex)
{
	/* create mutex */
#ifdef _WIN32
	InitializeCriticalSection(mutex);
#else
	pthread_mutex_init(mutex, NULL);
#endif
}

/** \internal
 *  \brief Mutex destructor.
 *  \param mutex mutex to destruct
 *  \ingroup threads
 */
void GIMutex_destruct(GIMutex *mutex)
{
	/* create mutex */
#ifdef _WIN32
	DeleteCriticalSection(mutex);
#else
	pthread_mutex_destroy(mutex);
#endif
}

/** \internal
 *  \brief Lock mutex.
 *  \param mutex mutex to lock
 *  \ingroup threads
 */
void GIMutex_lock(GIMutex *mutex)
{
	/* lock mutex */
#ifdef _WIN32
	EnterCriticalSection(mutex);
#else
	pthread_mutex_lock(mutex);
#endif
}

/** \internal
 *  \brief Unlock mutex.
 *  \param mutex mutex to unlock
 *  \ingroup threads
 */
void GIMutex_unlock(GIMutex *mutex)
{
	/* unlock mutex */
#ifdef _WIN32
	LeaveCriticalSection(mutex);
#else
	pthread_mutex_unlock(mutex);
#endif
}

/** \internal
 *  \brief Non-blocking mutex lock.
 *  \param mutex mutex to lock
 *  \retval GI_TRUE if mutex locked successfully
 *  \retval GI_FALSE if mutex allready locked
 *  \ingroup threads
 */
GIboolean GIMutex_trylock(GIMutex *mutex)
{
	/* try to lock mutex */
#ifdef _WIN32
	#if _WIN32_WINNT >= 0x0400
		return TryEnterCriticalSection(mutex) == TRUE;
	#else
		return GI_FALSE;
	#endif
#else
	return pthread_mutex_trylock(mutex) != EBUSY;
#endif
}

/** \internal
 *  \brief Condition variable constructor.
 *  \param cond condition variable to construct
 *  \ingroup threads
 */
void GICondVar_construct(GICondVar *cond)
{
	/* create condition variable */
#ifdef _WIN32
	cond->waiters = 0;
	InitializeCriticalSection(&cond->waiters_mutex);
	cond->events[GI_WAKE_ONE] = CreateEvent(NULL, FALSE, FALSE, NULL);
	cond->events[GI_WAKE_ALL] = CreateEvent(NULL, TRUE, FALSE, NULL);
#else
	pthread_cond_init(cond, NULL);
#endif
}

/** \internal
 *  \brief Condition variable destructor.
 *  \param cond condition variable to destruct
 *  \ingroup threads
 */
void GICondVar_destruct(GICondVar *cond)
{
	/* destroy condition variable */
#ifdef _WIN32
	GICondVar_wake_all(cond);
	DeleteCriticalSection(&cond->waiters_mutex);
	CloseHandle(cond->events[GI_WAKE_ONE]);
	CloseHandle(cond->events[GI_WAKE_ALL]);
#else
	pthread_cond_destroy(cond);
#endif
}

/** \internal
 *  \brief Wait for condition variable.
 *  \param cond condition variable to wait for
 *  \param mutex locked mutex
 *  \ingroup threads
 */
void GICondVar_wait(GICondVar *cond, GIMutex *mutex)
{
	/* wait for condition variable */
#ifdef _WIN32
	GICondVar_timedwait(cond, mutex, INFINITE);
#else
	pthread_cond_wait(cond, mutex);
#endif
}

/** \internal
 *  \brief Wait for condition variable or timeout.
 *  \param cond condition variable to wait for
 *  \param mutex locked mutex
 *  \param msecs timeout in milliseconds
 *  \retval GI_TRUE if condition signalled
 *  \retval GI_FALSE if waiting timed out
 *  \ingroup threads
 */
GIboolean GICondVar_timedwait(GICondVar *cond, GIMutex *mutex, GIuint msecs)
{
#ifdef _WIN32
	/* wait for events */
	DWORD dwResult;
	GIboolean bLastWaiter;
	EnterCriticalSection(&cond->waiters_mutex);
	++cond->waiters;
	LeaveCriticalSection(&cond->waiters_mutex);
	LeaveCriticalSection(mutex);
	dwResult = WaitForMultipleObjects(2, cond->events, FALSE, msecs);

	/* no more waiters -> manually reset broadcast event */
	EnterCriticalSection(&cond->waiters_mutex);
	--cond->waiters;
	bLastWaiter = dwResult == WAIT_OBJECT_0+GI_WAKE_ALL && !cond->waiters;
	LeaveCriticalSection(&cond->waiters_mutex);
	if(bLastWaiter)
		ResetEvent(cond->events[GI_WAKE_ALL]);
	EnterCriticalSection(mutex);
	return dwResult != WAIT_TIMEOUT;
#else
	/* wait for relative time */
	struct timespec timeout;
	clock_gettime(CLOCK_REALTIME, &timeout);
	timeout.tv_nsec += 1000*msecs;
	timeout.tv_sec += timeout.tv_nsec / 1000000;
	timeout.tv_nsec %= 1000000;
	return pthread_cond_timedwait(cond, mutex, &timeout) != ETIMEDOUT;
#endif
}

/** \internal
 *  \brief Signal a thread waiting for condition variable.
 *  \param cond condition variable to signal
 *  \ingroup threads
 */
void GICondVar_wake_one(GICondVar *cond)
{
	/* signal one thread */
#ifdef _WIN32
	GIboolean bWaiters;
	EnterCriticalSection(&cond->waiters_mutex);
	bWaiters = cond->waiters > 0;
	LeaveCriticalSection(&cond->waiters_mutex);
	if(bWaiters)
		SetEvent(cond->events[GI_WAKE_ONE]);
#else
	pthread_cond_signal(cond);
#endif
}

/** \internal
 *  \brief Signal all threads waiting for condition variable.
 *  \param cond condition variable to signal
 *  \ingroup threads
 */
void GICondVar_wake_all(GICondVar *cond)
{
	/* broadcast to all threads */
#ifdef _WIN32
	GIboolean bWaiters;
	EnterCriticalSection(&cond->waiters_mutex);
	bWaiters = cond->waiters > 0;
	LeaveCriticalSection(&cond->waiters_mutex);
	if(bWaiters)
		SetEvent(cond->events[GI_WAKE_ALL]);
#else
	pthread_cond_broadcast(cond);
#endif
}

/** \internal
 *  \brief Barrier constructor.
 *  \param barrier barrier to construct
 *  \param num_threads number of threads to wait for
 *  \ingroup threads
 */
void GIBarrier_construct(GIBarrier *barrier, GIuint num_threads)
{
	/* construct empty barrier */
	barrier->num_threads = num_threads;
	barrier->counter = 0;
	GIMutex_construct(&barrier->mutex);
	GICondVar_construct(&barrier->cond_var);
}

/** \internal
 *  \brief Barrier destructor.
 *  \param barrier barrier to destruct
 *  \ingroup threads
 */
void GIBarrier_destruct(GIBarrier *barrier)
{
	/* wake and clean up */
	GICondVar_wake_all(&barrier->cond_var);
	GIMutex_destruct(&barrier->mutex);
	GICondVar_destruct(&barrier->cond_var);
}

/** \internal
 *  \brief Enter barrier.
 *  \param barrier barrier to enter
 *  \ingroup threads
 */
void GIBarrier_enter(GIBarrier *barrier)
{
	/* enter barrier and wait till full */
#ifdef _WIN32
	EnterCriticalSection(&barrier->mutex);
#else
	pthread_mutex_lock(&barrier->mutex);
#endif
	if(++barrier->counter == barrier->num_threads)
	{
		barrier->counter = 0;
		GICondVar_wake_all(&barrier->cond_var);
	}
	else
		GICondVar_wait(&barrier->cond_var, &barrier->mutex);
#ifdef _WIN32
	LeaveCriticalSection(&barrier->mutex);
#else
	pthread_mutex_unlock(&barrier->mutex);
#endif
}

#endif	/* OPENGI_NUM_THREADS > 1 */
