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
 *  \brief Declaration of structures and functions for memory management.
 */

#ifndef __GI_MEMORY_H__
#define __GI_MEMORY_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_thread.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef OPENGI_SSE
	#define OPENGI_SSE				0
#endif
#if OPENGI_SSE > 0
	#define GI_SSE_ALIGN_FLOAT		16
	#define GI_SSE_ALIGN_DOUBLE		16
	#define GI_SSE_SIZE(s)			(((s)+15)&(~15))
#else
	#define GI_SSE_ALIGN_FLOAT		4
	#define GI_SSE_ALIGN_DOUBLE		8
	#define GI_SSE_SIZE(s)			(s)
#endif


/*************************************************************************/
/* Macros */

/** \internal
 *  \brief Alignment of automatic variables.
 *  \ingroup memory
 */
#if defined(_MSC_VER) || defined(__ICL)
	#define GI_ALIGNED(v,a)		__declspec(align(a)) v
#elif defined(__GNU_C__)
	#define GI_ALIGNED(v,a)		v __attribute__((aligned(a)))
#else
	#define GI_ALIGNED(v,a)		v
#endif

/** \internal
 *  \brief Allocate memory for an array of objects.
 *  \ingroup memory
 */
#define GI_MALLOC_ARRAY(n,s)		malloc((n)*(s))

/** \internal
 *  \brief Allocate and clear memory for an array of objects.
 *  \ingroup memory
 */
#define GI_CALLOC_ARRAY(n,s)		calloc(n, s)

/** \internal
 *  \brief Reallocate memory for an array of objects.
 *  \ingroup memory
 */
#define GI_REALLOC_ARRAY(p,n,s)		realloc(p, (n)*(s))

/** \internal
 *  \brief Free memory of an array.
 *  \ingroup memory
 */
#define GI_FREE_ARRAY(p)			free(p)

/** \internal
 *  \brief Allocate memory for a single object.
 *  \ingroup memory
 */
//#define GI_MALLOC_SINGLE(s)			malloc(s)
#define GI_MALLOC_SINGLE(s)			GISmallObjectAllocator_allocate(&g_SmallObjAlloc, s)

/** \internal
 *  \brief Allocate and clear memory for a single object.
 *  \ingroup memory
 */
//#define GI_CALLOC_SINGLE(s)			calloc(1, s)
#define GI_CALLOC_SINGLE(s)			GISmallObjectAllocator_callocate(&g_SmallObjAlloc, s)

/** \internal
 *  \brief Free memory of a single object.
 *  \ingroup memory
 */
//#define GI_FREE_SINGLE(p,s)			free(p)
#define GI_FREE_SINGLE(p,s)			GISmallObjectAllocator_deallocate(&g_SmallObjAlloc, p, s)

/** \internal
 *  \brief Allocate memory for a mesh element.
 *  \ingroup memory
 */
//#define GI_MALLOC_PERSISTENT(s)		malloc(s)
#define GI_MALLOC_PERSISTENT(s)		GISmallObjectAllocator_allocate(&g_PersistentAlloc, s)

/** \internal
 *  \brief Allocate and clear memory for a mesh element.
 *  \ingroup memory
 */
//#define GI_CALLOC_PERSISTENT(s)		calloc(1, s)
#define GI_CALLOC_PERSISTENT(s)		GISmallObjectAllocator_callocate(&g_PersistentAlloc, s)

/** \internal
 *  \brief Free memory of a mesh element.
 *  \ingroup memory
 */
//#define GI_FREE_PERSISTENT(p,s)		free(p)
#define GI_FREE_PERSISTENT(p,s)		GISmallObjectAllocator_deallocate(&g_PersistentAlloc, p, s)

/** \internal
 *  \brief Allocate memory aligned to specififc byte boundary.
 *  \ingroup memory
 */
#define GI_MALLOC_ALIGNED(s,a)	GIMemory_allocate_aligned(s, a)

/** \internal
 *  \brief Allocate and clear memory aligned to specififc byte boundary.
 *  \ingroup memory
 */
#define GI_CALLOC_ALIGNED(s,a)	GIMemory_callocate_aligned(s, a)

/** \internal
 *  \brief Free memory aligned to specififc byte boundary.
 *  \ingroup memory
 */
#define GI_FREE_ALIGNED(p)		GIMemory_deallocate_aligned(p)


/*************************************************************************/
/* Typedefs */

typedef size_t GIusize;


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Chunk of memory.
 *  \details This class represents a fixed sized chunk of memory containing blocks 
 *  of fixed size. It is actually a C-port from Andrei Alexandrescu's Loki library.
 *  \ingroup memory
 */
typedef struct _GIChunk
{
	GIubyte	*data;							/**< Actual memory chunk. */
	GIubyte	first_free_block;				/**< Head of free block list. */
	GIubyte	num_free_blocks;				/**< Number of free blocks in chunk. */
} GIChunk;

/** \internal
 *  \brief Allocator for fixed sized blocks.
 *  \details This class represents an allocator for memory blocks of a fixed 
 *  size. It is actually a C-port from Andrei Alexandrescu's Loki library.
 *  \ingroup memory
 */
typedef struct _GIFixedAllocator
{
	GIuint	block_size;						/**< Block size this allocator serves. */
	GIubyte	num_blocks;						/**< Number of blocks per chunk. */
	GIusize	num_chunks;						/**< Number of chunks. */
	GIusize	max_chunks;						/**< Maximum number of chunks (till array grow). */
	GIusize	first_free_chunk;				/**< First chunk with free blocks (actually <= first free chunk). */
	GIChunk	*chunks;						/**< Chunk array. */
	GIChunk	*alloc_chunk;					/**< Last chunk used for allocation. */
	GIChunk	*dealloc_chunk;					/**< Last chunk used for deallocation. */
	GIChunk *empty_chunk;					/**< Only empty or NULL if none empty. */
} GIFixedAllocator;

/** \internal
 *  \brief Allocator for small objects.
 *  \details This class represents an allocator for single small objects. 
 *  It is actually a C-port from Andrei Alexandrescu's Loki library.
 *  \ingroup memory
 */
typedef struct _GISmallObjectAllocator
{
	GIFixedAllocator	*pool;				/**< Allocators for served block sizes. */
	GIMutex				mutex;				/**< Mutex for thread-safety */
} GISmallObjectAllocator;

extern GISmallObjectAllocator g_SmallObjAlloc;

extern GISmallObjectAllocator g_PersistentAlloc;


/*************************************************************************/
/* Functions */

/** \name Aligned allocation
 *  \{
 */
GIvoid* GIMemory_allocate_aligned(GIusize size, GIusize align);
GIvoid* GIMemory_callocate_aligned(GIusize size, GIusize align);
void GIMemory_deallocate_aligned(GIvoid *address);
/** \} */

/** \name Chunk methods
 *  \{
 */
GIboolean GIChunk_construct(GIChunk *chunk, GIuint block_size, GIuint blocks);
void GIChunk_destruct(GIChunk *chunk);
GIvoid* GIChunk_allocate(GIChunk *chunk, GIuint block_size);
void GIChunk_deallocate(GIChunk *chunk, GIvoid *address, GIuint block_size);
/** \} */

/** \name Fixed allocator methods
 *  \{
 */
void GIFixedAllocator_construct(GIFixedAllocator *alloc, GIuint block_size);
void GIFixedAllocator_destruct(GIFixedAllocator *alloc);
GIvoid* GIFixedAllocator_allocate(GIFixedAllocator *alloc);
void GIFixedAllocator_deallocate(GIFixedAllocator *alloc, GIvoid *address);
GIChunk* GIFixedAllocator_find(GIFixedAllocator *alloc, GIvoid *address);
GIboolean GIFixedAllocator_truncate(GIFixedAllocator *alloc);
/** \} */

/** \name Small object allocator methods
 *  \{
 */
void GISmallObjectAllocator_construct(GISmallObjectAllocator *alloc);
void GISmallObjectAllocator_destruct(GISmallObjectAllocator *alloc);
GIvoid* GISmallObjectAllocator_allocate(GISmallObjectAllocator *alloc, GIusize size);
GIvoid* GISmallObjectAllocator_callocate(GISmallObjectAllocator *alloc, GIusize size);
void GISmallObjectAllocator_deallocate(GISmallObjectAllocator *alloc, GIvoid *address, GIusize size);
GIboolean GISmallObjectAllocator_truncate(GISmallObjectAllocator *alloc);
void GISmallObjectAllocator_print(GISmallObjectAllocator *alloc, FILE *file);
/** \} */


#endif
