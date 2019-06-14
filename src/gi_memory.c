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
 *  \brief Implementation of structures and functions for memory management.
 */

#include "gi_memory.h"
#include "gi_math.h"

#include <assert.h>
#include <limits.h>

#define GI_CHUNKSIZE				2048
#define GI_MAX_BLOCKSIZE			256
#define GI_ALIGNMENT				4
#define GI_MIN_OBJPERCHUNK			8
#define GI_MAX_OBJPERCHUNK			(UCHAR_MAX+1)

#define GI_ALIGN_OFFSET(n)			((n+GI_ALIGNMENT-1)/GI_ALIGNMENT)

#define GI_CHUNK_CONTAINS(c,p,s)	((GIubyte*)p>=c->data && (GIubyte*)p<(c->data+s))


/** \internal
 *  \brief Global allocator object.
 *  \ingroup memory
 */
GISmallObjectAllocator g_SmallObjAlloc = { NULL };

/** \internal
 *  \brief Global allocator object for persistent.
 *  \ingroup memory
 */
GISmallObjectAllocator g_PersistentAlloc = { NULL };


/** \internal
 *  \brief Allocate memory aligned to specific boundary.
 *  \param size needed size of memory
 *  \param align alignment boundary (power of 2)
 *  \return pointer to aligned memory (do not call \c free on this)
 *  \ingroup memory
 */
GIvoid* GIMemory_allocate_aligned(GIusize size, GIusize align)
{
	GIvoid *pData, *pUnaligned;
	assert(align && GI_POWER_OF_2(align));
	if(!(pUnaligned=malloc(size+align-1+sizeof(GIvoid*))))
		return NULL;
	--align;
	pData = (GIvoid*)(((GIusize)pUnaligned+align+sizeof(GIvoid*)) & ~align);
	*((GIvoid**)pData-1) = pUnaligned;
	memset((GIubyte*)pData+size-align, 0, align);
	return pData;
}

/** \internal
 *  \brief Allocate and clear memory aligned to specific boundary.
 *  \param size needed size of memory
 *  \param align alignment boundary (power of 2)
 *  \return pointer to aligned memory (do not call \c free on this)
 *  \ingroup memory
 */
GIvoid* GIMemory_callocate_aligned(GIusize size, GIusize align)
{
	GIvoid *pData = GIMemory_allocate_aligned(size, align);
	if(pData)
		memset(pData, 0, size);
	return pData;
}

/** \internal
 *  \brief Free aligned memory block.
 *  \param address memory block to free
 *  \ingroup memory
 */
void GIMemory_deallocate_aligned(GIvoid *address)
{
	free(*((GIvoid**)address-1));
}

/** \internal
 *  \brief Chunk constructor.
 *  \param chunk chunk to construct
 *  \param block_size size of blocks
 *  \param blocks number of blocks
 *  \retval GI_TRUE if constructed successfully
 *  \retval GI_FALSE if out of memory
 *  \ingroup memory
 */
GIboolean GIChunk_construct(GIChunk *chunk, GIuint block_size, GIuint blocks)
{
	GIuint i;
	GIubyte *p = chunk->data = (GIubyte*)malloc(block_size*blocks*sizeof(GIubyte));
	if(!chunk->data)
		return GI_FALSE;

	/* init list of free blocks */
	chunk->first_free_block = 0;
	chunk->num_free_blocks = blocks;
	for(i=0; i<blocks; p+=block_size)
		*p = ++i;
	return GI_TRUE;
}

/** \internal
 *  \brief Chunk destructor.
 *  \param chunk chunk to destruct
 *  \ingroup memory
 */
void GIChunk_destruct(GIChunk *chunk)
{
	/* free chunk memory */
	assert(chunk->data);
	free(chunk->data);
}

/** \internal
 *  \brief Allocate memory block from chunk of fixed blocks.
 *  \param chunk chunk to take memory from
 *  \param block_size size of memory block to allocate
 *  \return pointer to allocated memory
 *  \ingroup memory
 */
GIvoid* GIChunk_allocate(GIChunk *chunk, GIuint block_size)
{
	GIubyte *pResult;

	/* get first free block and remove from list of free blocks */
	assert(chunk->num_free_blocks);
	pResult = chunk->data + chunk->first_free_block*block_size;
	chunk->first_free_block = *pResult;
	--chunk->num_free_blocks;
	return pResult;
}

/** \internal
 *  \brief Deallocate memory block from chunk of fixed blocks.
 *  \param chunk chunk to deallocate memory from
 *  \param address address of memory block to deallocate
 *  \param block_size size of memory block
 *  \ingroup memory
 */
void GIChunk_deallocate(GIChunk *chunk, GIvoid *address, GIuint block_size)
{
	GIubyte *pFree = (GIubyte*)address;

	/* add block to list of free blocks */
	assert(address >= chunk->data);
	assert((pFree-chunk->data)%block_size==0);
	*pFree = chunk->first_free_block;
	chunk->first_free_block = (pFree-chunk->data) / block_size;
	++chunk->num_free_blocks;
}

/** \internal
 *  \brief Fixed allocator constructor.
 *  \param alloc allocator to construct
 *  \param block_size size of fixed blocks
 *  \ingroup memory
 */
void GIFixedAllocator_construct(GIFixedAllocator *alloc, GIuint block_size)
{
	/* init allocator */
	alloc->block_size = block_size;
	alloc->num_blocks = GI_CHUNKSIZE / block_size;
	alloc->num_blocks = GI_CLAMP(alloc->num_blocks, 
		GI_MIN_OBJPERCHUNK, GI_MAX_OBJPERCHUNK);
	alloc->num_chunks = alloc->max_chunks = 0;
	alloc->first_free_chunk = 0;
	alloc->chunks = NULL;
	alloc->alloc_chunk = alloc->dealloc_chunk = alloc->empty_chunk = NULL;
}

/** \internal
 *  \brief Fixed allocator destructor.
 *  \param alloc allocator to destruct
 *  \ingroup memory
 */
void GIFixedAllocator_destruct(GIFixedAllocator *alloc)
{
	GIusize i;

	/* clean up */
	for(i=0; i<alloc->num_chunks; ++i)
		GIChunk_destruct(alloc->chunks+i);
	free(alloc->chunks);
	memset(alloc, 0, sizeof(GIFixedAllocator));
}

/** \internal
 *  \brief Allocate memory block from allocator of fixed block chunks.
 *  \param alloc allocator to take memory from
 *  \return pointer to allocated memory
 *  \ingroup memory
 */
GIvoid* GIFixedAllocator_allocate(GIFixedAllocator *alloc)
{
	/* find free chunk and allocate */
	if(!alloc->alloc_chunk || !alloc->alloc_chunk->num_free_blocks)
	{
		if(alloc->empty_chunk)
		{
			/* take empty chunk */
			alloc->alloc_chunk = alloc->empty_chunk;
			alloc->empty_chunk = NULL;
		}
		else
		{
			GIusize i = alloc->first_free_chunk;
			for(;; ++i)
			{
				if(i == alloc->num_chunks)
				{
					/* grow array if neccessary */
					if(alloc->num_chunks == alloc->max_chunks)
					{
						GIusize uiMaxChunks = alloc->num_chunks ? (alloc->num_chunks<<1) : 8;
						GIChunk *pChunks = (GIChunk*)realloc(alloc->chunks, uiMaxChunks*sizeof(GIChunk));
						if(!pChunks)
							return NULL;
						alloc->max_chunks = uiMaxChunks;
						alloc->chunks = pChunks;
					}

					/* init chunk and update cache */
					alloc->dealloc_chunk = alloc->chunks;
					alloc->alloc_chunk = alloc->chunks + i;
					if(!GIChunk_construct(alloc->alloc_chunk, alloc->block_size, alloc->num_blocks))
						return alloc->alloc_chunk = NULL;
					++alloc->num_chunks;
					break;
				}
				if(alloc->chunks[i].num_free_blocks)
				{
					/* update cache */
					alloc->alloc_chunk = alloc->chunks + i;
					break;
				}
			}
			alloc->first_free_chunk = i;
		}
	}
	else if(alloc->alloc_chunk == alloc->empty_chunk)
		alloc->empty_chunk = NULL;
	return GIChunk_allocate(alloc->alloc_chunk, alloc->block_size);
}

/** \internal
 *  \brief Deallocate memory block from allocator of fixed block chunks.
 *  \param alloc allocator to deallocate memory from
 *  \param address address of memory block to deallocate
 *  \ingroup memory
 */
void GIFixedAllocator_deallocate(GIFixedAllocator *alloc, GIvoid *address)
{

	/* find chunk with block and deallocate */
	assert(alloc->dealloc_chunk>=alloc->chunks && 
		alloc->dealloc_chunk<=(alloc->chunks+(alloc->num_chunks-1)));
	alloc->dealloc_chunk = GIFixedAllocator_find(alloc, address);
	assert(alloc->dealloc_chunk && alloc->dealloc_chunk->num_free_blocks<alloc->num_blocks);
	GIChunk_deallocate(alloc->dealloc_chunk, address, alloc->block_size);

	/* chunk now empty? */
	if(alloc->dealloc_chunk->num_free_blocks == alloc->num_blocks)
	{
		if(alloc->empty_chunk)
		{
			/* 2 empty chunks -> move one to end and release */
			GIChunk *pLast = alloc->chunks + (alloc->num_chunks-1);
			if(pLast == alloc->dealloc_chunk)
				alloc->dealloc_chunk = alloc->empty_chunk;
			else if(pLast != alloc->empty_chunk)
			{
				GIChunk temp;
				GI_SWAP(*pLast, *alloc->empty_chunk, temp);
			}
			assert(pLast->num_free_blocks == alloc->num_blocks);
			GIChunk_destruct(pLast);
			--alloc->num_chunks;
			if(alloc->alloc_chunk == pLast || !alloc->alloc_chunk->num_free_blocks)
				alloc->alloc_chunk = alloc->dealloc_chunk;
		}
		alloc->empty_chunk = alloc->dealloc_chunk;
	}
	else if(alloc->dealloc_chunk->num_free_blocks == 1)
	{
		/* chunk was full -> could be first free chunk now */
		GIusize i = alloc->dealloc_chunk - alloc->chunks;
		if(i < alloc->first_free_chunk)
			alloc->first_free_chunk = i;
	}
}

/** \internal
 *  \brief Find chunk to given memory address
 *  \param alloc fixed allocator to search in
 *  \param address address to look for
 *  \return chunk containing specified block or NULL if not found
 *  \ingroup memory
 */
GIChunk* GIFixedAllocator_find(GIFixedAllocator *alloc, GIvoid *address)
{
	GIChunk *pLo = alloc->dealloc_chunk, *pHi = alloc->dealloc_chunk+1, 
		*pLoBound = alloc->chunks, *pHiBound = alloc->chunks+alloc->num_chunks;
	GIuint uiChunkSize = alloc->block_size * alloc->num_blocks;

	/* bidirectional search starting at last deallocation */
	if(pHi == pHiBound)
		pHi = NULL;
	for(;;)
	{
		if(pLo)
		{
			if(GI_CHUNK_CONTAINS(pLo, address, uiChunkSize))
				return pLo;
			if(pLo == pLoBound)
			{
				pLo = NULL;
				if(!pHi)
					break;
			}
			else
				--pLo;
		}
		if(pHi)
		{
			if(GI_CHUNK_CONTAINS(pHi, address, uiChunkSize))
				return pHi;
			if(++pHi == pHiBound)
			{
				pHi = NULL;
				if(!pLo)
					break;
			}
		}
	}
	return NULL;
}

/** \internal
 *  \brief Truncate memory.
 *  \param alloc fixed allocator to truncate
 *  \retval GI_TRUE if truncated successfully
 *  \retval GI_FALSE if nothing to truncate
 *  \ingroup memory
 */
GIboolean GIFixedAllocator_truncate(GIFixedAllocator *alloc)
{
	if(alloc->empty_chunk)
	{
		/* release empty chunk */
		GIChunk *pLast = alloc->chunks + (alloc->num_chunks-1);
		if(alloc->empty_chunk != pLast)
		{
			GIChunk temp;
			GI_SWAP(*alloc->empty_chunk, *pLast, temp);
		}
		assert(pLast->num_free_blocks == alloc->num_blocks);
		GIChunk_destruct(pLast);
		--alloc->num_chunks;
		alloc->empty_chunk = NULL;
	}

	/* truncate if neccessary */
	if(alloc->num_chunks == alloc->max_chunks)
		return GI_FALSE;
	if(alloc->num_chunks)
	{
		/* truncate chunk array */
		GIChunk *pChunks = (GIChunk*)realloc(alloc->chunks, alloc->num_chunks*sizeof(GIChunk));
		if(pChunks)
			alloc->chunks = pChunks;
		alloc->dealloc_chunk = alloc->chunks;
		alloc->alloc_chunk = alloc->chunks + (alloc->num_chunks-1);
		if(!pChunks)
			return GI_FALSE;
	}
	else
	{
		/* free chunk array */
		free(alloc->chunks);
		alloc->chunks = alloc->alloc_chunk = alloc->dealloc_chunk = NULL;
	}
	alloc->max_chunks = alloc->num_chunks;
	return GI_TRUE;
}

/** \internal
 *  \brief Small object allocator constructor.
 *  \param alloc allocator to construct
 *  \ingroup memory
 */
void GISmallObjectAllocator_construct(GISmallObjectAllocator *alloc)
{
	GIuint i, uiNumAllocs = GI_ALIGN_OFFSET(GI_MAX_BLOCKSIZE);

	/* create allocators for fixed block sizes */
#if OPENGI_NUM_THREADS > 1
	GIMutex_construct(&alloc->mutex);
#endif
	alloc->pool = (GIFixedAllocator*)malloc(uiNumAllocs*sizeof(GIFixedAllocator));
	for(i=0; i<uiNumAllocs; ++i)
		GIFixedAllocator_construct(alloc->pool+i, (i+1)*GI_ALIGNMENT);
}

/** \internal
 *  \brief Small object allocator destructor.
 *  \param alloc allocator to destruct
 *  \ingroup memory
 */
void GISmallObjectAllocator_destruct(GISmallObjectAllocator *alloc)
{
	GIuint i, uiNumAllocs = GI_ALIGN_OFFSET(GI_MAX_BLOCKSIZE);

	/* clean up */
	for(i=0; i<uiNumAllocs; ++i)
		GIFixedAllocator_destruct(alloc->pool+i);
	free(alloc->pool);
	alloc->pool = NULL;
#if OPENGI_NUM_THREADS > 1
	GIMutex_destruct(&alloc->mutex);
#endif
}

/** \internal
 *  \brief Allocate memory for small object.
 *  \param alloc allocator to take memory from
 *  \param size size of memory block to allocate
 *  \return pointer to allocated memory or NULL if out of memory
 *  \ingroup memory
 */
GIvoid* GISmallObjectAllocator_allocate(GISmallObjectAllocator *alloc, 
										GIusize size)
{
	GIvoid *pResult;
	GIFixedAllocator *pAlloc;
	if(size > GI_MAX_BLOCKSIZE)
		return malloc(size);

	/* find allocator for block size and allocate */
	assert(size);
	pAlloc = alloc->pool + (GI_ALIGN_OFFSET(size)-1);
/*#if OPENGI_NUM_THREADS > 1
	if(g_uiActiveThreads > 1)
	{
		GIMutex_lock(&alloc->mutex);
		pResult = GIFixedAllocator_allocate(pAlloc);
		if(!pResult && GISmallObjectAllocator_truncate(alloc))
			pResult = GIFixedAllocator_allocate(pAlloc);
		GIMutex_unlock(&alloc->mutex);
	}
	else
#endif*/
	{
		pResult = GIFixedAllocator_allocate(pAlloc);
		if(!pResult && GISmallObjectAllocator_truncate(alloc))
			pResult = GIFixedAllocator_allocate(pAlloc);
	}
	return pResult;
}

/** \internal
 *  \brief Allocate and clear memory for small object.
 *  \param alloc allocator to take memory from
 *  \param size size of memory block to allocate
 *  \return pointer to allocated memory or NULL if out of memory
 *  \ingroup memory
 */
GIvoid* GISmallObjectAllocator_callocate(GISmallObjectAllocator *alloc, 
										 GIusize size)
{
	GIvoid *pResult;
	if(size > GI_MAX_BLOCKSIZE)
		return calloc(1, size);

	/* find allocate and init */
	pResult = GISmallObjectAllocator_allocate(alloc, size);
	memset(pResult, 0, size);
	return pResult;
}

/** \internal
 *  \brief Dellocate memory of small object.
 *  \param alloc allocator to deallocate memory from
 *  \param address address of memory block to deallocate
 *  \param size size of memory block to deallocate
 *  \ingroup memory
 */
void GISmallObjectAllocator_deallocate(GISmallObjectAllocator *alloc, 
									   GIvoid *address, GIusize size)
{
	GIFixedAllocator *pAlloc;
	assert(address);
	if(size > GI_MAX_BLOCKSIZE)
		free(address);
	else
	{
		/* find allocator for block size and deallocate */
		assert(size);
		pAlloc = alloc->pool + (GI_ALIGN_OFFSET(size)-1);
/*#if OPENGI_NUM_THREADS > 1
		if(g_uiActiveThreads > 1)
		{
			GIMutex_lock(&alloc->mutex);
			GIFixedAllocator_deallocate(pAlloc, address);
			GIMutex_unlock(&alloc->mutex);
		}
		else
#endif*/
			GIFixedAllocator_deallocate(pAlloc, address);
	}
}

/** \internal
 *  \brief Truncate memory.
 *  \param alloc allocator to truncate
 *  \retval GI_TRUE if truncated successfully
 *  \retval GI_FALSE if nothing to truncate
 *  \ingroup memory
 */
GIboolean GISmallObjectAllocator_truncate(GISmallObjectAllocator *alloc)
{
	GIuint i, uiNumAllocs = GI_ALIGN_OFFSET(GI_MAX_BLOCKSIZE);
	GIboolean bFound = GI_FALSE;

	/* truncate all fixed allocators */
/*#if OPENGI_NUM_THREADS > 1
	if(g_uiActiveThreads > 1)
	{
		GIMutex_lock(&alloc->mutex);
		for(i=0; i<uiNumAllocs; ++i)
			if(GIFixedAllocator_truncate(alloc->pool+i))
				bFound = GI_TRUE;
		GIMutex_unlock(&alloc->mutex);
	}
	else
#endif*/
		for(i=0; i<uiNumAllocs; ++i)
			if(GIFixedAllocator_truncate(alloc->pool+i))
				bFound = GI_TRUE;
	return bFound;
}

/** \internal
 *  \brief Print memory usage.
 *  \param alloc allocator to print
 *  \param file file to print to
 *  \ingroup memory
 */
void GISmallObjectAllocator_print(GISmallObjectAllocator *alloc, FILE *file)
{
	GIusize i, j, uiMem = 0, uiMaxMem = 0, uiNumAllocs = GI_ALIGN_OFFSET(GI_MAX_BLOCKSIZE);

#if OPENGI_NUM_THREADS > 1
//	GIMutex_lock(&alloc->mutex);
#endif
	/* print memory usage */
	for(i=0; i<uiNumAllocs; ++i)
	{
		GIFixedAllocator *pAlloc = alloc->pool + i;
		GIusize uiBlocks = 0;
		for(j=0; j<pAlloc->num_chunks; ++j)
			uiBlocks += pAlloc->num_blocks - pAlloc->chunks[j].num_free_blocks;
#ifdef _DEBUG
		if(pAlloc->num_chunks)
			fprintf(file, "%3dB: %d blocks, %f\n", pAlloc->block_size, uiBlocks, 
				pAlloc->num_chunks ? ((GIfloat)uiBlocks/
				(GIfloat)(pAlloc->num_chunks*pAlloc->num_blocks)) : 0);
#endif
		uiMem += uiBlocks * pAlloc->block_size;
		uiMaxMem += pAlloc->num_chunks * pAlloc->num_blocks * pAlloc->block_size;
	}
#if OPENGI_NUM_THREADS > 1
//	GIMutex_unlock(&alloc->mutex);
#endif
#ifdef _DEBUG
	fprintf(file, "overall: %dB, %f\n", uiMem, 
		uiMaxMem ? ((GIfloat)uiMem/(GIfloat)uiMaxMem) : 0.0f);
#endif
}
