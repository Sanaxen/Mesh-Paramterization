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
 *  \brief Implementation of generic container types.
 */

#include "gi_container.h"
#include "gi_math.h"

#include <stdlib.h>
#include <string.h>
#include <float.h>

/** \internal
 *  \brief nearest greater prime numbers for powers of two.
 *  \ingroup container
 */
static const GIuint g_primes[32] = {
	2, 
	5, 
	11, 
	17, 
	37, 
	67, 
	131, 
	257, 
	521, 
	1031, 
	2053, 
	4099, 
	8209, 
	16411, 
	32771, 
	65537, 
	131101, 
	262147, 
	524309, 
	1048583, 
	2097169, 
	4194319, 
	8388617, 
	16777259, 
	33554467, 
	67108879, 
	134217757, 
	268435459, 
	536870923, 
	1073741827, 
	2147483659U, 
	4294967291U
};


/** \internal
 *  \brief return next prime number.
 *  \param number number to search for
 *  \return nearest greater or equal POT prime number
 *  \ingroup container
 */
static GIuint next_prime(GIuint number)
{
	GIuint i;

	/* find nearest prime */
	for(i=0; i<32; ++i)
		if(g_primes[i] >= number)
			return g_primes[i];
	return number;
}

/** \internal
 *  \brief Hash constructor.
 *  \param hash hash to initialize
 *  \param size size of hash table
 *  \param key_size size of key type in bytes
 *  \param load maximum occupancy
 *  \param hash_func hash function to use
 *  \param comp_func function for comparing keys
 *  \param copy_func function for copying keys
 *  \ingroup container
 */
void GIHash_construct(GIHash *hash, GIuint size, GIfloat load, GIuint key_size, 
						   GIhashfunc hash_func, GIcompfunc comp_func, 
						   GIcopyfunc copy_func)
{
	if(load <= 0.0f)
		load = 0.8f;

	/* initialize hash */
	hash->key_size = key_size;
	hash->size = next_prime(size);
	hash->threshold = load * (GIfloat)hash->size;
	hash->count = 0;
	hash->load_factor = load;
	hash->data = (GIHashNode**)GI_CALLOC_ARRAY(hash->size, sizeof(GIHashNode*));
	hash->hash = hash_func;
	hash->comp = comp_func;
	hash->copy = copy_func;
}

/** \internal
 *  \brief Hash destructor.
 *  \param hash hash to destroy
 *  \param value_size size of value data to free or 0 for not freeing
 *  \ingroup container
 */
void GIHash_destruct(GIHash *hash, GIuint value_size)
{
	/* clean up */
	if(!hash->data)
		return;
	GIHash_clear(hash, value_size);
	GI_FREE_ARRAY(hash->data);
}

/** \internal
 *  \brief Resize hash.
 *  \param hash hash to resize
 *  \param size new size
 *  \ingroup container
 */
void GIHash_resize(GIHash *hash, GIuint size)
{
	GIHashNode **pData = hash->data, *pTemp;
	GIuint i, uiOldSize = hash->size;

	/* resize hash */
	hash->size = size;
	hash->threshold = hash->load_factor * (GIfloat)hash->size;
	hash->count = 0;
	hash->data = (GIHashNode**)GI_CALLOC_ARRAY(hash->size, sizeof(GIHashNode*));

	/* copy and free data */
	for(i=0; i<uiOldSize; ++i)
	{
		while(pData[i])
		{
			pTemp = pData[i];
			pData[i] = pTemp->next;
			GIHash_insert(hash, pTemp+1, pTemp->value);
			GI_FREE_SINGLE(pTemp, sizeof(GIHashNode)+hash->key_size);
		}
	}
	GI_FREE_ARRAY(pData);
}

/** \internal
 *  \brief Insert item into hash.
 *  \param hash hash to insert into
 *  \param key key of item
 *  \param value item to insert
 *  \retval GI_TRUE if inserted successfully
 *  \retval GI_FALSE if key already existing, value overwritten
 *  \ingroup container
 */
GIboolean GIHash_insert(GIHash *hash, const GIvoid *key, GIvoid *value)
{
	GIuint i;
	GIHashNode *pItem;
	if(!hash->data)
		return GI_FALSE;
	i = hash->hash(key, hash->size);
	pItem = hash->data[i];

	/* allready in hash -> overwrite */
	while(pItem)
	{
		if((*hash->comp)(pItem+1, key))
		{
			pItem->value = value;
			return GI_FALSE;
		}
		pItem = pItem->next;
	}

	/* create item and insert */
	pItem = (GIHashNode*)GI_MALLOC_SINGLE(sizeof(GIHashNode)+hash->key_size);
	(*hash->copy)(pItem+1, key);
	pItem->value = value;
	pItem->next = hash->data[i];
	hash->data[i] = pItem;

	/* grow if neccessary */
	if(++hash->count >= hash->threshold)
		GIHash_resize(hash, next_prime(hash->size+1));
	return GI_TRUE;
}

/** \internal
 *  \brief Remove item from hash.
 *  \param hash hash to remove from
 *  \param key key of item to remove
 *  \return removed item or NULL if not in hash
 *  \ingroup container
 */
GIvoid* GIHash_remove(GIHash *hash, const GIvoid *key)
{
	GIuint i;
	GIHashNode *pItem, *pNext;
	GIvoid *pValue;
	if(!hash->data)
		return NULL;
	i = hash->hash(key, hash->size);
	pItem = hash->data[i];
	if(!pItem)
		return NULL;

	/* remove at front */
	if((*hash->comp)(pItem+1, key))
	{
		hash->data[i] = pItem->next;
		pValue = pItem->value;
		GI_FREE_SINGLE(pItem, sizeof(GIHashNode)+hash->key_size);
		--hash->count;
		return pValue;
	}

	/* search item and remove */
	pNext = pItem->next;
	while(pNext && !(*hash->comp)(pNext+1, key))
	{
		pItem = pNext;
		pNext = pNext->next;
	}
	if(pNext)
	{
		pItem->next = pNext->next;
		pValue = pNext->value;
		GI_FREE_SINGLE(pNext, sizeof(GIHashNode)+hash->key_size);
		--hash->count;
		return pValue;
	}
	return NULL;
}

/** \internal
 *  \brief Remove all items from hash
 *  \param hash hash to clear
 *  \param value_size size of value data to free or 0 for not freeing
 *  \ingroup container
 */
void GIHash_clear(GIHash *hash, GIuint value_size)
{
	GIuint i;
	GIHashNode *pTemp;
	if(!hash->data)
		return;

	/* remove all items */
	for(i=0; i<hash->size && hash->count; ++i)
	{
		while(hash->data[i])
		{
			pTemp = hash->data[i];
			hash->data[i] = hash->data[i]->next;
			if(value_size)
				GI_FREE_SINGLE(pTemp->value, value_size);
			GI_FREE_SINGLE(pTemp, sizeof(GIHashNode)+hash->key_size);
			--hash->count;
		}
	}
}

/** \internal
 *  \brief Find item in hash.
 *  \param hash hash to search in
 *  \param key key of item to find
 *  \return item or NULL if not in hash
 *  \ingroup container
 */
GIvoid* GIHash_find(const GIHash *hash, const GIvoid *key)
{
	GIuint i;
	GIHashNode *pItem;
	if(!hash->data)
		return NULL;
	i = hash->hash(key, hash->size);
	pItem = hash->data[i];

	/* /search for item */
	while(pItem)
	{
		if((*hash->comp)(pItem+1, key))
			return pItem->value;
		pItem = pItem->next;
	}
	return NULL;
}

/** \internal
 *  \brief Heap constructor.
 *  \param heap heap to construct
 *  \param size maximum number of elements
 *  \param comp function to return GI_TRUE if first argument has higher priority
 *  \param sup lower bound for all keys (e.g. DBL_MIN)
 *  \param store_pos store positions of items
 *  \ingroup container
 */
void GIHeap_construct(GIHeap *heap, GIuint size, 
					  GIdcfunc comp, GIdouble sup, GIboolean store_pos)
{
	if(!comp)
	{
		comp = lessd;
		sup = -DBL_MAX;
	}

	/* initialize heap */
	heap->size = size;
	heap->count = 0;
	heap->comp = comp;
	heap->items = (GIHeapItem*)GI_MALLOC_ARRAY(size+1, sizeof(GIHeapItem));
	heap->items[0].priority = sup;
	heap->items[0].data = NULL;
	heap->store_pos = store_pos;
	if(store_pos)
		GIHash_construct(&heap->pos_map, size, 0.0f, 
			sizeof(GIvoid*), hash_pointer, compare_pointer, copy_pointer);
}

/** \internal
 *  \brief Heap destructor.
 *  \param heap heap to destruct
 *  \ingroup container
 */
void GIHeap_destruct(GIHeap *heap)
{
	/* clean up */
	if(heap->items)
		GI_FREE_ARRAY(heap->items);
	if(heap->store_pos)
		GIHash_destruct(&heap->pos_map, 0);
	memset(heap, 0, sizeof(GIHeap));
}

/** \internal
 *  \brief Resize heap.
 *  \param heap heap to resize
 *  \param size new size
 *  \ingroup container
 */
void GIHeap_resize(GIHeap *heap, GIuint size)
{
	if(!size)
		++size;
	else if(size <= heap->size)
		return;

	/* resize data */
	heap->items = (GIHeapItem*)GI_REALLOC_ARRAY(
		heap->items, size+1, sizeof(GIHeapItem));
	heap->size = size;
}

/** \internal
 *  \brief Enqueue element to heap.
 *  \param heap to insert into
 *  \param data data of item
 *  \param priority priority of item
 *  \retval GI_TRUE if item enqueued successfully
 *  \retval GI_FALSE if item already in heap
 *  \ingroup container
 */
GIboolean GIHeap_enqueue(GIHeap *heap, GIvoid *data, GIdouble priority)
{
	GIHeapItem *pItems = heap->items;
	GIdcfunc pfnComp = heap->comp;
	GIuint i, p;
	GIboolean bNew = GI_TRUE;

	/* heap full? */
	if(heap->count == heap->size)
		GIHeap_resize(heap, heap->size<<1);

	/* insert element and move up */
	if(heap->store_pos)
	{
		i = (GIuint)GIHash_find(&heap->pos_map, &data);
		if(i)
		{
			if(pfnComp(pItems[i].priority, priority) < 0)
			{
				GIHeap_remove(heap, data);
				i = ++heap->count;
			}
			bNew = GI_FALSE;
		}
		else
			i = ++heap->count;
		for(; pfnComp(pItems[i>>1].priority, priority)>0; i=p)
		{
			p = i >> 1;
			pItems[i] = pItems[p];
			GIHash_insert(&heap->pos_map, &pItems[i].data, (GIvoid*)i);
		}
		GIHash_insert(&heap->pos_map, &data, (GIvoid*)i);
	}
	else
	{
		for(i=++heap->count; pfnComp(pItems[i>>1].priority, priority)>0; i=p)
		{
			p = i >> 1;
			pItems[i] = pItems[p];
		}
	}
	pItems[i].data = data;
	pItems[i].priority = priority;
	return bNew;
}

/** \internal
 *  \brief Dequeue element with highest priority from heap.
 *  \param heap to get item from
 *  \param priority address to take priority of item or NULL if not wanted
 *  \return data of item with highest priority or NULL if heap empty
 *  \ingroup container
 */
GIvoid* GIHeap_dequeue(GIHeap *heap, GIdouble *priority)
{
	GIHeapItem *pItems = heap->items;
	GIdcfunc pfnComp = heap->comp;
	GIvoid *result;
	GIuint i, c, count = heap->count;
	GIdouble lastP = pItems[count].priority;

	/* heap empty? */
	if(!count)
	{
		if(priority)
			*priority = lastP;
		return pItems->data;
	}

	/* get first element */
	if(priority)
		*priority = pItems[1].priority;
	result = pItems[1].data;

	/* reinsert last element and move down */
	if(heap->store_pos)
	{
		GIHash_remove(&heap->pos_map, &result);
		for(i=1; (i<<1)<count; i=c)
		{
			c = i << 1;
			if(c+1 != count && pfnComp(pItems[c+1].priority, pItems[c].priority) < 0)
				++c;
			if(pfnComp(lastP, pItems[c].priority) <= 0)
				break;
			pItems[i] = pItems[c];
			GIHash_insert(&heap->pos_map, &pItems[i].data, (GIvoid*)i);
		}
		GIHash_insert(&heap->pos_map, &pItems[count].data, (GIvoid*)i);
	}
	else
	{
		for(i=1; (i<<1)<count; i=c)
		{
			c = i << 1;
			if(c+1 != count && pfnComp(pItems[c+1].priority, pItems[c].priority) < 0)
				++c;
			if(pfnComp(lastP, pItems[c].priority) <= 0)
				break;
			pItems[i] = pItems[c];
		}
	}
	pItems[i] = pItems[heap->count--];
	return result;
}

/** \internal
 *  \brief Query element with highest priority without dequeueing.
 *  \param heap to query
 *  \param priority address to take priority of item or NULL if not wanted
 *  \return data of item with highest priority or NULL if heap empty
 *  \ingroup container
 */
GIvoid* GIHeap_front(GIHeap *heap, GIdouble *priority)
{
	/* heap not empty? */
	if(heap->size)
	{
		if(priority)
			*priority = heap->items[1].priority;
		return heap->items[1].data;
	}

	/* heap empty */
	if(priority)
		*priority = heap->items->priority;
	return NULL;
}

/** \internal
 *  \brief Check if item is in heap.
 *  \param heap heap to query
 *  \param data item to look for
 *  \retval GI_TRUE if item in heap
 *  \retval GI_FALSE if item not in heap or position tracking not supported
 *  \ingroup container
 */
GIboolean GIHeap_contains(GIHeap *heap, GIvoid *data)
{
	/* look up item */
	if(heap->store_pos)
		return GIHash_find(&heap->pos_map, &data) != NULL;
	return GI_FALSE;
}

/** \internal
 *  \brief Remove item from heap.
 *  \param heap heap to remove from
 *  \param data item to remove
 *  \retval GI_TRUE if item removed successfully
 *  \retval GI_FALSE if not in heap or position tracking not supported
 *  \ingroup container
 */
GIboolean GIHeap_remove(GIHeap *heap, GIvoid *data)
{
	GIHeapItem *pItems = heap->items;
	GIdcfunc pfnComp = heap->comp;
	GIuint i, c, count = heap->count;
	GIdouble lastP = pItems[count].priority;

	/* item in heap? */
	if(!heap->store_pos || !(c=(GIuint)GIHash_remove(&heap->pos_map, &data)))
		return GI_FALSE;

	/* reinsert last element and move down */
	for(i=c; (i<<1)<count; i=c)
	{
		c = i << 1;
		if(c+1 != count && pfnComp(pItems[c+1].priority, pItems[c].priority) < 0)
			++c;
		if(pfnComp(lastP, pItems[c].priority) <= 0)
			break;
		pItems[i] = pItems[c];
		GIHash_insert(&heap->pos_map, &pItems[i].data, (GIvoid*)i);
	}
	GIHash_insert(&heap->pos_map, &pItems[count].data, (GIvoid*)i);
	pItems[i] = pItems[heap->count--];
	return GI_TRUE;
}

/** \internal
 *  \brief Remove all items from heap.
 *  \param heap heap to clear
 *  \ingroup container
 */
void GIHeap_clear(GIHeap *heap)
{
	/* no elements in heap */
	heap->count = 0;
	if(heap->store_pos)
		GIHash_clear(&heap->pos_map, 0);
}

/** \internal
 *  \brief Fibonacci tree destructor.
 *  \param node tree to destruct
 *  \ingroup container
 */
static void GIFibonacciNode_destruct(GIFibonacciNode *node)
{
	GIFibonacciNode *pNode;

	/* delete children and free */
	GI_LIST_FOREACH(node->children, pNode)
		GIFibonacciNode_destruct(pNode);
	GI_LIST_NEXT(node->children, pNode)
	GI_LIST_CLEAR(node->children, sizeof(GIFibonacciNode));
}

/** \internal
 *  \brief Cut and reinsert tree in Fibonacci heap.
 *  \param heap heap to work on
 *  \param node tree to cut out
 *  \ingroup container
 */
static void GIFibonacciHeap_cut(GIFibonacciHeap *heap, GIFibonacciNode *node)
{
	GIFibonacciNode *pParent = node->parent;

	/* cut subtree and meld lazy */
	GI_LIST_REMOVE(pParent->children, node);
	GI_LIST_ADD(heap->trees, node);
	node->parent = NULL;
	node->flag = 0;
	--pParent->ccount;

	/* cut parent if neccessary */
	if(pParent->parent)
	{
		if(pParent->flag)
			GIFibonacciHeap_cut(heap, pParent);
		else
			pParent->flag = 1;
	}
}

/** \internal
 *  \brief Fibonacci heap constructor.
 *  \param heap heap to construct
 *  \param size approximate maximum number of elements
 *  \param comp function to return GI_TRUE if first argument has higher priority
 *  \ingroup container
 */
void GIFibonacciHeap_construct(GIFibonacciHeap *heap, GIuint size, GIdcfunc comp)
{
	if(!comp)
		comp = lessd;

	/* initialize heap */
	heap->trees = heap->min_tree = NULL;
	heap->comp = comp;
	GIHash_construct(&heap->pos_map, size, 0.0f, sizeof(GIvoid*), 
		hash_pointer, compare_pointer, copy_pointer);
}

/** \internal
 *  \brief Fibonacci heap destructor.
 *  \param heap heap to destruct
 *  \ingroup container
 */
void GIFibonacciHeap_destruct(GIFibonacciHeap *heap)
{
	/* clean up */
	GIFibonacciHeap_clear(heap);
	GIHash_destruct(&heap->pos_map, 0);
	memset(heap, 0, sizeof(GIFibonacciHeap));
}

/** \internal
 *  \brief Enqueue element to Fibonacci heap.
 *  \param heap to insert into
 *  \param data data of item
 *  \param priority priority of item
 *  \retval GI_TRUE if item enqueued successfully
 *  \retval GI_FALSE if item already in heap
 *  \ingroup container
 */
GIboolean GIFibonacciHeap_enqueue(GIFibonacciHeap *heap, 
								  GIvoid *data, GIdouble priority)
{
	GIFibonacciNode *pNode = (GIFibonacciNode*)GIHash_find(&heap->pos_map, &data);
	GIdcfunc pfnComp = heap->comp;
	GIboolean bNew = GI_TRUE;

	/* item already in heap */
	if(pNode)
	{
		if(pfnComp(pNode->priority, priority) < 0)
		{
			/* remove and reinsert */
			GIFibonacciHeap_remove(heap, data);
			pNode = NULL;
		}
		else
		{
			/* update priority */
			pNode->priority = priority;
			if(pNode->parent && pfnComp(priority, pNode->parent->priority) < 0)
				GIFibonacciHeap_cut(heap, pNode);
		}
		bNew = GI_FALSE;
	}

	/* create and add new tree */
	if(pNode)
	{
		pNode = (GIFibonacciNode*)GI_CALLOC_SINGLE(sizeof(GIFibonacciNode));
		GI_LIST_ADD(heap->trees, pNode);
		pNode->data = data;
		pNode->priority = priority;
		GIHash_insert(&heap->pos_map, &data, pNode);
	}

	/* update min */
	if(!heap->min_tree || pfnComp(priority, heap->min_tree->priority) < 0)
		heap->min_tree = pNode;
	return bNew;
}

/** \internal
 *  \brief Dequeue element with highest priority from Fibonacci heap.
 *  \param heap to get item from
 *  \param priority address to take priority of item or NULL if not wanted
 *  \return data of item with highest priority or NULL if heap empty
 *  \ingroup container
 */
GIvoid* GIFibonacciHeap_dequeue(GIFibonacciHeap *heap, GIdouble *priority)
{
	GIFibonacciNode *pNode, *pNode2, *pMin = heap->min_tree;
	GIFibonacciNode *pTrees[32];
	GIdcfunc pfnComp = heap->comp;
	GIvoid *pData;
	GIuint i;
	if(!pMin)
		return NULL;

	/* get data and delete min */
	pData = pMin->data;
	if(priority)
		*priority = pMin->priority;
	GI_LIST_FOREACH(pMin->children, pNode)
		GI_LIST_ADD(heap->trees, pNode)
	GI_LIST_NEXT(pMin->children, pNode)
	GI_LIST_DELETE(heap->trees, pMin, sizeof(GIFibonacciNode));
	GIHash_remove(&heap->pos_map, &pData);

	/* compress heap */
	memset(pTrees, 0, 32*sizeof(GIFibonacciNode*));
	GI_LIST_FOREACH(heap->trees, pNode)
		pNode2 = pTrees[pNode->ccount];
		while(pNode2)
		{
			if(pfnComp(pNode->priority, pNode2->priority) > 0)
			{
				GI_SWAP(pNode, pNode2, pMin);
			}
			GI_LIST_ADD(pNode->children, pNode2)
			pNode2->parent = pNode;
			pTrees[pNode->ccount] = NULL;
			pNode2 = pTrees[++pNode->ccount];
		}
		pTrees[pNode->ccount] = pNode;
	GI_LIST_NEXT(heap->trees, pNode)

	/* relist and find min */
	heap->trees = heap->min_tree = NULL;
	for(i=0; i<32; ++i)
	{
		pNode = pTrees[i];
		if(pNode)
		{
			GI_LIST_ADD(heap->trees, pNode)
			if(!heap->min_tree || pfnComp(pNode->priority, heap->min_tree->priority) < 0)
				heap->min_tree = pNode;
		}
	}

	return pData;
}

/** \internal
 *  \brief Query element with highest priority without dequeueing.
 *  \param heap to query
 *  \param priority address to take priority of item or NULL if not wanted
 *  \return data of item with highest priority or NULL if heap empty
 *  \ingroup container
 */
GIvoid* GIFibonacciHeap_front(GIFibonacciHeap *heap, GIdouble *priority)
{
	/* return minimum if existent */
	if(heap->min_tree)
	{
		if(priority)
			*priority = heap->min_tree->priority;
		return heap->min_tree->data;
	}
	return NULL;
}

/** \internal
 *  \brief Check if item is in Fibonacci heap.
 *  \param heap heap to query
 *  \param data item to look for
 *  \retval GI_TRUE if item in heap
 *  \retval GI_FALSE if item not in heap or position tracking not supported
 *  \ingroup container
 */
GIboolean GIFibonacciHeap_contains(GIFibonacciHeap *heap, GIvoid *data)
{
	/* find item */
	return GIHash_find(&heap->pos_map, &data) != NULL;
}

/** \internal
 *  \brief Remove item from Fibonacci heap.
 *  \param heap heap to remove from
 *  \param data item to remove
 *  \retval GI_TRUE if item removed successfully
 *  \retval GI_FALSE if not in heap or position tracking not supported
 *  \ingroup container
 */
GIboolean GIFibonacciHeap_remove(GIFibonacciHeap *heap, GIvoid *data)
{
	GIFibonacciNode *pChild, *pNode = (GIFibonacciNode*)GIHash_remove(&heap->pos_map, &data);
	if(!pNode)
		return GI_FALSE;

	/* remove node and reinsert children */
	if(pNode->parent)
		GIFibonacciHeap_cut(heap, pNode);
	GI_LIST_FOREACH(pNode->children, pChild)
		GI_LIST_ADD(heap->trees, pChild)
	GI_LIST_NEXT(pNode->children, pChild)
	GI_LIST_DELETE(heap->trees, pNode, sizeof(GIFibonacciNode));
	return GI_TRUE;
}

/** \internal
 *  \brief Remove all items from Fibonacci heap.
 *  \param heap heap to clear
 *  \ingroup container
 */
void GIFibonacciHeap_clear(GIFibonacciHeap *heap)
{
	GIFibonacciNode *pNode;

	/* delete children and clear hash */
	GI_LIST_FOREACH(heap->trees, pNode)
		GIFibonacciNode_destruct(pNode);
	GI_LIST_NEXT(heap->trees, pNode)
	GI_LIST_CLEAR(heap->trees, sizeof(GIFibonacciNode));
	GIHash_clear(&heap->pos_map, 0);
}

/** \internal
 *  \brief Dynamic queue constructor.
 *  \param queue queue to construct
 *  \ingroup container
 */
void GIDynamicQueue_construct(GIDynamicQueue *queue)
{
	/* initialize empty queue */
	queue->size = 0;
	queue->head = queue->tail = NULL;
}

/** \internal
 *  \brief Insert item at end of queue
 *  \param queue queue to insert in
 *  \param data item to insert
 *  \ingroup container
 */
void GIDynamicQueue_enqueue(GIDynamicQueue *queue, GIvoid *data)
{
	GIQueueNode *pNode = (GIQueueNode*)GI_MALLOC_SINGLE(sizeof(GIQueueNode));
	pNode->data = data;
	pNode->next = NULL;

	/* insert element at end */
	if(queue->tail)
	{
		queue->tail->next = pNode;
		queue->tail = pNode;
	}
	else
		queue->head = queue->tail = pNode;
	++queue->size;
}

/** \internal
 *  \brief Insert item at front of queue
 *  \param queue queue to insert in
 *  \param data item to insert
 *  \ingroup container
 */
void GIDynamicQueue_push(GIDynamicQueue *queue, GIvoid *data)
{
	GIQueueNode *pNode = (GIQueueNode*)GI_MALLOC_SINGLE(sizeof(GIQueueNode));
	pNode->data = data;
	pNode->next = queue->head;

	/* insert element at front */
	queue->head = pNode;
	if(!queue->tail)
		queue->tail = pNode;
	++queue->size;
}

/** \internal
 *  \brief Dequeue item at front of queue
 *  \param queue queue to remove at
 *  \return item at front of queue
 *  \ingroup container
 */
GIvoid* GIDynamicQueue_dequeue(GIDynamicQueue *queue)
{
	GIQueueNode *pNode = queue->head;
	GIvoid *result = NULL;

	/* remove item if not empty */
	if(pNode)
	{
		result = pNode->data;
		queue->head = pNode->next;
		GI_FREE_SINGLE(pNode, sizeof(GIQueueNode));
		--queue->size;

		/* now empty? */
		if(!queue->head)
			queue->tail = NULL;
	}
	return result;
}

/** \internal
 *  \brief Query element at front without dequeueing.
 *  \param queue queue to query
 *  \return data of item at front of queue or NULL if queue empty
 *  \ingroup container
 */
GIvoid* GIDynamicQueue_front(GIDynamicQueue *queue)
{
	/* return head data */
	if(queue->head)
		return queue->head->data;
	return NULL;
}

/** \internal
 *  \brief Remove all items in queue
 *  \param queue queue to clear
 *  \ingroup container
 */
void GIDynamicQueue_clear(GIDynamicQueue *queue)
{
	GIQueueNode *pNode;

	/* remove all items */
	while(queue->head)
	{
		pNode = queue->head;
		queue->head = pNode->next;
		GI_FREE_SINGLE(pNode, sizeof(GIQueueNode));
	}
	queue->size = 0;
	queue->tail = NULL;
}

/** \internal
 *  \brief Compare doubles.
 *  \param a first argument
 *  \param b second argument
 *  \return negative value if a < b, positive value if a > b, 0 else
 *  \ingroup container
 */
GIint lessd(GIdouble a, GIdouble b)
{
	if(a < b)
		return -1;
	if(a > b)
		return 1;
	return 0;
}

/** \internal
 *  \brief Compare doubles.
 *  \param a first argument
 *  \param b second argument
 *  \return negative value if a > b, positive value if a < b, 0 else
 *  \ingroup container
 */
GIint greaterd(GIdouble a, GIdouble b)
{
	if(a > b)
		return -1;
	if(a < b)
		return 1;
	return 0;
}

/** \internal
 *  \brief Hash unsigned int.
 *  \param i value to hash
 *  \param size size of hash table
 *  \return computed hash value
 *  \ingroup container
 */
GIuint hash_uint(const GIvoid *i, GIuint size)
{
	return *(const GIuint*)i % size;
}

/** \internal
 *  \brief Copy unsigned ints.
 *  \param d destination
 *  \param s source
 */
void copy_uint(GIvoid *d, const GIvoid *s)
{
	*((GIuint*)d) = *((const GIuint*)s);
}

/** \internal
 *  \brief Compare unsigned ints.
 *  \param i first value
 *  \param j second value
 *  \retval GI_TRUE if uints are equal
 *  \retval GI_FALSE if uints are not equal
 *  \ingroup container
 */
GIboolean compare_uint(const GIvoid *i, const GIvoid *j)
{
	return *(const GIuint*)i == *(const GIuint*)j;
}

/** \internal
 *  \brief Hash pointer.
 *  \param pointer pointer to hash
 *  \param size size of hash table
 *  \return computed hash value
 *  \ingroup container
 */
GIuint hash_pointer(const GIvoid *pointer, GIuint size)
{
	return (*(const GIuint*)pointer >> 2) % size;
}

/** \internal
 *  \brief Copy pointer.
 *  \param d destination pointer
 *  \param s source pointer
 *  \ingroup container
 */
void copy_pointer(GIvoid *d, const GIvoid *s)
{
	*(GIvoid**)d = *(GIvoid**)s;
}

/** \internal
 *  \brief Compare unsigned ints.
 *  \param p first uint
 *  \param q second uint
 *  \retval GI_TRUE if uints are equal
 *  \retval GI_FALSE if uints are not equal
 *  \ingroup container
 */
GIboolean compare_pointer(const GIvoid *p, const GIvoid *q)
{
	return *(const GIvoid**)p == *(const GIvoid**)q;
}

/** \internal
 *  \brief Hash 3D vector.
 *  \param v vector to hash
 *  \param size size of hash table
 *  \return computed hash value
 *  \ingroup container
 */
GIuint hash_vec3f(const GIvoid *v, GIuint size)
{
	return (*((const GIuint*)v) ^ *((const GIuint*)v+1) ^ 
		*((const GIuint*)v+2)) % size;
}

/** \internal
 *  \brief Copy 3D vectors.
 *  \param d destination vector
 *  \param v source vector
 *  \ingroup container
 */
void copy_vec3f(GIvoid *d, const GIvoid *v)
{
	GI_VEC3_COPY((GIfloat*)d, (const GIfloat*)v);
}

/** \internal
 *  \brief Compare 3D vectors.
 *  \param v first vector
 *  \param w second vector
 *  \retval GI_TRUE if vectors are equal
 *  \retval GI_FALSE if vectors are not equal
 *  \ingroup container
 */
GIboolean compare_vec3f(const GIvoid *v, const GIvoid *w)
{
	return GI_VEC3_EQUAL((const GIfloat*)v, (const GIfloat*)w);
}

/** \internal
 *  \brief Hash string.
 *  \param str string to hash
 *  \param size size of hash table
 *  \return computed hash value
 *  \ingroup container
 */
GIuint hash_string(const GIvoid *str, GIuint size)
{
	const GIchar *c = (GIchar*)str;
	GIuint uiHash = 5381;
	while(*c)
		uiHash += (uiHash<<5) + *(c++);
/*	GIuint uiHash = 2166136261U;
	while(*c)
	{
		uiHash ^= (GIuint)(*(c++));
		uiHash *= 16777619;
	}
*/	return uiHash % size;
}

/** \internal
 *  \brief Copy strings.
 *  \param d destination string
 *  \param s source string
 *  \ingroup container
 */
void copy_string(GIvoid *d, const GIvoid *s)
{
	strcpy((char*)d, (const char*)s);
}

/** \internal
 *  \brief Compare strings.
 *  \param s first string
 *  \param t second string
 *  \retval GI_TRUE if strings are equal
 *  \retval GI_FALSE if strings are not equal
 *  \ingroup container
 */
GIboolean compare_string(const GIvoid *s, const GIvoid *t)
{
	return !strcmp((const char*)s, (const char*)t);
}

/** \internal
 *  \brief Hash unsigned int pair.
 *  \param pair pair to hash
 *  \param size size of hash table
 *  \return computed hash value
 *  \ingroup container
 */
GIuint hash_uintpair(const GIvoid *pair, GIuint size)
{
	return (((const GIUIntPair*)pair)->first + 
		((const GIUIntPair*)pair)->second) % size;
}

/** \internal
 *  \brief Copy unsigned int pair.
 *  \param d destination pair
 *  \param s source pair
 *  \ingroup container
 */
void copy_uintpair(GIvoid *d, const GIvoid *s)
{
	*(GIUIntPair*)d = *(const GIUIntPair*)s;
}

/** \internal
 *  \brief Compare 3D vectors.
 *  \param p first pair
 *  \param q second pair
 *  \retval GI_TRUE if pairs are equal
 *  \retval GI_FALSE if pairs are not equal
 *  \ingroup container
 */
GIboolean compare_uintpair(const GIvoid *p, const GIvoid *q)
{
	return (((const GIUIntPair*)p)->first == ((const GIUIntPair*)q)->first && 
		((const GIUIntPair*)p)->second == ((const GIUIntPair*)q)->second) || 
		(((const GIUIntPair*)p)->first == ((const GIUIntPair*)q)->second && 
		((const GIUIntPair*)p)->second == ((const GIUIntPair*)q)->first);
}
