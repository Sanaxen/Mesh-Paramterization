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
 *  \brief Declaration of generic container types.
 */

#ifndef __GI_CONTAINER_H__
#define __GI_CONTAINER_H__

#if HAVE_CONFIG_H
	#include <config.h>
#endif
#include <GI/gi.h>

#include "gi_memory.h"


/*************************************************************************/
/* Macros */

/** \internal
 *  \brief Insert item into doubly linked list.
 *  \ingroup container
 */
#define GI_LIST_INSERT(h,n,p)		if(h) { \
										(p)->prev=(n)->prev; (p)->next=(n); (n)->prev=(p); (p)->prev->next=(p); \
									} else { \
										(h)=(p)->next=(p)->prev=(p); }

/** \internal
 *  \brief Add item to back of doubly linked list.
 *  \ingroup container
 */
#define GI_LIST_ADD(h,p)			GI_LIST_INSERT(h,h,p)

/** \internal
 *  \brief Remove item from doubly linked list.
 *  \ingroup container
 */
#define GI_LIST_REMOVE(h,p)			if((h)==(h)->next) \
										(h)=NULL; \
									else if((p)==(h)) \
										(h)=(h)->next; \
									(p)->next->prev=(p)->prev; (p)->prev->next=(p)->next;

/** \internal
 *  \brief Remove item from doubly linked list and free its memory.
 *  \ingroup container
 */
#define GI_LIST_DELETE(h,p,s)		GI_LIST_REMOVE(h,p); GI_FREE_SINGLE(p,s);

/** \internal
 *  \brief Remove item from doubly linked list of persistent data and free its memory.
 *  \ingroup container
 */
#define GI_LIST_DELETE_PERSISTENT(h,p,s)	GI_LIST_REMOVE(h,p); GI_FREE_PERSISTENT(p,s);

/** \internal
 *  \brief Remove all items from doubly linked list and free their data.
 *  \ingroup container
 */
#define GI_LIST_CLEAR(h,s)			if(h) { \
										void *t, *e=(h); do { \
											t=(h); (h)=(h)->next; GI_FREE_SINGLE(t,s); \
										} while((h)!=e); (h)=NULL; }

/** \internal
 *  \brief Remove all items from doubly linked list of persistent data and free their data.
 *  \ingroup container
 */
#define GI_LIST_CLEAR_PERSISTENT(h,s)	if(h) { \
										void *t, *e=(h); do { \
											t=(h); (h)=(h)->next; GI_FREE_PERSISTENT(t,s); \
										} while((h)!=e); (h)=NULL; }

/** \internal
 *  \brief Begin iteration loop over doubly linked list.
 *  \ingroup container
 */
#define GI_LIST_FOREACH(h,p)		if(h) { (p)=(h); do {

/** \internal
 *  \brief End iteration loop over doubly linked list.
 *  \ingroup container
 */
#define GI_LIST_NEXT(h,p)			(p)=(p)->next; } while((p)!=(h)); }


/*************************************************************************/
/* Typedefs */

/** \internal
 *  \brief Generic hash function.
 *  \ingroup container
 */
typedef GIuint (*GIhashfunc)(const GIvoid*, GIuint);

/** \internal
 *  \brief Function for comparing keys in hash.
 *  \ingroup container
 */
typedef GIboolean (*GIcompfunc)(const GIvoid*, const GIvoid*);

/** \internal
 *  \brief Function for copying keys in hash.
 *  \ingroup container
 */
typedef void (*GIcopyfunc)(GIvoid*, const GIvoid*);

/** \internal
 *  \brief Function for comparing doubles.
 *  \ingroup container
 */
typedef GIint (*GIdcfunc)(GIdouble, GIdouble);


/*************************************************************************/
/* Structures */

/** \internal
 *  \brief Generic hash table list node.
 *  \details This structure represents a node in the collision list of a hash entry.
 *  \ingroup container
 */
typedef struct _GIHashNode
{
	GIvoid				*value;					/**< Value of item. */
	struct _GIHashNode	*next;					/**< Next in list. */
} GIHashNode;

/** \internal
 *  \brief Generic hash table.
 *  \details This structure represents a hash table associating generic data with pointers.
 *  \ingroup container
 */
typedef struct _GIHash
{
	GIHashNode	**data;							/**< Hashed items. */
	GIuint		key_size;						/**< Size of keys. */
	GIuint		size;							/**< Size of hash table. */
	GIuint		threshold;						/**< Size threshold. */
	GIuint		count;							/**< Number of items currently in hash. */
	GIfloat		load_factor;					/**< Maximum occupancy. */
	GIhashfunc	hash;							/**< Hash function to use. */
	GIcompfunc	comp;							/**< Function to compare keys with. */
	GIcopyfunc	copy;							/**< Function for copying keys. */
} GIHash;

/** \internal
 *  \brief Generic heap item.
 *  \details This structure represents an item of a heap-based priority queue.
 *  \ingroup container
 */
typedef struct _GIHeapItem
{
	GIdouble	priority;						/**< Priority of item. */
	GIvoid		*data;							/**< Data of item. */
} GIHeapItem;

/** \internal
 *  \brief Generic heap.
 *  \details This structure represents a generic heap based priority queue.
 *  \ingroup container
 */
typedef struct _GIHeap
{
	GIuint		size;							/**< Maximum number of elements. */
	GIuint		count;							/**< Number of elements in priority queue. */
	GIdcfunc	comp;							/**< Compare function. */
	GIHeapItem	*items;							/**< Items of priority queue. */
	GIboolean	store_pos;						/**< Store positions. */
	GIHash		pos_map;						/**< Positions of items in heap. */
} GIHeap;

/** \internal
 *  \brief Fibonacci tree.
 *  \details This structure represents a node of a Fibonacci tree.
 *  \ingroup container
 */
typedef struct _GIFibonacciNode
{
	GIdouble				priority;			/**< Priority of item. */
	GIvoid					*data;				/**< Data of item. */
	GIuint					ccount : 31;		/**< Number of children. */
	GIuint					flag : 1;			/**< Child removed flag. */
	struct _GIFibonacciNode	*parent;			/**< Parent node. */
	struct _GIFibonacciNode *children;			/**< List of children. */
	struct _GIFibonacciNode *next;				/**< Next child in list. */
	struct _GIFibonacciNode *prev;				/**< Previous child in list. */
} GIFibonacciNode;

/** \internal
 *  \brief Fibonacci heap.
 *  \details This structure represents a Fibonacci heap.
 */
typedef struct _GIFibonacciHeap
{
	GIFibonacciNode	*trees;						/**< List of Fibonacci trees. */
	GIFibonacciNode	*min_tree;						/**< Tree with minimum as root. */
	GIdcfunc		comp;						/**< Compare function. */
	GIHash			pos_map;					/**< Positions of items in heap. */
} GIFibonacciHeap;

/** \internal
 *  \brief Generic dynamic queue list node.
 *  \details This structure represents a node in the singly linked list of a LIFO/FIFO queue.
 *  \ingroup container
 */
typedef struct _GIQueueNode
{
	GIvoid				*data;					/**< Data of item. */
	struct _GIQueueNode	*next;					/**< Next list item. */
} GIQueueNode;

/** \internal
 *  \brief Generic dynamic queue.
 *  \details This structure represents a generic dynamic LIFO/FIFO queue.
 *  \ingroup container
 */
typedef struct _GIDynamicQueue
{
	GIuint		size;							/**< Number of elements in queue. */
	GIQueueNode	*head;							/**< First element in queue. */
	GIQueueNode	*tail;							/**< Last element in queue. */
} GIDynamicQueue;

/** \internal
 *  \brief Pair of unsigned ints.
 *  \details This structure represents a pair of unsigned integers, 
 *  mainly used for hashing.
 *  \ingroup container
 */
typedef struct _GIUintPair
{
	GIuint first;								/**< First element. */
	GIuint second;								/**< Second element. */
} GIUIntPair;


/*************************************************************************/
/* Functions */

/** \name Hash methods
 *  \{
 */
void GIHash_construct(GIHash *hash, GIuint size, GIfloat load, GIuint key_size, 
	GIhashfunc hash_func, GIcompfunc comp_func, GIcopyfunc copy_func);
void GIHash_destruct(GIHash *hash, GIuint value_size);
void GIHash_resize(GIHash *hash, GIuint size);
GIboolean GIHash_insert(GIHash *hash, const GIvoid *key, GIvoid *value);
GIvoid* GIHash_remove(GIHash *hash, const GIvoid *key);
void GIHash_clear(GIHash *hash, GIuint value_size);
GIvoid* GIHash_find(const GIHash *hash, const GIvoid *key);
/** \} */

/** \name Heap methods
 *  \{
 */
void GIHeap_construct(GIHeap *heap, GIuint size, 
	GIdcfunc comp, GIdouble sup, GIboolean store_pos);
void GIHeap_destruct(GIHeap *heap);
void GIHeap_resize(GIHeap *heap, GIuint size);
GIboolean GIHeap_enqueue(GIHeap *heap, GIvoid *data, GIdouble priority);
GIvoid* GIHeap_dequeue(GIHeap *heap, GIdouble *priority);
GIvoid* GIHeap_front(GIHeap *heap, GIdouble *priority);
GIboolean GIHeap_contains(GIHeap *heap, GIvoid *data);
GIboolean GIHeap_remove(GIHeap *heap, GIvoid *data);
void GIHeap_clear(GIHeap *heap);
/** \} */

/** \name Dynamic queue methods
 *  \{
 */
void GIDynamicQueue_construct(GIDynamicQueue *queue);
void GIDynamicQueue_enqueue(GIDynamicQueue *queue, GIvoid *data);
void GIDynamicQueue_push(GIDynamicQueue *queue, GIvoid *data);
GIvoid* GIDynamicQueue_dequeue(GIDynamicQueue *queue);
GIvoid* GIDynamicQueue_front(GIDynamicQueue *queue);
void GIDynamicQueue_clear(GIDynamicQueue *queue);
#define GIDynamicQueue_pop		GIDynamicQueue_dequeue
#define GIDynamicQueue_top		GIDynamicQueue_front
#define GIDynamicQueue_destruct	GIDynamicQueue_clear
/** \} */

/** \name Fibonacci heap methods
 *  \{
 */
void GIFibonacciHeap_construct(GIFibonacciHeap *heap, GIuint size, GIdcfunc comp);
void GIFibonacciHeap_destruct(GIFibonacciHeap *heap);
GIboolean GIFibonacciHeap_enqueue(GIFibonacciHeap *heap, GIvoid *data, GIdouble priority);
GIvoid* GIFibonacciHeap_dequeue(GIFibonacciHeap *heap, GIdouble *priority);
GIvoid* GIFibonacciHeap_front(GIFibonacciHeap *heap, GIdouble *priority);
GIboolean GIFibonacciHeap_contains(GIFibonacciHeap *heap, GIvoid *data);
GIboolean GIFibonacciHeap_remove(GIFibonacciHeap *heap, GIvoid *data);
void GIFibonacciHeap_clear(GIFibonacciHeap *heap);
/** \} */

/** \name Comparison functions
 *  \{
 */
GIint lessd(GIdouble a, GIdouble b);
GIint greaterd(GIdouble a, GIdouble b);
/** \} */

/** \name Functions for hashing
 *  \{
 */
GIuint hash_uint(const GIvoid *i, GIuint size);
void copy_uint(GIvoid *d, const GIvoid *s);
GIboolean compare_uint(const GIvoid *i, const GIvoid *j);
GIuint hash_pointer(const GIvoid *pointer, GIuint size);
void copy_pointer(GIvoid *d, const GIvoid *s);
GIboolean compare_pointer(const GIvoid *p, const GIvoid *q);
GIuint hash_vec3f(const GIvoid *v, GIuint size);
void copy_vec3f(GIvoid *d, const GIvoid *v);
GIboolean compare_vec3f(const GIvoid *v, const GIvoid *w);
GIuint hash_string(const GIvoid *str, GIuint size);
void copy_string(GIvoid *d, const GIvoid *s);
GIboolean compare_string(const GIvoid *s, const GIvoid *t);
GIuint hash_uintpair(const GIvoid *str, GIuint size);
void copy_uintpair(GIvoid *d, const GIvoid *s);
GIboolean compare_uintpair(const GIvoid *s, const GIvoid *t);
/** \} */


#endif
