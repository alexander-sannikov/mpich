/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2017 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef MPIC_TYPES_H_INCLUDED
#define MPIC_TYPES_H_INCLUDED

#ifdef HAVE_SYS_QUEUE_H
#include "sys/queue.h"
#endif /* HAVE_SYS_QUEUE_H */

#define MPIC_FLAG_REDUCE_L 1
#define MPIC_FLAG_REDUCE_R 0

struct MPIR_Comm;
struct MPIR_Op;

#ifndef container_of
#define container_of(P,T,F) \
    ((T*) ((char *)P - offsetof(T,F)))
#endif

#define MPIC_Assert MPIR_Assert

#define MPIC_COMM_DECL MPIC_comm_t coll;
#define MPIC_REQ_DECL MPIC_req_t coll;

#define MPIC_COMM(comm) (&((comm)->coll))
#define MPIC_REQ(req) (&((req)->coll))

#ifdef MPIC_DEBUG
#define MPIC_DBG(...) do {MPL_DBG_MSG_FMT(MPIR_DBG_COLL,VERBOSE,(MPL_DBG_FDEST,__VA_ARGS__));} while (0)
#define MPIC_INLINE static
#else
#define MPIC_DBG(...)
#define MPIC_INLINE MPL_STATIC_INLINE_PREFIX
#endif

#include "coll_sched_db.h"

/* Queue element data structure in collectives queue.
 * Only non-blocking collectives are put on the queue */
typedef struct MPIC_queue_elem_t {
    TAILQ_ENTRY(MPIC_queue_elem_t) list_data;
    int (*kick_fn) (struct MPIC_queue_elem_t *); /* fn to make progress on the collective */
} MPIC_queue_elem_t;

/* Collectives request data structure */
typedef struct MPIC_req_t {
    MPIC_queue_elem_t elem;
    void *sched; /* pointer to the schedule */
} MPIC_req_t;

typedef struct MPIC_progress_global_t {
    TAILQ_HEAD(MPIC_queue_t, MPIC_queue_elem_t) head;
    int (*progress_fn) (void);  /* function to make progress on the collectives queue */
} MPIC_progress_global_t;


/* Coll 'pre'-definitions */
#include "types_decl.h"

/* Generic datatypes */
/* collective algorithm communicators */
typedef struct {
    int use_tag;
    MPIC_STUB_STUB_comm_t stub_stub;
    MPIC_STUB_TREE_comm_t stub_tree;
    MPIC_MPICH_STUB_comm_t mpich_stub;
    MPIC_MPICH_TREE_comm_t mpich_tree;
} MPIC_comm_t;

/* global data for every algorithm communicator */
typedef struct {
    MPIC_STUB_global_t tsp_stub;
    MPIC_MPICH_global_t tsp_mpich;

    MPIC_STUB_STUB_global_t stub_stub;
    MPIC_STUB_TREE_global_t stub_tree;
    MPIC_MPICH_STUB_global_t mpich_stub;
    MPIC_MPICH_TREE_global_t mpich_tree;
} MPIC_global_t;

extern MPIC_global_t MPIC_global_instance;

#endif /* MPIC_TYPES_H_INCLUDED */
