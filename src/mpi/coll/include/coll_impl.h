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

#ifndef MPIC_COLL_IMPL_H_INCLUDED
#define MPIC_COLL_IMPL_H_INCLUDED

#include <sys/queue.h>

extern MPIC_progress_global_t MPIC_progress_global;
extern MPIC_global_t MPIC_global_instance;

#define GLOBAL_NAME    MPIC_

#include "../transports/stub/transport.h"
#include "../transports/mpich/transport.h"

#include "../transports/stub/api_def.h"
#include "../algorithms/tree/kary_post.h"
#include "../algorithms/tree/knomial_post.h"
#include "../src/tsp_namespace_undef.h"

#include "../transports/mpich/api_def.h"
#include "../algorithms/tree/kary_post.h"
#include "../algorithms/tree/knomial_post.h"
#include "../src/tsp_namespace_undef.h"

#undef GLOBAL_NAME

#define MPIC_INIT_DT_BUILTIN(X)              \
    MPID_Datatype_get_ptr(X, dt_ptr);   \
    MPIC_dt_init(dt_ptr);

#define MPIC_INIT_CXX_DT_BUILTIN(X)              \
    MPID_Datatype_get_ptr(X, dt_ptr);   \
    dt_ptr->handle=X;                   \
    MPIC_dt_init(dt_ptr);


#define MPIC_NUM_ENTRIES (8)
struct MPIR_Request;

MPL_STATIC_INLINE_PREFIX int MPIR_Coll_cycle_algorithm(MPIR_Comm *comm_ptr,int pick[], int num) {
    int idx;
    if(comm_ptr->comm_kind == MPIR_COMM_KIND__INTERCOMM)
        return 0;
    MPIC_COMM(comm_ptr)->issued_collectives++;
    idx = MPIC_COMM(comm_ptr)->issued_collectives%num;
    return pick[idx];
}

static inline int MPIC_Progress(int n, void *cq[])
{
    int i = 0;
    COLL_queue_elem_t *s = MPIC_progress_global.head.tqh_first;
    for ( ; ((s != NULL) && (i < n)); s = s->list_data.tqe_next) {
        if(s->kick_fn(s))
            cq[i++] = (void *) s;
    }
    return i;
}

static inline int MPIC_progress_hook()
{
    void *coll_entries[MPIC_NUM_ENTRIES];
    int i, coll_count, mpi_errno;
    mpi_errno = MPI_SUCCESS;
    coll_count = MPIC_Progress(MPIC_NUM_ENTRIES, coll_entries);
    for (i = 0; i < coll_count; i++) {
        MPIC_req_t *base = (MPIC_req_t *) coll_entries[i];
        MPIR_Request *req = container_of(base, MPIR_Request, coll);
        MPID_Request_complete(req);
    }
    return mpi_errno;
}

extern int MPIC_comm_counter;
extern MPIC_sched_entry_t *MPIC_sched_table;
MPL_STATIC_INLINE_PREFIX int MPIC_comm_init(struct MPIR_Comm *comm)
{
    int rank, size, mpi_errno = MPI_SUCCESS;
    int *tag = &(MPIC_COMM(comm)->use_tag);
    *tag = 0;
    MPIC_COMM(comm)->issued_collectives = 0;
    rank = comm->rank;
    size = comm->local_size;
    /*initialize communicators for ch4 collectives
    * each communicator is assigned a unique id 
    */
    MPIC_STUB_STUB_comm_init(&(MPIC_COMM(comm)->stub_stub), MPIC_comm_counter++, tag,  rank, size);
    MPIC_STUB_KARY_comm_init(&(MPIC_COMM(comm)->stub_kary), MPIC_comm_counter++, tag, rank, size);
    MPIC_STUB_KNOMIAL_comm_init(&(MPIC_COMM(comm)->stub_knomial), MPIC_comm_counter++, tag,  rank, size);
    MPIC_STUB_RECEXCH_comm_init(&(MPIC_COMM(comm)->stub_recexch), MPIC_comm_counter++, tag,  rank, size);
    MPIC_STUB_DISSEM_comm_init(&(MPIC_COMM(comm)->stub_dissem), MPIC_comm_counter++, tag,  rank, size);
    MPIC_MPICH_STUB_comm_init(&(MPIC_COMM(comm)->mpich_stub), MPIC_comm_counter++, tag, rank, size);
    MPIC_MPICH_KARY_comm_init(&(MPIC_COMM(comm)->mpich_kary), MPIC_comm_counter++, tag, rank, size);
    MPIC_MPICH_KNOMIAL_comm_init(&(MPIC_COMM(comm)->mpich_knomial), MPIC_comm_counter++, tag, rank, size);
    MPIC_MPICH_RECEXCH_comm_init(&(MPIC_COMM(comm)->mpich_recexch), MPIC_comm_counter++, tag,  rank, size);
    MPIC_MPICH_DISSEM_comm_init(&(MPIC_COMM(comm)->mpich_dissem), MPIC_comm_counter++, tag,  rank, size);
    MPIC_BMPICH_KARY_comm_init(&(MPIC_COMM(comm)->bmpich_kary), MPIC_comm_counter++, tag, rank, size);
    MPIC_BMPICH_KNOMIAL_comm_init(&(MPIC_COMM(comm)->bmpich_knomial), MPIC_comm_counter++, tag, rank, size);
    MPIC_X_TREEBASIC_comm_init(&(MPIC_COMM(comm)->x_treebasic), MPIC_comm_counter++, tag, rank, size);

#ifdef MPID_COLL_COMM_INIT_HOOK
    MPID_COLL_COMM_INIT_HOOK;
#endif
    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX int MPIC_comm_cleanup(struct MPIR_Comm *comm)
{
    int mpi_errno = MPI_SUCCESS;
    /*cleanup netmod collective communicators*/
#ifdef MPID_COLL_COMM_CLEANUP_HOOK 
    MPID_COLL_COMM_CLEANUP_HOOK;
#endif
    /*cleanup all ch4 collective communicators*/
    MPIC_STUB_STUB_comm_cleanup(&(MPIC_COMM(comm)->stub_stub));
    MPIC_STUB_KARY_comm_cleanup(&(MPIC_COMM(comm)->stub_kary));
    MPIC_STUB_KNOMIAL_comm_cleanup(&(MPIC_COMM(comm)->stub_knomial));
    MPIC_STUB_RECEXCH_comm_cleanup(&(MPIC_COMM(comm)->stub_recexch));
    MPIC_STUB_DISSEM_comm_cleanup(&(MPIC_COMM(comm)->stub_dissem));

    MPIC_MPICH_STUB_comm_cleanup(&(MPIC_COMM(comm)->mpich_stub));
    MPIC_MPICH_KARY_comm_cleanup(&(MPIC_COMM(comm)->mpich_kary));
    MPIC_MPICH_KNOMIAL_comm_cleanup(&(MPIC_COMM(comm)->mpich_knomial));
    MPIC_MPICH_RECEXCH_comm_cleanup(&(MPIC_COMM(comm)->mpich_recexch));
    MPIC_MPICH_DISSEM_comm_cleanup(&(MPIC_COMM(comm)->mpich_dissem));

    MPIC_BMPICH_KARY_comm_cleanup(&(MPIC_COMM(comm)->bmpich_kary));
    MPIC_BMPICH_KNOMIAL_comm_cleanup(&(MPIC_COMM(comm)->bmpich_knomial));

    MPIC_X_TREEBASIC_comm_cleanup(&(MPIC_COMM(comm)->x_treebasic));

    return mpi_errno;
}


/* There is no initialization required at the MPI layer 
 * but we should still call these hooks in case the transports need to do 
 * their own initialization (e.g. UCX-netmod, OFI netmod).*/
MPL_STATIC_INLINE_PREFIX int MPIC_dt_init(MPIR_Datatype * dt)
{
    int mpi_errno = MPI_SUCCESS;
    /*initialize ch4 collective datatypes*/
#ifdef MPID_COLL_DT_INIT_HOOK
    MPID_COLL_DT_INIT_HOOK;
#endif

    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX int MPIC_op_init(struct MPIR_Op *op)
{
    int mpi_errno = MPI_SUCCESS;
    /*initialize ch4 collective operations*/
#ifdef MPID_COLL_OP_INIT_HOOK
    MPID_COLL_OP_INIT_HOOK;
#endif

    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX int MPIC_init_builtin_dt()
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Datatype *dt_ptr;

    MPIC_INIT_DT_BUILTIN(MPI_CHAR);
    MPIC_INIT_DT_BUILTIN(MPI_UNSIGNED_CHAR);
    MPIC_INIT_DT_BUILTIN(MPI_SIGNED_CHAR);
    MPIC_INIT_DT_BUILTIN(MPI_BYTE);
    MPIC_INIT_DT_BUILTIN(MPI_WCHAR);
    MPIC_INIT_DT_BUILTIN(MPI_SHORT);
    MPIC_INIT_DT_BUILTIN(MPI_UNSIGNED_SHORT);
    MPIC_INIT_DT_BUILTIN(MPI_INT);
    MPIC_INIT_DT_BUILTIN(MPI_UNSIGNED);
    MPIC_INIT_DT_BUILTIN(MPI_LONG);
    MPIC_INIT_DT_BUILTIN(MPI_UNSIGNED_LONG);       /* 10 */

    MPIC_INIT_DT_BUILTIN(MPI_FLOAT);
    MPIC_INIT_DT_BUILTIN(MPI_DOUBLE);
    MPIC_INIT_DT_BUILTIN(MPI_LONG_DOUBLE);
    MPIC_INIT_DT_BUILTIN(MPI_LONG_LONG);
    MPIC_INIT_DT_BUILTIN(MPI_UNSIGNED_LONG_LONG);
    MPIC_INIT_DT_BUILTIN(MPI_PACKED);
    MPIC_INIT_DT_BUILTIN(MPI_LB);
    MPIC_INIT_DT_BUILTIN(MPI_UB);
    MPIC_INIT_DT_BUILTIN(MPI_2INT);

    /* C99 types */
    MPIC_INIT_DT_BUILTIN(MPI_INT8_T);      /* 20 */
    MPIC_INIT_DT_BUILTIN(MPI_INT16_T);
    MPIC_INIT_DT_BUILTIN(MPI_INT32_T);
    MPIC_INIT_DT_BUILTIN(MPI_INT64_T);
    MPIC_INIT_DT_BUILTIN(MPI_UINT8_T);
    MPIC_INIT_DT_BUILTIN(MPI_UINT16_T);
    MPIC_INIT_DT_BUILTIN(MPI_UINT32_T);
    MPIC_INIT_DT_BUILTIN(MPI_UINT64_T);
    MPIC_INIT_DT_BUILTIN(MPI_C_BOOL);
    MPIC_INIT_DT_BUILTIN(MPI_C_FLOAT_COMPLEX);
    MPIC_INIT_DT_BUILTIN(MPI_C_DOUBLE_COMPLEX);    /* 30 */
    MPIC_INIT_DT_BUILTIN(MPI_C_LONG_DOUBLE_COMPLEX);
    /* address/offset/count types */
    MPIC_INIT_DT_BUILTIN(MPI_AINT);
    MPIC_INIT_DT_BUILTIN(MPI_OFFSET);
    MPIC_INIT_DT_BUILTIN(MPI_COUNT);
#ifdef HAVE_FORTRAN_BINDING
    MPIC_INIT_DT_BUILTIN(MPI_COMPLEX);
    MPIC_INIT_DT_BUILTIN(MPI_DOUBLE_COMPLEX);
    MPIC_INIT_DT_BUILTIN(MPI_LOGICAL);
    MPIC_INIT_DT_BUILTIN(MPI_REAL);
    MPIC_INIT_DT_BUILTIN(MPI_DOUBLE_PRECISION);
    MPIC_INIT_DT_BUILTIN(MPI_INTEGER);     /* 40 */
    MPIC_INIT_DT_BUILTIN(MPI_2INTEGER);
#ifdef MPICH_DEFINE_2COMPLEX
    MPIC_INIT_DT_BUILTIN(MPI_2COMPLEX);
    MPIC_INIT_DT_BUILTIN(MPI_2DOUBLE_COMPLEX);
#endif
    MPIC_INIT_DT_BUILTIN(MPI_2REAL);
    MPIC_INIT_DT_BUILTIN(MPI_2DOUBLE_PRECISION);
    MPIC_INIT_DT_BUILTIN(MPI_CHARACTER);
    MPIC_INIT_DT_BUILTIN(MPI_REAL4);
    MPIC_INIT_DT_BUILTIN(MPI_REAL8);
    MPIC_INIT_DT_BUILTIN(MPI_REAL16);
    MPIC_INIT_DT_BUILTIN(MPI_COMPLEX8);    /* 50 */
    MPIC_INIT_DT_BUILTIN(MPI_COMPLEX16);
    MPIC_INIT_DT_BUILTIN(MPI_COMPLEX32);
    MPIC_INIT_DT_BUILTIN(MPI_INTEGER1);
    MPIC_INIT_DT_BUILTIN(MPI_INTEGER2);
    MPIC_INIT_DT_BUILTIN(MPI_INTEGER4);
    MPIC_INIT_DT_BUILTIN(MPI_INTEGER8);

    if (MPI_INTEGER16 != MPI_DATATYPE_NULL) {
        MPIC_INIT_DT_BUILTIN(MPI_INTEGER16);
    }
#endif

    MPIC_INIT_DT_BUILTIN(MPI_FLOAT_INT);
    MPIC_INIT_DT_BUILTIN(MPI_DOUBLE_INT);
    MPIC_INIT_DT_BUILTIN(MPI_LONG_INT);
    MPIC_INIT_DT_BUILTIN(MPI_SHORT_INT);   /* 60 */
    MPIC_INIT_DT_BUILTIN(MPI_LONG_DOUBLE_INT);

#ifdef HAVE_CXX_BINDING
    MPIC_INIT_CXX_DT_BUILTIN(MPI_CXX_BOOL);
    MPIC_INIT_CXX_DT_BUILTIN(MPI_CXX_FLOAT_COMPLEX);
    MPIC_INIT_CXX_DT_BUILTIN(MPI_CXX_DOUBLE_COMPLEX);
    MPIC_INIT_CXX_DT_BUILTIN(MPI_CXX_LONG_DOUBLE_COMPLEX);
#endif
    return mpi_errno;

}

MPL_STATIC_INLINE_PREFIX int MPIC_init_builtin_ops()
{
    int mpi_errno = MPI_SUCCESS;
    int i;
    for (i = 0; i < MPIR_OP_N_BUILTIN; ++i) {
        MPIR_Op_builtin[i].handle = MPI_MAX + i;
        MPIC_op_init(&MPIR_Op_builtin[i]);
    }
    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX int MPIC_Op_get_ptr(MPI_Op op, MPIR_Op ** ptr)
{
    int mpi_errno = MPI_SUCCESS;
    if(HANDLE_GET_KIND(op) == HANDLE_KIND_BUILTIN)
    {
        *ptr = &MPIR_Op_builtin[(op & 0xFF) -1];
    } else {
        MPIR_Op *tmp;
        MPIR_Op_get_ptr(op, tmp);
        *ptr = tmp;
    }
    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX int MPIC_init()
{
    int mpi_errno = MPI_SUCCESS;
    MPIC_progress_global.progress_fn = MPID_Progress_test;
    MPIC_comm_counter=0;
    MPIC_sched_table = NULL;
    MPIC_STUB_init();
    MPIC_MPICH_init();

    MPIC_STUB_STUB_init();
    MPIC_STUB_KARY_init();
    MPIC_STUB_KNOMIAL_init();
    MPIC_STUB_RECEXCH_init();
    MPIC_STUB_DISSEM_init();
    MPIC_MPICH_STUB_init();
    MPIC_MPICH_KARY_init();
    MPIC_MPICH_KNOMIAL_init();
    MPIC_MPICH_RECEXCH_init();
    MPIC_MPICH_DISSEM_init();
    MPIC_BMPICH_KARY_init();
    MPIC_BMPICH_KNOMIAL_init();
    MPIC_X_TREEBASIC_init();

#ifdef MPIDI_NM_COLL_INIT_HOOK
    MPIDI_NM_COLL_INIT_HOOK;
#endif
    MPIC_init_builtin_dt();
    MPIC_init_builtin_ops();
    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX int MPIC_finalize()
{
    int mpi_errno = MPI_SUCCESS;
    MPIC_delete_sched_table();
    return mpi_errno;
}
#undef MPIC_INIT_DT_BUILTIN
#undef MPIC_INIT_CXX_DT_BUILTIN
#endif /* MPIC_COLL_IMPL_H_INCLUDED */
