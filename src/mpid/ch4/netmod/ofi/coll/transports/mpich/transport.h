/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2012 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef MPICHTRANSPORT_H_INCLUDED
#define MPICHTRANSPORT_H_INCLUDED

#include "../../include/tsp_namespace_pre.h"

#ifndef TSP_NAMESPACE
#error "TSP_NAMESPACE must be defined before including a collective transport"
#endif

static inline int TSP_init()
{
    return 0;
}

static inline void TSP_sched_init(TSP_sched_t *sched)
{
    sched->cntr       = 0;
    sched->cur_thresh = 0;
    sched->posted     = 0;
    sched->completed  = 0;
}

static inline void TSP_sched_commit(TSP_sched_t *sched)
{

}

static inline void TSP_sched_start(TSP_sched_t *sched)
{

}

static inline void TSP_sched_finalize(TSP_sched_t *sched)
{

}


static inline int TSP_init_control_dt(TSP_dt_t *dt)
{
    dt->mpi_dt = MPI_CHAR;
    return 0;
}

static inline void TSP_opinfo(TSP_op_t  *op,
                              int       *is_commutative)
{
    MPIR_Op *op_ptr;

    if(HANDLE_GET_KIND(op->mpi_op) == HANDLE_KIND_BUILTIN)
        *is_commutative=1;
    else {
        MPIR_Op_get_ptr(op->mpi_op, op_ptr);

        if(op_ptr->kind == MPIR_OP_KIND__USER_NONCOMMUTE)
            *is_commutative = 0;
        else
            *is_commutative = 1 ;
    }
}

static inline int TSP_isinplace(const void *buf)
{
    if(buf == MPI_IN_PLACE) return 1;

    return 0;
}
static inline void TSP_dtinfo(TSP_dt_t *dt,
                              int      *iscontig,
                              size_t   *size,
                              size_t   *out_extent,
                              size_t   *lower_bound)
{
    MPI_Aint true_lb,extent,true_extent,type_size;
    MPID_Datatype_get_size_macro(dt->mpi_dt, type_size);
    MPIR_Type_get_true_extent_impl(dt->mpi_dt,&true_lb,&true_extent);
    MPID_Datatype_get_extent_macro(dt->mpi_dt,extent);
    MPIDI_Datatype_check_contig(dt->mpi_dt,*iscontig);
    *size = type_size;
    *out_extent  = MPL_MAX(extent,true_extent);
    *lower_bound = true_lb;
}

static inline void TSP_addref_dt(TSP_dt_t *dt,
                                 int       up)
{
    MPIR_Datatype *dt_ptr;

    if(HANDLE_GET_KIND(dt->mpi_dt) != HANDLE_KIND_BUILTIN) {
        MPID_Datatype_get_ptr(dt->mpi_dt, dt_ptr);

        if(up)
            MPID_Datatype_add_ref(dt_ptr);
        else
            MPID_Datatype_release(dt_ptr);
    }
}

static inline void TSP_fence(TSP_sched_t *sched)
{
    sched->cur_thresh+=sched->posted;

    if(0) fprintf(stderr, "TSP(mpich) : sched [fence] posted=%ld cur_thresh=%ld\n",
                      sched->posted, sched->cur_thresh);

    sched->posted = 0;
}

static inline void TSP_addref_op(TSP_op_t *op,
                                 int       up)
{
    MPIR_Op     *op_ptr;

    if(HANDLE_GET_KIND(op->mpi_op) != HANDLE_KIND_BUILTIN) {
        MPIR_Op_get_ptr(op->mpi_op, op_ptr);

        if(up)
            MPIR_Op_add_ref(op_ptr);
        else
            MPIR_Op_release(op_ptr);
    }
}

static inline void TSP_addref_dt_nb(TSP_dt_t    *dt,
                                    int          up,
                                    TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                = TSP_KIND_ADDREF_DT;
    req->state               = TSP_STATE_INIT;
    req->completion_cntr     = &sched->cntr;
    req->trigger_cntr        = &sched->cntr;
    req->threshold           = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.addref_dt.dt = dt;
    req->nbargs.addref_dt.up = up;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [addref_dt]\n",sched->total);
}

static inline void TSP_addref_op_nb(TSP_op_t    *op,
                                    int          up,
                                    TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                = TSP_KIND_ADDREF_OP;
    req->state               = TSP_STATE_INIT;
    req->completion_cntr     = &sched->cntr;
    req->trigger_cntr        = &sched->cntr;
    req->threshold           = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.addref_op.op = op;
    req->nbargs.addref_op.up = up;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [addref_op]\n",sched->total);
}
static inline int TSP_test(TSP_sched_t *sched);
static inline void TSP_send(const void  *buf,
                            int          count,
                            TSP_dt_t    *dt,
                            int          dest,
                            int          tag,
                            TSP_comm_t  *comm,
                            TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                  = TSP_KIND_SEND;
    req->state                 = TSP_STATE_INIT;
    req->completion_cntr       = &sched->cntr;
    req->trigger_cntr          = &sched->cntr;
    req->threshold             = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.sendrecv.buf   = (void *)buf;
    req->nbargs.sendrecv.count = count;
    req->nbargs.sendrecv.dt    = dt;
    req->nbargs.sendrecv.dest  = dest;
    req->nbargs.sendrecv.tag   = tag;
    req->nbargs.sendrecv.comm  = comm;
    req->mpid_req[1]           = NULL;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [send]\n",sched->total);

    TSP_test(sched);
}

static inline void TSP_send_accumulate(const void  *buf,
                                       int          count,
                                       TSP_dt_t    *dt,
                                       TSP_op_t    *op,
                                       int          dest,
                                       int          tag,
                                       TSP_comm_t  *comm,
                                       TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                  = TSP_KIND_SEND;
    req->state                 = TSP_STATE_INIT;
    req->completion_cntr       = &sched->cntr;
    req->trigger_cntr          = &sched->cntr;
    req->threshold             = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.sendrecv.buf   = (void *)buf;
    req->nbargs.sendrecv.count = count;
    req->nbargs.sendrecv.dt    = dt;
    req->nbargs.sendrecv.dest  = dest;
    req->nbargs.sendrecv.tag   = tag;
    req->nbargs.sendrecv.comm  = comm;
    req->mpid_req[1]           = NULL;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [send_accumulate]\n",sched->total);

    TSP_test(sched);
}

static inline void TSP_recv(void        *buf,
                            int          count,
                            TSP_dt_t    *dt,
                            int          source,
                            int          tag,
                            TSP_comm_t  *comm,
                            TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                  = TSP_KIND_RECV;
    req->state                 = TSP_STATE_INIT;
    req->completion_cntr       = &sched->cntr;
    req->trigger_cntr          = &sched->cntr;
    req->threshold             = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.sendrecv.buf   = buf;
    req->nbargs.sendrecv.count = count;
    req->nbargs.sendrecv.dt    = dt;
    req->nbargs.sendrecv.dest  = source;
    req->nbargs.sendrecv.tag   = tag;
    req->nbargs.sendrecv.comm  = comm;
    req->mpid_req[1]           = NULL;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [recv]\n",sched->total);

    TSP_test(sched);
}

static inline int TSP_NAMESPACE(queryfcn)(
    void       *data,
    MPI_Status *status)
{
    TSP_recv_reduce_arg_t *rr
        = (TSP_recv_reduce_arg_t *)data;

    if(rr->req->mpid_req[0] == NULL && !rr->done) {
        MPI_Datatype dt = rr->datatype->mpi_dt;
        MPI_Op       op = rr->op->mpi_op;

        if(rr->flags==0 || rr->flags & TSP_FLAG_REDUCE_L) {
            if(0) fprintf(stderr, "  --> MPICH transport (recv_reduce L) complete to %p\n",
                              rr->inoutbuf);

            MPIR_Reduce_local_impl(rr->inbuf,rr->inoutbuf,rr->count,dt,op);
        } else {
            if(0) fprintf(stderr, "  --> MPICH transport (recv_reduce R) complete to %p\n",
                              rr->inoutbuf);

            MPIR_Reduce_local_impl(rr->inoutbuf,rr->inbuf,rr->count,dt,op);
            MPIR_Localcopy(rr->inbuf,rr->count,dt,
                           rr->inoutbuf,rr->count,dt);
        }

        MPL_free(rr->inbuf);
        MPIR_Grequest_complete_impl(rr->req->mpid_req[1]);
        rr->done = 1;
    }

    status->MPI_SOURCE = MPI_UNDEFINED;
    status->MPI_TAG    = MPI_UNDEFINED;
    MPI_Status_set_cancelled(status, 0);
    MPI_Status_set_elements(status, MPI_BYTE, 0);
    return MPI_SUCCESS;
}

static inline void TSP_recv_reduce(void        *buf,
                                   int          count,
                                   TSP_dt_t    *datatype,
                                   TSP_op_t    *op,
                                   int          source,
                                   int          tag,
                                   TSP_comm_t  *comm,
                                   TSP_sched_t *sched,
                                   uint64_t     flags)
{
    int        iscontig;
    size_t     type_size, out_extent, lower_bound;
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                = TSP_KIND_RECV_REDUCE;
    req->state               = TSP_STATE_INIT;
    req->completion_cntr     = &sched->cntr;
    req->trigger_cntr        = &sched->cntr;
    req->threshold           = sched->cur_thresh;
    sched->posted++;
    sched->total++;

    TSP_dtinfo(datatype,&iscontig,&type_size,&out_extent,&lower_bound);
    req->nbargs.recv_reduce.inbuf    = MPL_malloc(count*out_extent);
    req->nbargs.recv_reduce.inoutbuf = buf;
    req->nbargs.recv_reduce.count    = count;
    req->nbargs.recv_reduce.datatype = datatype;
    req->nbargs.recv_reduce.op       = op;
    req->nbargs.recv_reduce.source   = source;
    req->nbargs.recv_reduce.tag      = tag;
    req->nbargs.recv_reduce.comm     = comm;
    req->nbargs.recv_reduce.req      = req;
    req->nbargs.recv_reduce.done     = 0;
    req->nbargs.recv_reduce.flags    = flags;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [recv_reduce]\n",sched->total);

    TSP_test(sched);
}

static inline int TSP_rank(TSP_comm_t  *comm)
{
    return comm->mpid_comm->rank;
}

static inline int TSP_size(TSP_comm_t *comm)
{
    return comm->mpid_comm->local_size;
}

static inline void TSP_reduce_local(const void  *inbuf,
                                    void        *inoutbuf,
                                    int          count,
                                    TSP_dt_t    *datatype,
                                    TSP_op_t    *operation,
                                    TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                         = TSP_KIND_REDUCE_LOCAL;
    req->state                        = TSP_STATE_INIT;
    req->completion_cntr              = &sched->cntr;
    req->trigger_cntr                 = &sched->cntr;
    req->threshold                    = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.reduce_local.inbuf    = inbuf;
    req->nbargs.reduce_local.inoutbuf = inoutbuf;
    req->nbargs.reduce_local.count    = count;
    req->nbargs.reduce_local.dt       = datatype;
    req->nbargs.reduce_local.op       = operation;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [reduce_local]\n",sched->total);
}

static inline int TSP_dtcopy(void       *tobuf,
                             int         tocount,
                             TSP_dt_t   *totype,
                             const void *frombuf,
                             int         fromcount,
                             TSP_dt_t   *fromtype)
{
    return MPIR_Localcopy(frombuf,   /* yes, parameters are reversed        */
                          fromcount, /* MPICH forgot what memcpy looks like */
                          fromtype->mpi_dt,
                          tobuf,
                          tocount,
                          totype->mpi_dt);
}

static inline void TSP_dtcopy_nb(void        *tobuf,
                                 int          tocount,
                                 TSP_dt_t    *totype,
                                 const void  *frombuf,
                                 int          fromcount,
                                 TSP_dt_t    *fromtype,
                                 TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind   = TSP_KIND_DTCOPY;
    req->state  = TSP_STATE_INIT;

    req->completion_cntr          = &sched->cntr;
    req->trigger_cntr             = &sched->cntr;
    req->threshold                = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.dtcopy.tobuf      = tobuf;
    req->nbargs.dtcopy.tocount    = tocount;
    req->nbargs.dtcopy.totype     = totype;
    req->nbargs.dtcopy.frombuf    = frombuf;
    req->nbargs.dtcopy.fromcount  = fromcount;
    req->nbargs.dtcopy.fromtype   = fromtype;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [dt_copy]\n",sched->total);
}

static inline void *TSP_allocate_mem(size_t size)
{
    return MPL_malloc(size);
}

static inline void TSP_free_mem(void *ptr)
{
    MPL_free(ptr);
}

static inline void TSP_free_mem_nb(void        *ptr,
                                   TSP_sched_t *sched)
{
    TSP_req_t *req;
    MPIR_Assert(sched->total < 32);
    req = &sched->requests[sched->total];
    req->kind                = TSP_KIND_FREE_MEM;
    req->state               = TSP_STATE_INIT;
    req->completion_cntr     = &sched->cntr;
    req->trigger_cntr        = &sched->cntr;
    req->threshold           = sched->cur_thresh;
    sched->posted++;
    sched->total++;
    req->nbargs.free_mem.ptr = ptr;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [free_mem]\n",sched->total);
}


static inline int TSP_test(TSP_sched_t *sched)
{
    TSP_req_t *req, *rp;
    int i;
    req = &sched->requests[0];
    /* First check for ops that need to be issued */
    int consecutiveComplete=1;

    for(i=sched->completed; i<sched->total; i++) {
        rp = &req[i];

        if(rp->state == TSP_STATE_COMPLETE && consecutiveComplete) {
            sched->completed++;
            continue;
        } else
            consecutiveComplete=0;

        if(rp->state == TSP_STATE_INIT && *(rp->trigger_cntr) >= rp->threshold) {
            switch(rp->kind) {
                case TSP_KIND_SEND: {
                    if(0) fprintf(stderr, "  --> MPICH transport (isend) issued\n");

                    MPIR_Errflag_t  errflag  = MPIR_ERR_NONE;
                    MPIC_Isend(rp->nbargs.sendrecv.buf,
                               rp->nbargs.sendrecv.count,
                               rp->nbargs.sendrecv.dt->mpi_dt,
                               rp->nbargs.sendrecv.dest,
                               rp->nbargs.sendrecv.tag,
                               rp->nbargs.sendrecv.comm->mpid_comm,
                               &rp->mpid_req[0],
                               &errflag);
                    rp->state = TSP_STATE_ISSUED;
                }
                break;

                case TSP_KIND_RECV: {
                    if(0) fprintf(stderr, "  --> MPICH transport (irecv) issued\n");

                    MPIC_Irecv(rp->nbargs.sendrecv.buf,
                               rp->nbargs.sendrecv.count,
                               rp->nbargs.sendrecv.dt->mpi_dt,
                               rp->nbargs.sendrecv.dest,
                               rp->nbargs.sendrecv.tag,
                               rp->nbargs.sendrecv.comm->mpid_comm,
                               &rp->mpid_req[0]);
                    rp->state = TSP_STATE_ISSUED;
                }
                break;

                case TSP_KIND_ADDREF_DT:
                    TSP_addref_dt(rp->nbargs.addref_dt.dt,
                                  rp->nbargs.addref_dt.up);

                    if(0) fprintf(stderr, "  --> MPICH transport (addref dt) complete\n");

                    rp->state = TSP_STATE_COMPLETE;
                    (*rp->completion_cntr)++;
                    break;

                case TSP_KIND_ADDREF_OP:
                    TSP_addref_op(rp->nbargs.addref_op.op,
                                  rp->nbargs.addref_op.up);

                    if(0) fprintf(stderr, "  --> MPICH transport (addref op) complete\n");

                    rp->state = TSP_STATE_COMPLETE;
                    (*rp->completion_cntr)++;
                    break;

                case TSP_KIND_DTCOPY:
                    TSP_dtcopy(rp->nbargs.dtcopy.tobuf,
                               rp->nbargs.dtcopy.tocount,
                               rp->nbargs.dtcopy.totype,
                               rp->nbargs.dtcopy.frombuf,
                               rp->nbargs.dtcopy.fromcount,
                               rp->nbargs.dtcopy.fromtype);

                    if(0) fprintf(stderr, "  --> MPICH transport (dtcopy) complete\n");

                    rp->state = TSP_STATE_COMPLETE;
                    (*rp->completion_cntr)++;
                    break;

                case TSP_KIND_FREE_MEM:
                    if(0) fprintf(stderr, "  --> MPICH transport (freemem) complete\n");

                    TSP_free_mem(rp->nbargs.free_mem.ptr);
                    rp->state = TSP_STATE_COMPLETE;
                    (*rp->completion_cntr)++;
                    break;

                case TSP_KIND_RECV_REDUCE: {
                    MPIC_Irecv(rp->nbargs.recv_reduce.inbuf,
                               rp->nbargs.recv_reduce.count,
                               rp->nbargs.recv_reduce.datatype->mpi_dt,
                               rp->nbargs.recv_reduce.source,
                               rp->nbargs.recv_reduce.tag,
                               rp->nbargs.recv_reduce.comm->mpid_comm,
                               &rp->mpid_req[0]);
                    MPIR_Grequest_start_impl(TSP_NAMESPACE(queryfcn),
                                             NULL,
                                             NULL,
                                             &rp->nbargs.recv_reduce,
                                             &rp->mpid_req[1]);

                    if(0) fprintf(stderr, "  --> MPICH transport (recv_reduce) issued\n");

                    rp->state = TSP_STATE_ISSUED;
                }
                break;

                case TSP_KIND_REDUCE_LOCAL:
                    MPIR_Reduce_local_impl(rp->nbargs.reduce_local.inbuf,
                                           rp->nbargs.reduce_local.inoutbuf,
                                           rp->nbargs.reduce_local.count,
                                           rp->nbargs.reduce_local.dt->mpi_dt,
                                           rp->nbargs.reduce_local.op->mpi_op);

                    if(0) fprintf(stderr, "  --> MPICH transport (reduce local) complete\n");

                    rp->state = TSP_STATE_COMPLETE;
                    (*rp->completion_cntr)++;
                    break;
            }
        }
    }

    /* Now check for issued ops that have been completed */
    for(i=sched->completed; i<sched->total; i++) {
        rp = &req[i];

        if(rp->state == TSP_STATE_ISSUED) {
            MPI_Status     status;
            MPIR_Request  *mpid_req0 = rp->mpid_req[0];
            MPIR_Request  *mpid_req1 = rp->mpid_req[1];

            if(mpid_req1) {
                (mpid_req1->u.ureq.greq_fns->query_fn)
                (mpid_req1->u.ureq.greq_fns->grequest_extra_state,
                 &status);
            }

            switch(rp->kind) {
                case TSP_KIND_SEND:
                case TSP_KIND_RECV:
                case TSP_KIND_RECV_REDUCE:
                    if(mpid_req0 && MPIR_Request_is_complete(mpid_req0)) {
                        MPIR_Request_free(mpid_req0);
                        rp->mpid_req[0] = NULL;
                    }

                    if(mpid_req1 && MPIR_Request_is_complete(mpid_req1)) {
                        MPIR_Request_free(mpid_req1);
                        rp->mpid_req[1] = NULL;
                    }

                    if(!rp->mpid_req[0] && !rp->mpid_req[1]) {
                        if(0) {
                            fprintf(stderr, "  --> MPICH transport (kind=%d) complete\n",
                                    rp->kind);
                            fprintf(stderr, "data send/recvd: %d\n", *(int *)(rp->nbargs.sendrecv.buf));
                        }

                        rp->state = TSP_STATE_COMPLETE;
                        (*rp->completion_cntr)++;
                    }

                    break;

                default:
                    break;
            }
        }
    }

    if(sched->cntr >= sched->cur_thresh) {
        if(0) fprintf(stderr, "  --> MPICH transport (test) complete:  cntr->value=%ld threshold=%ld\n",
                          sched->cntr, sched->cur_thresh);

        return 1;
    } else
        return 0;
}

#endif
