/*
 *
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

static inline void *MPIC_MPICH_allocate_mem(size_t size)
{
    return MPL_malloc(size);
}

static inline void MPIC_MPICH_free_mem(void *ptr)
{
    MPL_free(ptr);
}

static inline void MPIC_MPICH_issue_vtx(int,MPIC_MPICH_vtx_t*,MPIC_MPICH_sched_t*);

static inline void MPIC_MPICH_decrement_num_unfinished_dependecies(MPIC_MPICH_vtx_t *vtxp, MPIC_MPICH_sched_t *sched){
    /*for each outgoing vertex of vertex vtxp, decrement number of unfinished dependencies*/
    MPIC_MPICH_int_array *outvtcs = &vtxp->outvtcs;
    if(0) fprintf(stderr,"num outvtcs of %d = %d\n", vtxp->id, outvtcs->used);
    int i;
    for(i=0; i<outvtcs->used; i++){
        int num_unfinished_dependencies = --(sched->vtcs[outvtcs->array[i]].num_unfinished_dependencies);
        /*if all dependencies are complete, issue the op*/
        if(num_unfinished_dependencies == 0){
            if(0) fprintf(stderr,"issuing request number: %d\n",outvtcs->array[i]);
            MPIC_MPICH_issue_vtx(outvtcs->array[i], &sched->vtcs[outvtcs->array[i]], sched);
        }
    }
}

static inline void MPIC_MPICH_record_vtx_completion(MPIC_MPICH_vtx_t *vtxp, MPIC_MPICH_sched_t *sched){
    vtxp->state = MPIC_MPICH_STATE_COMPLETE;
    /*update the dependencies*/
    sched->num_completed++;
    if(0) fprintf(stderr, "num completed vertices: %d\n", sched->num_completed);
    MPIC_MPICH_decrement_num_unfinished_dependecies(vtxp,sched);
}

static inline void MPIC_MPICH_record_vtx_issue(MPIC_MPICH_vtx_t *vtxp, MPIC_MPICH_sched_t *sched){
    vtxp->state = MPIC_MPICH_STATE_ISSUED;
    
    if(sched->issued_head == NULL){
        sched->issued_head = sched->last_issued = vtxp;
        sched->last_issued->next_issued = sched->vtx_iter;
    }else if(sched->last_issued->next_issued==vtxp){
        sched->last_issued = vtxp;
    }
    else{
        sched->last_issued->next_issued = vtxp;
        vtxp->next_issued = sched->vtx_iter;
        sched->last_issued = vtxp;
    }

    /*print issued task list*/
    if(0){
        MPIC_MPICH_vtx_t *vtx = sched->issued_head;
        fprintf(stderr, "issued list: ");
        while(vtx){
            fprintf(stderr, "%d ", vtx->id);
            vtx = vtx->next_issued;
        }
        if(sched->vtx_iter)
            fprintf(stderr, ", vtx_iter: %d\n", sched->vtx_iter->id);
        else
            fprintf(stderr,"\n");
    }
}

static inline int MPIC_MPICH_init()
{
    MPIC_global_instance.tsp_mpich.control_dt = MPI_BYTE;
    return 0;
}

static inline int MPIC_MPICH_comm_init(MPIC_MPICH_comm_t * comm, void* base)
{
    MPIR_Comm * mpi_comm = container_of(base, MPIR_Comm, ch4_coll);
    comm->mpid_comm = mpi_comm;
}

static inline int MPIC_MPICH_comm_cleanup(MPIC_MPICH_comm_t * comm)
{
        comm->mpid_comm = NULL;
}

static inline void MPIC_MPICH_reset_issued_list(MPIC_MPICH_sched_t *sched){
    sched->issued_head = NULL;
    int i;
    /*If schedule is being reused, reset only used vtcs
    else it is a new schedule being generated and therefore
    reset all the vtcs**/
    int nvtcs = (sched->total==0)?sched->max_vtcs:sched->total;
    
    for(i=0; i<nvtcs; i++){
        sched->vtcs[i].next_issued = NULL;
    }
    sched->vtx_iter = NULL;
}

static inline void MPIC_MPICH_vtx_init(MPIC_MPICH_vtx_t* vtx, int max_edges){
}

static inline void MPIC_MPICH_sched_init(MPIC_MPICH_sched_t *sched, int tag)
{
    sched->total      = 0;
    sched->num_completed  = 0;
    sched->last_wait = -1;
    sched->tag       = tag;
    sched->max_vtcs = MPIC_MPICH_MAX_REQUESTS;
    sched->max_edges_per_vtx = MPIC_MPICH_MAX_EDGES;
    

    /*allocate memory for vertices*/
    sched->vtcs = MPIC_MPICH_allocate_mem(sizeof(MPIC_MPICH_vtx_t)*sched->max_vtcs);
    /*initialize array for storing memory buffer addresses*/
    sched->buf_array.array = (void**)MPIC_MPICH_allocate_mem(sizeof(void*)*sched->max_vtcs);
    sched->buf_array.size = sched->max_vtcs;
    sched->buf_array.used = 0;

    MPIC_MPICH_reset_issued_list(sched);
}

static inline void MPIC_MPICH_sched_reset(MPIC_MPICH_sched_t *sched, int tag){
    int i;
    sched->num_completed=0;
    sched->tag = tag;
    for(i=0; i<sched->total; i++){
        MPIC_MPICH_vtx_t* vtx = &sched->vtcs[i];
        vtx->state = MPIC_MPICH_STATE_INIT;
        vtx->num_unfinished_dependencies = vtx->invtcs.used;
        if(vtx->kind==MPIC_MPICH_KIND_RECV_REDUCE)
            vtx->nbargs.recv_reduce.done=0;
    }
    MPIC_MPICH_reset_issued_list(sched);
}

static inline void MPIC_MPICH_sched_commit(MPIC_MPICH_sched_t *sched)
{

}

static inline void MPIC_MPICH_sched_start(MPIC_MPICH_sched_t *sched)
{

}

static inline void MPIC_MPICH_sched_finalize(MPIC_MPICH_sched_t *sched)
{

}


static inline void MPIC_MPICH_opinfo(MPIC_MPICH_op_t  op,
                              int       *is_commutative)
{
    MPIR_Op *op_ptr;

    if(HANDLE_GET_KIND(op) == HANDLE_KIND_BUILTIN)
        *is_commutative=1;
    else {
        MPIR_Op_get_ptr(op, op_ptr);

        if(op_ptr->kind == MPIR_OP_KIND__USER_NONCOMMUTE)
            *is_commutative = 0;
        else
            *is_commutative = 1 ;
    }
}

static inline int MPIC_MPICH_isinplace(const void *buf)
{
    if(buf == MPI_IN_PLACE) return 1;

    return 0;
}
static inline void MPIC_MPICH_dtinfo(MPIC_MPICH_dt_t dt,
                              int      *iscontig,
                              size_t   *size,
                              size_t   *out_extent,
                              size_t   *lower_bound)
{
    MPI_Aint true_lb,extent,true_extent,type_size;
    MPID_Datatype_get_size_macro(dt, type_size);
    MPIR_Type_get_true_extent_impl(dt,&true_lb,&true_extent);
    MPID_Datatype_get_extent_macro(dt,extent);
    MPIDI_Datatype_check_contig(dt,*iscontig);
    *size = type_size;
    *out_extent  = MPL_MAX(extent,true_extent);
    *lower_bound = true_lb;
}

static inline void MPIC_MPICH_addref_dt(MPIC_MPICH_dt_t dt,
                                 int       up)
{
    MPIR_Datatype *dt_ptr;

    if(HANDLE_GET_KIND(dt) != HANDLE_KIND_BUILTIN) {
        MPID_Datatype_get_ptr(dt, dt_ptr);

        if(up)
            MPID_Datatype_add_ref(dt_ptr);
        else
            MPID_Datatype_release(dt_ptr);
    }
}

static inline void MPIC_MPICH_allocate_vtx(MPIC_MPICH_vtx_t *vtx, int invtcs_array_size, int outvtcs_array_size){
    vtx->invtcs.array = MPIC_MPICH_allocate_mem(sizeof(int)*invtcs_array_size);
    vtx->invtcs.size = invtcs_array_size;
    vtx->outvtcs.array = MPIC_MPICH_allocate_mem(sizeof(int)*outvtcs_array_size);
    vtx->outvtcs.size = outvtcs_array_size;
}

static inline void MPIC_MPICH_init_vtx(MPIC_MPICH_sched_t *sched, MPIC_MPICH_vtx_t *vtx, int id){
    /*allocate array for storing incoming and outgoing edges*/
    MPIC_MPICH_allocate_vtx(vtx, sched->max_edges_per_vtx, sched->max_edges_per_vtx);
    vtx->invtcs.used                 = 0;
    vtx->outvtcs.used                = 0;

    vtx->state                       = MPIC_MPICH_STATE_INIT;
    vtx->id                          = id;   
    vtx->num_unfinished_dependencies = 0;
}

static inline void MPIC_MPICH_free_vtx(MPIC_MPICH_vtx_t *vtx){
    MPIC_MPICH_free_mem(vtx->invtcs.array);
    MPIC_MPICH_free_mem(vtx->outvtcs.array);
}

static inline void MPIC_MPICH_add_elems_int_array(MPIC_MPICH_int_array* in, int n_elems, int *elems){
    if(in->size+n_elems > in->size){
        int old_size = in->size;
        int* old_array = in->array;
        in->size *= 2;
        /*reallocate array*/
        in->array = MPIC_MPICH_allocate_mem(sizeof(int)*in->size);
        /*copy old elements*/
        memcpy(in->array, old_array, sizeof(int)*in->used);
        /*free old array*/
        MPIC_MPICH_free_mem(old_array);
    }
    /*add new elements*/
    memcpy(in->array+in->used, elems, sizeof(int)*n_elems);
    in->used += n_elems;
}

/*this vertex sets the incoming edges (inedges) to vtx
and also add vtx to the outgoing edge list of vertices in inedges*/
static inline void MPIC_MPICH_add_vtx_dependencies(MPIC_MPICH_sched_t *sched, int vtx_id, int n_invtcs, int *invtcs){
    int i;
    MPIC_MPICH_vtx_t *vtx = &sched->vtcs[vtx_id];
    MPIC_MPICH_int_array *in = &vtx->invtcs;
    if(0) fprintf(stderr,"updating invtcs of vtx %d, kind %d, in->used %d, n_invtcs %d\n",vtx_id,vtx->kind,in->used, n_invtcs);
    /*insert the incoming edges*/
    MPIC_MPICH_add_elems_int_array(in, n_invtcs, invtcs);

    MPIC_MPICH_int_array *outvtcs;
    /*update the outgoing edges of incoming vertices*/
    for(i=0; i<n_invtcs; i++){
        if(0)fprintf(stderr,"invtx: %d\n", invtcs[i]);
        outvtcs = &sched->vtcs[invtcs[i]].outvtcs;
        MPIC_MPICH_add_elems_int_array(outvtcs, 1, &vtx_id);

        if(sched->vtcs[invtcs[i]].state != MPIC_MPICH_STATE_COMPLETE)
            vtx->num_unfinished_dependencies++;
    }

    /*check if there was any MPIC_MPICH_wait operation and add appropriate dependencies*/
    if(sched->last_wait != -1 && sched->last_wait != vtx_id){
        /*add incoming edge to vtx*/
        MPIC_MPICH_add_elems_int_array(in, 1, &(sched->last_wait));

        /*add vtx as outgoing vtx of last_wait*/
        outvtcs = &sched->vtcs[sched->last_wait].outvtcs;
        MPIC_MPICH_add_elems_int_array(outvtcs, 1, &vtx_id);

        if(sched->vtcs[sched->last_wait].state != MPIC_MPICH_STATE_COMPLETE)
            vtx->num_unfinished_dependencies++;
    }
}

static inline int MPIC_MPICH_get_new_vtx(MPIC_MPICH_sched_t *sched, MPIC_MPICH_vtx_t **vtx){
    if(sched->total == sched->max_vtcs){/*increase array size*/
        int old_size = sched->max_vtcs;
        MPIC_MPICH_vtx_t* old_vtcs = sched->vtcs;
        sched->max_vtcs *= 2;
        sched->vtcs = MPIC_MPICH_allocate_mem(sizeof(MPIC_MPICH_vtx_t)*sched->max_vtcs);
        memcpy(sched->vtcs, old_vtcs, sizeof(MPIC_MPICH_vtx_t)*old_size);
        int i;
        for(i=0; i<old_size; i++){
            MPIC_MPICH_allocate_vtx(&sched->vtcs[i], old_vtcs[i].invtcs.size, old_vtcs[i].outvtcs.size);
            /*used values get copied in the memcpy above*/
            memcpy(sched->vtcs[i].invtcs.array, old_vtcs[i].invtcs.array, sizeof(int)*old_vtcs[i].invtcs.used);
            memcpy(sched->vtcs[i].outvtcs.array, old_vtcs[i].outvtcs.array, sizeof(int)*old_vtcs[i].outvtcs.used);
            MPIC_MPICH_free_vtx(&old_vtcs[i]);
        }
        MPIC_MPICH_free_mem(old_vtcs);
    }
    *vtx = &sched->vtcs[sched->total];
    return sched->total++;
}

/*This function should go away, keeping it there for smooth transition
from phase array to completely DAG based collectives*/
static inline int MPIC_MPICH_fence(MPIC_MPICH_sched_t *sched)
{
    //assert(sched->total < 32);
    if(0) fprintf(stderr, "TSP(mpich) : sched [fence] total=%ld \n", sched->total);
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    vtxp->kind = MPIC_MPICH_KIND_NOOP;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    int *invtcs = (int*)MPL_malloc(sizeof(int)*vtx_id);
    int i, n_invtcs=0;
    for(i=vtx_id-1; i>=0; i--){
        if(sched->vtcs[i].kind == MPIC_MPICH_KIND_NOOP)
            break;
        else{
            invtcs[vtx_id-1-i]=i;
            n_invtcs++;
        }
    }

    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    MPL_free(invtcs);
    return vtx_id;
}

/*MPIC_MPICH_wait waits for all the operations posted before it to complete
before issuing any operations posted after it. This is useful in composing
multiple schedules, for example, allreduce can be written as
COLL_sched_reduce(s)
MPIC_MPICH_wait(s)
COLL_sched_bcast(s)
This is different from the fence operation in the sense that fence requires
any vertex to post dependencies on it while MPIC_MPICH_wait is used internally
by the transport to add it as a dependency to any operations poster after it
*/

static inline int MPIC_MPICH_wait(MPIC_MPICH_sched_t *sched)
{
    if(0) fprintf(stderr, "scheduling a wait\n");
    sched->last_wait = sched->total;
    return MPIC_MPICH_fence(sched);
}

static inline void MPIC_MPICH_addref_op(MPIC_MPICH_op_t op,
                                 int       up)
{
    MPIR_Op     *op_ptr;

    if(HANDLE_GET_KIND(op) != HANDLE_KIND_BUILTIN) {
        MPIR_Op_get_ptr(op, op_ptr);

        if(up)
            MPIR_Op_add_ref(op_ptr);
        else
            MPIR_Op_release(op_ptr);
    }
}

static inline int MPIC_MPICH_addref_dt_nb(MPIC_MPICH_dt_t    dt,
                                    int          up,
                                    MPIC_MPICH_sched_t *sched,
                                    int          n_invtcs,
                                    int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    vtxp->kind                = MPIC_MPICH_KIND_ADDREF_DT;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    vtxp->nbargs.addref_dt.dt = dt;
    vtxp->nbargs.addref_dt.up = up;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [addref_dt]\n",vtx_id);
    return vtx_id;
}

static inline int MPIC_MPICH_addref_op_nb(MPIC_MPICH_op_t    op,
                                    int          up,
                                    MPIC_MPICH_sched_t *sched,
                                    int          n_invtcs,
                                    int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    vtxp->kind                = MPIC_MPICH_KIND_ADDREF_OP;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    vtxp->nbargs.addref_op.op = op;
    vtxp->nbargs.addref_op.up = up;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [addref_op]\n",vtx_id);
    return vtx_id;
}


static inline int MPIC_MPICH_test(MPIC_MPICH_sched_t *sched);
static inline int MPIC_MPICH_send(const void  *buf,
                            int          count,
                            MPIC_MPICH_dt_t    dt,
                            int          dest,
                            int          tag,
                            MPIC_MPICH_comm_t  *comm,
                            MPIC_MPICH_sched_t *sched,
                            int          n_invtcs,
                            int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    sched->tag = tag;
    vtxp->kind                  = MPIC_MPICH_KIND_SEND;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    vtxp->nbargs.sendrecv.buf   = (void *)buf;
    vtxp->nbargs.sendrecv.count = count;
    vtxp->nbargs.sendrecv.dt    = dt;
    vtxp->nbargs.sendrecv.dest  = dest;
    vtxp->nbargs.sendrecv.comm  = comm;
    vtxp->mpid_req[1]           = NULL;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [send]\n",vtx_id);
    return vtx_id;
}

static inline int MPIC_MPICH_send_accumulate(const void  *buf,
                                       int          count,
                                       MPIC_MPICH_dt_t    dt,
                                       MPIC_MPICH_op_t    op,
                                       int          dest,
                                       int          tag,
                                       MPIC_MPICH_comm_t  *comm,
                                       MPIC_MPICH_sched_t *sched,
                                       int          n_invtcs,
                                       int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    sched->tag = tag;
    vtxp->kind                  = MPIC_MPICH_KIND_SEND;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    vtxp->nbargs.sendrecv.buf   = (void *)buf;
    vtxp->nbargs.sendrecv.count = count;
    vtxp->nbargs.sendrecv.dt    = dt;
    vtxp->nbargs.sendrecv.dest  = dest;
    vtxp->nbargs.sendrecv.comm  = comm;
    vtxp->mpid_req[1]           = NULL;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [send_accumulate]\n",vtx_id);
    return vtx_id;
}

static inline int MPIC_MPICH_recv(void        *buf,
                            int          count,
                            MPIC_MPICH_dt_t    dt,
                            int          source,
                            int          tag,
                            MPIC_MPICH_comm_t  *comm,
                            MPIC_MPICH_sched_t *sched,
                            int          n_invtcs,
                            int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    sched->tag = tag;
    vtxp->kind                  = MPIC_MPICH_KIND_RECV;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    vtxp->nbargs.sendrecv.buf   = buf;
    vtxp->nbargs.sendrecv.count = count;
    vtxp->nbargs.sendrecv.dt    = dt;
    vtxp->nbargs.sendrecv.dest  = source;
    vtxp->nbargs.sendrecv.comm  = comm;
    vtxp->mpid_req[1]           = NULL;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [recv]\n",vtx_id);
    return vtx_id;
}

static inline int MPIC_MPICH_queryfcn(
    void       *data,
    MPI_Status *status)
{
    MPIC_MPICH_recv_reduce_arg_t *rr
        = (MPIC_MPICH_recv_reduce_arg_t *)data;
    if(rr->vtxp->mpid_req[0] == NULL && !rr->done) {
        MPI_Datatype dt = rr->datatype;
        MPI_Op       op = rr->op;

        if(rr->flags==-1 || rr->flags & MPIC_FLAG_REDUCE_L) {
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

        //MPL_free(rr->inbuf);
        MPIR_Grequest_complete_impl(rr->vtxp->mpid_req[1]);
        rr->done = 1;
    }

    status->MPI_SOURCE = MPI_UNDEFINED;
    status->MPI_TAG    = MPI_UNDEFINED;
    MPI_Status_set_cancelled(status, 0);
    MPI_Status_set_elements(status, MPI_BYTE, 0);
    return MPI_SUCCESS;
}

static inline int MPIC_MPICH_recv_reduce(void        *buf,
                                   int          count,
                                   MPIC_MPICH_dt_t    datatype,
                                   MPIC_MPICH_op_t    op,
                                   int          source,
                                   int          tag,
                                   MPIC_MPICH_comm_t  *comm,
                                   uint64_t     flags,
                                   MPIC_MPICH_sched_t *sched,
                                   int          n_invtcs,
                                   int         *invtcs)
{
    int        iscontig;
    size_t     type_size, out_extent, lower_bound;
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    sched->tag = tag;
    vtxp->kind                = MPIC_MPICH_KIND_RECV_REDUCE;
    MPIC_MPICH_init_vtx(sched, vtxp, vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);

    MPIC_MPICH_dtinfo(datatype,&iscontig,&type_size,&out_extent,&lower_bound);
    vtxp->nbargs.recv_reduce.inbuf    = MPL_malloc(count*out_extent);
    vtxp->nbargs.recv_reduce.inoutbuf = buf;
    vtxp->nbargs.recv_reduce.count    = count;
    vtxp->nbargs.recv_reduce.datatype = datatype;
    vtxp->nbargs.recv_reduce.op       = op;
    vtxp->nbargs.recv_reduce.source   = source;
    vtxp->nbargs.recv_reduce.comm     = comm;
    vtxp->nbargs.recv_reduce.vtxp      = vtxp;
    vtxp->nbargs.recv_reduce.done     = 0;
    vtxp->nbargs.recv_reduce.flags    = flags;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [recv_reduce]\n",vtx_id);
    
    return vtx_id;
}

static inline int MPIC_MPICH_rank(MPIC_MPICH_comm_t  *comm)
{
    return comm->mpid_comm->rank;
}

static inline int MPIC_MPICH_size(MPIC_MPICH_comm_t *comm)
{
    return comm->mpid_comm->local_size;
}

static inline int MPIC_MPICH_reduce_local(const void  *inbuf,
                                    void        *inoutbuf,
                                    int          count,
                                    MPIC_MPICH_dt_t    datatype,
                                    MPIC_MPICH_op_t    operation,
                                    uint64_t     flags,
                                    MPIC_MPICH_sched_t *sched,
                                    int          n_invtcs,
                                    int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    vtxp->kind                         = MPIC_MPICH_KIND_REDUCE_LOCAL;
    vtxp->state                        = MPIC_MPICH_STATE_INIT;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    vtxp->nbargs.reduce_local.inbuf    = inbuf;
    vtxp->nbargs.reduce_local.inoutbuf = inoutbuf;
    vtxp->nbargs.reduce_local.count    = count;
    vtxp->nbargs.reduce_local.dt       = datatype;
    vtxp->nbargs.reduce_local.op       = operation;
    vtxp->nbargs.reduce_local.flags    = flags;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [reduce_local]\n",vtx_id);
    return vtx_id;
}

static inline int MPIC_MPICH_dtcopy(void       *tobuf,
                             int         tocount,
                             MPIC_MPICH_dt_t   totype,
                             const void *frombuf,
                             int         fromcount,
                             MPIC_MPICH_dt_t   fromtype)
{
    return MPIR_Localcopy(frombuf,   /* yes, parameters are reversed        */
                          fromcount, /* MPICH forgot what memcpy looks like */
                          fromtype,
                          tobuf,
                          tocount,
                          totype);
}

static inline int MPIC_MPICH_dtcopy_nb(void        *tobuf,
                                 int          tocount,
                                 MPIC_MPICH_dt_t    totype,
                                 const void  *frombuf,
                                 int          fromcount,
                                 MPIC_MPICH_dt_t    fromtype,
                                 MPIC_MPICH_sched_t *sched,
                                 int          n_invtcs,
                                 int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    vtxp->kind   = MPIC_MPICH_KIND_DTCOPY;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);

    vtxp->nbargs.dtcopy.tobuf      = tobuf;
    vtxp->nbargs.dtcopy.tocount    = tocount;
    vtxp->nbargs.dtcopy.totype     = totype;
    vtxp->nbargs.dtcopy.frombuf    = frombuf;
    vtxp->nbargs.dtcopy.fromcount  = fromcount;
    vtxp->nbargs.dtcopy.fromtype   = fromtype;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [dt_copy]\n",vtx_id);
    return vtx_id;
}

static inline void MPIC_MPICH_add_elem_ptr_array(MPIC_MPICH_ptr_array *buf_array, void *buf){
    if(buf_array->used == buf_array->size){
        int old_size = buf_array->size;
        int old_array = buf_array->array;
        buf_array->size *= 2;
        buf_array->array = (void**)MPIC_MPICH_allocate_mem(sizeof(void*)*buf_array->size);
        memcpy(buf_array->array, old_array, sizeof(void*)*old_size);
        MPIC_MPICH_free_mem(old_array);
    }
    buf_array->array[buf_array->used++] = buf;
}

/*This function allocates memory required for schedule execution.
 *This is recorded in the schedule so that this memory can be 
 *freed when the schedule is destroyed
 */
static inline void *MPIC_MPICH_allocate_buffer(size_t size, MPIC_MPICH_sched_t *s){
    void *buf = MPIC_MPICH_allocate_mem(size);
    /*record memory allocation*/
    MPIC_MPICH_add_elem_ptr_array(&s->buf_array, buf);
    return buf;
}


static inline int MPIC_MPICH_free_mem_nb(void        *ptr,
                                   MPIC_MPICH_sched_t *sched,
                                   int          n_invtcs,
                                   int         *invtcs)
{
    MPIC_MPICH_vtx_t *vtxp;
    int vtx_id = MPIC_MPICH_get_new_vtx(sched, &vtxp);
    vtxp->kind                = MPIC_MPICH_KIND_FREE_MEM;
    MPIC_MPICH_init_vtx(sched,vtxp,vtx_id);
    MPIC_MPICH_add_vtx_dependencies(sched, vtx_id, n_invtcs, invtcs);
    vtxp->nbargs.free_mem.ptr = ptr;

    if(0) fprintf(stderr, "TSP(mpich) : sched [%ld] [free_mem]\n",vtx_id);
    return vtx_id;
}

static inline void MPIC_MPICH_issue_vtx(int vtxid, MPIC_MPICH_vtx_t *rp, MPIC_MPICH_sched_t *sched){
    if(rp->state == MPIC_MPICH_STATE_INIT && rp->num_unfinished_dependencies==0) {
        if(0) fprintf(stderr, "issuing request: %d\n", vtxid);
        switch(rp->kind) {
            case MPIC_MPICH_KIND_SEND: {
                if(0) fprintf(stderr, "  --> MPICH transport (isend) issued\n");

                MPIR_Errflag_t  errflag  = MPIR_ERR_NONE;
                MPIC_Isend(rp->nbargs.sendrecv.buf,
                           rp->nbargs.sendrecv.count,
                           rp->nbargs.sendrecv.dt,
                           rp->nbargs.sendrecv.dest,
                           sched->tag,
                           rp->nbargs.sendrecv.comm->mpid_comm,
                           &rp->mpid_req[0],
                           &errflag);
                MPIC_MPICH_record_vtx_issue(rp, sched);
            }
            break;

            case MPIC_MPICH_KIND_RECV: {
                if(0) fprintf(stderr, "  --> MPICH transport (irecv) issued\n");

                MPIC_Irecv(rp->nbargs.sendrecv.buf,
                           rp->nbargs.sendrecv.count,
                           rp->nbargs.sendrecv.dt,
                           rp->nbargs.sendrecv.dest,
                           sched->tag,
                           rp->nbargs.sendrecv.comm->mpid_comm,
                           &rp->mpid_req[0]);
                MPIC_MPICH_record_vtx_issue(rp, sched);
            }
            break;

            case MPIC_MPICH_KIND_ADDREF_DT:
                MPIC_MPICH_addref_dt(rp->nbargs.addref_dt.dt,
                              rp->nbargs.addref_dt.up);

                if(0) fprintf(stderr, "  --> MPICH transport (addref dt) complete\n");
                MPIC_MPICH_record_vtx_completion(rp, sched);
                break;

            case MPIC_MPICH_KIND_ADDREF_OP:
                MPIC_MPICH_addref_op(rp->nbargs.addref_op.op,
                              rp->nbargs.addref_op.up);


                if(0) fprintf(stderr, "  --> MPICH transport (addref op) complete\n");

                MPIC_MPICH_record_vtx_completion(rp, sched);
                break;

            case MPIC_MPICH_KIND_DTCOPY:
                MPIC_MPICH_dtcopy(rp->nbargs.dtcopy.tobuf,
                           rp->nbargs.dtcopy.tocount,
                           rp->nbargs.dtcopy.totype,
                           rp->nbargs.dtcopy.frombuf,
                           rp->nbargs.dtcopy.fromcount,
                           rp->nbargs.dtcopy.fromtype);

                if(0) fprintf(stderr, "  --> MPICH transport (dtcopy) complete\n");
                MPIC_MPICH_record_vtx_completion(rp, sched);

                break;

            case MPIC_MPICH_KIND_FREE_MEM:
                if(0) fprintf(stderr, "  --> MPICH transport (freemem) complete\n");

                MPIC_MPICH_free_mem(rp->nbargs.free_mem.ptr);
                MPIC_MPICH_record_vtx_completion(rp, sched);
                break;

            case MPIC_MPICH_KIND_NOOP:
                if(0) fprintf(stderr, "  --> MPICH transport (noop) complete\n");

                MPIC_MPICH_record_vtx_completion(rp, sched);
                break;

            case MPIC_MPICH_KIND_RECV_REDUCE: {
                MPIC_Irecv(rp->nbargs.recv_reduce.inbuf,
                           rp->nbargs.recv_reduce.count,
                           rp->nbargs.recv_reduce.datatype,
                           rp->nbargs.recv_reduce.source,
                           sched->tag,
                           rp->nbargs.recv_reduce.comm->mpid_comm,
                           &rp->mpid_req[0]);
                MPIR_Grequest_start_impl(MPIC_MPICH_queryfcn,
                                         NULL,
                                         NULL,
                                         &rp->nbargs.recv_reduce,
                                         &rp->mpid_req[1]);

                if(0) fprintf(stderr, "  --> MPICH transport (recv_reduce) issued\n");

                MPIC_MPICH_record_vtx_issue(rp, sched);
            }
            break;

            case MPIC_MPICH_KIND_REDUCE_LOCAL:
                  if(rp->nbargs.reduce_local.flags==-1 || rp->nbargs.reduce_local.flags & MPIC_FLAG_REDUCE_L) {
                          if(0) fprintf("rp->nbargs.reduce_local.flags %d \n",rp->nbargs.reduce_local.flags);
                          MPIR_Reduce_local_impl(rp->nbargs.reduce_local.inbuf,
                                                 rp->nbargs.reduce_local.inoutbuf,
                                                 rp->nbargs.reduce_local.count,
                                                 rp->nbargs.reduce_local.dt,
                                                 rp->nbargs.reduce_local.op);
                          if(0) fprintf(stderr, "  --> MPICH transport (reduce local_L) complete\n");
                   } else {
                          printf("Right reduction rp->nbargs.reduce_local.flags %d \n",rp->nbargs.reduce_local.flags);
                          MPIR_Reduce_local_impl(rp->nbargs.reduce_local.inoutbuf,
                                                 rp->nbargs.reduce_local.inbuf,
                                                 rp->nbargs.reduce_local.count,
                                                 rp->nbargs.reduce_local.dt,
                                                 rp->nbargs.reduce_local.op);

                          MPIR_Localcopy(rp->nbargs.reduce_local.inbuf,
                                         rp->nbargs.reduce_local.count,
                                         rp->nbargs.reduce_local.dt,
                                         rp->nbargs.reduce_local.inoutbuf,
                                         rp->nbargs.reduce_local.count,
                                         rp->nbargs.reduce_local.dt);
                          if(0) fprintf(stderr, "  --> MPICH transport (reduce local_R) complete\n");
                   }

                MPIC_MPICH_record_vtx_completion(rp, sched);
                break;
        }
    }
}


static inline int MPIC_MPICH_test(MPIC_MPICH_sched_t *sched)
{
    MPIC_MPICH_vtx_t *req, *rp;
    int i;
    req = &sched->vtcs[0];
    //if(0)fprintf(stderr, "in TSP_test, num_completed=%d, total=%d\n", sched->num_completed, sched->total);
    /*if issued list is empty, generate it*/
    if(sched->issued_head == NULL){
        if(0) fprintf(stderr, "issued list is empty, issue ready vtcs\n");
        if(sched->total > 0 && sched->num_completed != sched->total){
            for(i=0; i<sched->total; i++)
                MPIC_MPICH_issue_vtx(i, &sched->vtcs[i], sched);
            if(0) fprintf(stderr, "completed traversal of vtcs, sched->total: %d, sched->num_completed: %d\n", sched->total, sched->num_completed);
            return 0;
        }
        else
            return 1;
    }
    if(sched->total == sched->num_completed){
        return 1;
    }
        /*fprintf(stderr, "issued list: ");
        TSP_req_t *tmp = sched->issued_head;
        while(tmp){fprintf(stderr, "%d ", req->kind); tmp = tmp->next_issued;}
        fprintf(stderr,"\n");*/
    
    assert(sched->issued_head != NULL);
    sched->vtx_iter = sched->issued_head;
    sched->issued_head = NULL;
    /* Check for issued ops that have been completed */
    while(sched->vtx_iter!=NULL) {
        rp = sched->vtx_iter;
        sched->vtx_iter = sched->vtx_iter->next_issued;

        if(rp->state == MPIC_MPICH_STATE_ISSUED) {
            MPI_Status     status;
            MPIR_Request  *mpid_req0 = rp->mpid_req[0];
            MPIR_Request  *mpid_req1 = rp->mpid_req[1];

            if(mpid_req1) {
                (mpid_req1->u.ureq.greq_fns->query_fn)
                (mpid_req1->u.ureq.greq_fns->grequest_extra_state,
                 &status);
            }

            switch(rp->kind) {
                case MPIC_MPICH_KIND_SEND:
                case MPIC_MPICH_KIND_RECV:
                    if(mpid_req0 && MPIR_Request_is_complete(mpid_req0)) {
                        MPIR_Request_free(mpid_req0);
                        rp->mpid_req[0] = NULL;
                    }
                    if(!rp->mpid_req[0]) {
                        if(0){ fprintf(stderr, "  --> MPICH transport (kind=%d) complete\n",
                                          rp->kind);
                             if(rp->nbargs.sendrecv.count>=1) fprintf(stderr, "data send/recvd: %d\n", *(int *)(rp->nbargs.sendrecv.buf));
                        }
                        MPIC_MPICH_record_vtx_completion(rp,sched);
                    }else
                        MPIC_MPICH_record_vtx_issue(rp,sched); /*record it again as issued*/
                    break;
                case MPIC_MPICH_KIND_RECV_REDUCE:
                    if(mpid_req0 && MPIR_Request_is_complete(mpid_req0)) {
                        MPIR_Request_free(mpid_req0);
                        if(0)fprintf(stderr,"recv in recv_reduce completed\n");
                        rp->mpid_req[0] = NULL;
                    }

                    if(mpid_req1 && MPIR_Request_is_complete(mpid_req1)) {
                        MPIR_Request_free(mpid_req1);
                        rp->mpid_req[1] = NULL;
                    }

                    if(!rp->mpid_req[0] && !rp->mpid_req[1]) {
                        if(0){ fprintf(stderr, "  --> MPICH transport (kind=%d) complete\n",
                                          rp->kind);
                             if(rp->nbargs.sendrecv.count>=1) fprintf(stderr, "data send/recvd: %d\n", *(int *)(rp->nbargs.sendrecv.buf));
                        }
                        MPIC_MPICH_record_vtx_completion(rp,sched);
                    }else
                        MPIC_MPICH_record_vtx_issue(rp,sched);/*record it again as issued*/

                    break;

                default:
                    break;
            }
        }
        //sched->vtx_iter = sched->vtx_iter->next_issued;
    }
    sched->last_issued->next_issued = NULL;

    if(0){

        if(sched->num_completed==sched->total) {
            if(0) fprintf(stderr, "  --> MPICH transport (test) complete:  sched->total=%ld\n",
                              sched->total);
        }
    }
    return sched->num_completed==sched->total;
}

 /*frees any memory allocated for execution of this schedule*/
static inline void MPIC_MPICH_free_buffers(MPIC_MPICH_sched_t *sched){
    int i;
    for(i=0; i<sched->total; i++){
        /*free the temporary memory allocated by recv_reduce call*/
        if(sched->vtcs[i].kind == MPIC_MPICH_KIND_RECV_REDUCE){
            MPIC_MPICH_free_mem(sched->vtcs[i].nbargs.recv_reduce.inbuf);
        }
    }
    /*free temporary buffers*/
    for(i=0; i<sched->buf_array.used; i++){
        MPIC_MPICH_free_mem(sched->buf_array.array[i]);
    }
    MPIC_MPICH_free_mem(sched->buf_array.array);

    /*free each vtx and then the list of vtcs*/
    for(i=0; i<sched->total; i++){ /*up to sched->total because we init vertices only when we need them*/
        MPIC_MPICH_free_vtx(&sched->vtcs[i]);
    }
    MPIC_MPICH_free_mem(sched->vtcs);
}



#endif
