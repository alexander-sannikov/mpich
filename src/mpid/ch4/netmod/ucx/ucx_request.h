/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2016 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Mellanox Technologies Ltd.
 *  Copyright (C) Mellanox Technologies Ltd. 2016. ALL RIGHTS RESERVED
 */
#ifndef UCX_REQUEST_H_INCLUDED
#define UCX_REQUEST_H_INCLUDED

#include "ucx_impl.h"

#undef FUNCNAME
#define FUNCNAME MPIDI_NM_am_request_init
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
MPL_STATIC_INLINE_PREFIX void MPIDI_NM_am_request_init(MPIR_Request * req)
{
    req->dev.ch4.am.netmod_am.ucx.pack_buffer = NULL;
}

#undef FUNCNAME
#define FUNCNAME MPIDI_NM_am_request_finalize
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
MPL_STATIC_INLINE_PREFIX void MPIDI_NM_am_request_finalize(MPIR_Request * req)
{
    MPL_free((req)->dev.ch4.am.netmod_am.ucx.pack_buffer);
    /* MPIR_Request_free(req); */
}

#endif /* UCX_REQUEST_H_INCLUDED */
