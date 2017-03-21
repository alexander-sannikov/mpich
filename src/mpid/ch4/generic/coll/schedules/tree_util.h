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


static inline int COLL_ilog(int k, int number)
{
    int i = 1, p = k - 1;

    for (; p - 1 < number; i++)
        p *= k;

    return i;
}

static inline int COLL_ipow(int base, int exp)
{
    int result = 1;

    while (exp) {
        if (exp & 1)
            result *= base;

        exp >>= 1;
        base *= base;
    }

    return result;
}

static inline int COLL_getdigit(int k, int number, int digit)
{
    return (number / COLL_ipow(k, digit)) % k;
}

static inline int COLL_setdigit(int k, int number, int digit, int newdigit)
{
    int res = number;
    int lshift = COLL_ipow(k, digit);
    res -= COLL_getdigit(k, number, digit) * lshift;
    res += newdigit * lshift;
    return res;
}

static inline void COLL_tree_add_child(COLL_tree_t * t, int rank)
{
    if (t->numRanges > 0 && t->children[t->numRanges - 1].endRank == rank - 1)
        t->children[t->numRanges - 1].endRank = rank;
    else {
        t->numRanges++;
        COLL_Assert(t->numRanges < COLL_MAX_TREE_BREADTH);
        t->children[t->numRanges - 1].startRank = rank;
        t->children[t->numRanges - 1].endRank = rank;
    }
}

static inline void COLL_tree_kary_init(int rank, int nranks, int k, int root, COLL_tree_t * ct)
{
    ct->rank = rank;
    ct->nranks = nranks;
    ct->numRanges = 0;
    ct->parent = -1;

    if (nranks <= 0)
        return;

    int lrank = (rank + (nranks - root)) % nranks;
    int child;
    COLL_Assert(k >= 1);

    ct->parent = (lrank <= 0) ? -1 : (((lrank - 1) / k) + root) % nranks;

    for (child = 1; child <= k; child++) {
        int val = lrank * k + child;
        if (val >= nranks)
            break;
        val = (val + root) % nranks;
        COLL_tree_add_child(ct, val);
    }
}

static inline void COLL_tree_knomial_init(int rank, int nranks, int k, int root, COLL_tree_t * ct)
{
    ct->rank = rank;
    ct->nranks = nranks;
    ct->numRanges = 0;
    ct->parent = -1;

    if (nranks <= 0)
        return;

    int lrank = (rank + (nranks - root)) % nranks;
    int basek, i, j;
    COLL_Assert(k >= 2);

    /* Parent Calculation */
    if (lrank <= 0)
        ct->parent = -1;
    else {
        basek = COLL_ilog(k, nranks - 1);

        for (i = 0; i < basek; i++) {
            if (COLL_getdigit(k, lrank, i)) {
                ct->parent = (COLL_setdigit(k, lrank, i, 0) + root) % nranks;
                break;
            }
        }
    }

    if (lrank >= nranks)
        return;

    /* Children Calculation */
    basek = COLL_ilog(k, nranks - 1);

    for (j = 0; j < basek; j++) {
        if (COLL_getdigit(k, lrank, j))
            break;

        for (i = 1; i < k; i++) {
            int child = COLL_setdigit(k, lrank, j, i);
            if (child < nranks)
                COLL_tree_add_child(ct, (child + root) % nranks);
        }
    }
}
