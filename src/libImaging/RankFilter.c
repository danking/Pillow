/*
 * The Python Imaging Library
 * $Id$
 *
 * min, max, median filters
 *
 * history:
 * 2002-06-08 fl    Created
 *
 * Copyright (c) Secret Labs AB 2002.  All rights reserved.
 *
 * See the README file for information on usage and redistribution.
 */

#include "Imaging.h"

/* Fast rank algorithm (due to Wirth), based on public domain code
   by Nicolas Devillard, available at http://ndevilla.free.fr */

#define SWAP(type, a, b)       \
    {                          \
        register type t = (a); \
        (a) = (b);             \
        (b) = t;               \
    }

#define MakeRankFunction(type)                       \
    static type Rank##type(type a[], int n, int k) { \
        register int i, j, l, m;                     \
        register type x;                             \
        l = 0;                                       \
        m = n - 1;                                   \
        while (l < m) {                              \
            x = a[k];                                \
            i = l;                                   \
            j = m;                                   \
            do {                                     \
                while (a[i] < x) {                   \
                    i++;                             \
                }                                    \
                while (x < a[j]) {                   \
                    j--;                             \
                }                                    \
                if (i <= j) {                        \
                    SWAP(type, a[i], a[j]);          \
                    i++;                             \
                    j--;                             \
                }                                    \
            } while (i <= j);                        \
            if (j < k) {                             \
                l = i;                               \
            }                                        \
            if (k < i) {                             \
                m = j;                               \
            }                                        \
        }                                            \
        return a[k];                                 \
    }

MakeRankFunction(UINT8) MakeRankFunction(INT32) MakeRankFunction(FLOAT32)

static inline UINT8
RankMedian3x3UINT8(
    UINT8 p0,
    UINT8 p1,
    UINT8 p2,
    UINT8 p3,
    UINT8 p4,
    UINT8 p5,
    UINT8 p6,
    UINT8 p7,
    UINT8 p8
) {
#define RANK_COMPARE_SWAP_UINT8(a, b)      \
    do {                                   \
        UINT8 lo_ = (a) < (b) ? (a) : (b); \
        UINT8 hi_ = (a) < (b) ? (b) : (a); \
        (a) = lo_;                         \
        (b) = hi_;                         \
    } while (0)

    RANK_COMPARE_SWAP_UINT8(p1, p2);
    RANK_COMPARE_SWAP_UINT8(p4, p5);
    RANK_COMPARE_SWAP_UINT8(p7, p8);
    RANK_COMPARE_SWAP_UINT8(p0, p1);
    RANK_COMPARE_SWAP_UINT8(p3, p4);
    RANK_COMPARE_SWAP_UINT8(p6, p7);
    RANK_COMPARE_SWAP_UINT8(p1, p2);
    RANK_COMPARE_SWAP_UINT8(p4, p5);
    RANK_COMPARE_SWAP_UINT8(p7, p8);
    RANK_COMPARE_SWAP_UINT8(p0, p3);
    RANK_COMPARE_SWAP_UINT8(p5, p8);
    RANK_COMPARE_SWAP_UINT8(p4, p7);
    RANK_COMPARE_SWAP_UINT8(p3, p6);
    RANK_COMPARE_SWAP_UINT8(p1, p4);
    RANK_COMPARE_SWAP_UINT8(p2, p5);
    RANK_COMPARE_SWAP_UINT8(p4, p7);
    RANK_COMPARE_SWAP_UINT8(p4, p2);
    RANK_COMPARE_SWAP_UINT8(p6, p4);
    RANK_COMPARE_SWAP_UINT8(p4, p2);

#undef RANK_COMPARE_SWAP_UINT8
    return p4;
}

static void
RankFilter3x3MedianUINT8Row(
    UINT8 *restrict out,
    const UINT8 *restrict row0,
    const UINT8 *restrict row1,
    const UINT8 *restrict row2,
    size_t xsize
) {
    for (size_t x = 0; x < xsize; x++) {
        out[x] = RankMedian3x3UINT8(
            row0[x],
            row0[x + 1],
            row0[x + 2],
            row1[x],
            row1[x + 1],
            row1[x + 2],
            row2[x],
            row2[x + 1],
            row2[x + 2]
        );
    }
}

static void
RankFilter3x3MedianUINT8(Imaging imOut, Imaging im) {
    // restrict is safe in the row helper: im is read-only and imOut is a fresh
    // allocation, while distinct input rows are only read.
    const size_t xsize = (size_t)imOut->xsize;
    for (int y = 0; y < imOut->ysize; y++) {
        RankFilter3x3MedianUINT8Row(
            (UINT8 *)imOut->image[y],
            (UINT8 *)im->image[y],
            (UINT8 *)im->image[y + 1],
            (UINT8 *)im->image[y + 2],
            xsize
        );
    }
}

    Imaging ImagingRankFilter(Imaging im, int size, int rank) {
    Imaging imOut = NULL;
    ImagingSectionCookie cookie;
    int x, y;
    int i, margin, size2;

    if (!im || im->bands != 1 || im->type == IMAGING_TYPE_I16) {
        return (Imaging)ImagingError_ModeError();
    }

    if (!(size & 1)) {
        return (Imaging)ImagingError_ValueError("bad filter size");
    }

    /* malloc check ok, for overflow in the define below */
    if (size > INT_MAX / (size * (int)sizeof(FLOAT32))) {
        return (Imaging)ImagingError_ValueError("filter size too large");
    }

    size2 = size * size;
    margin = (size - 1) / 2;

    if (rank < 0 || rank >= size2) {
        return (Imaging)ImagingError_ValueError("bad rank value");
    }

    // Every output pixel is written by the rank loop below
    imOut = ImagingNewDirty(im->mode, im->xsize - 2 * margin, im->ysize - 2 * margin);
    if (!imOut) {
        return NULL;
    }

    /* malloc check ok, checked above */
#define RANK_BODY(type)                                                           \
    do {                                                                          \
        type *buf = malloc(size2 * sizeof(type));                                 \
        if (!buf) {                                                               \
            goto nomemory;                                                        \
        }                                                                         \
        ImagingSectionEnter(&cookie);                                             \
        for (y = 0; y < imOut->ysize; y++) {                                      \
            for (x = 0; x < imOut->xsize; x++) {                                  \
                for (i = 0; i < size; i++) {                                      \
                    memcpy(                                                       \
                        buf + i * size,                                           \
                        &IMAGING_PIXEL_##type(im, x, y + i),                      \
                        size * sizeof(type)                                       \
                    );                                                            \
                }                                                                 \
                IMAGING_PIXEL_##type(imOut, x, y) = Rank##type(buf, size2, rank); \
            }                                                                     \
        }                                                                         \
        ImagingSectionLeave(&cookie);                                             \
        free(buf);                                                                \
    } while (0)

    if (im->image8 && size == 3 && rank == 4) {
        ImagingSectionEnter(&cookie);
        RankFilter3x3MedianUINT8(imOut, im);
        ImagingSectionLeave(&cookie);
    } else if (im->image8) {
        RANK_BODY(UINT8);
    } else if (im->type == IMAGING_TYPE_INT32) {
        RANK_BODY(INT32);
    } else if (im->type == IMAGING_TYPE_FLOAT32) {
        RANK_BODY(FLOAT32);
    } else {
        /* safety net (we shouldn't end up here) */
        ImagingDelete(imOut);
        return (Imaging)ImagingError_ModeError();
    }

    ImagingCopyPalette(imOut, im);

    return imOut;

nomemory:
    ImagingDelete(imOut);
    return (Imaging)ImagingError_MemoryError();
}
