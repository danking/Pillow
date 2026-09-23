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

    static inline UINT8 RankMedian3x3UINT8(
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

static inline UINT8
RankMedian5x5UINT8(const UINT8 p[25]) {
    UINT8 v[25];
    int i, j;

    memcpy(v, p, sizeof(v));

#define RANK_COMPARE_SWAP_UINT8_INDEX(a, b)    \
    do {                                       \
        UINT8 lo_ = v[a] < v[b] ? v[a] : v[b]; \
        UINT8 hi_ = v[a] < v[b] ? v[b] : v[a]; \
        v[a] = lo_;                            \
        v[b] = hi_;                            \
    } while (0)

    for (i = 0; i <= 12; i++) {
        for (j = i + 1; j < 25; j++) {
            RANK_COMPARE_SWAP_UINT8_INDEX(i, j);
        }
    }

#undef RANK_COMPARE_SWAP_UINT8_INDEX
    return v[12];
}

/*
 * This network is the scalar selection network specialized for five sorted
 * columns. Inputs use p[5 * horizontal_column + vertical_rank], and p12 is
 * the median. Keeping the operations scalar lets the compiler vectorize the
 * independent output pixels without architecture-specific types or intrinsics.
 */
#define RANK_MEDIAN_25_SORTED_COLUMNS_NETWORK(sort2) \
    sort2(4, 5);                                     \
    sort2(14, 15);                                   \
    sort2(5, 7);                                     \
    sort2(8, 10);                                    \
    sort2(9, 11);                                    \
    sort2(12, 14);                                   \
    sort2(0, 4);                                     \
    sort2(1, 5);                                     \
    sort2(2, 6);                                     \
    sort2(8, 12);                                    \
    sort2(10, 14);                                   \
    sort2(11, 15);                                   \
    sort2(16, 20);                                   \
    sort2(17, 21);                                   \
    sort2(18, 22);                                   \
    sort2(19, 23);                                   \
    sort2(0, 8);                                     \
    sort2(1, 9);                                     \
    sort2(2, 10);                                    \
    sort2(3, 11);                                    \
    sort2(4, 12);                                    \
    sort2(5, 13);                                    \
    sort2(6, 14);                                    \
    sort2(7, 15);                                    \
    sort2(3, 19);                                    \
    sort2(5, 21);                                    \
    sort2(6, 22);                                    \
    sort2(7, 23);                                    \
    sort2(8, 16);                                    \
    sort2(9, 17);                                    \
    sort2(10, 18);                                   \
    sort2(12, 20);                                   \
    sort2(1, 12);                                    \
    sort2(2, 9);                                     \
    sort2(3, 5);                                     \
    sort2(4, 10);                                    \
    sort2(11, 21);                                   \
    sort2(13, 22);                                   \
    sort2(14, 19);                                   \
    sort2(18, 20);                                   \
    sort2(9, 10);                                    \
    sort2(11, 13);                                   \
    sort2(7, 13);                                    \
    sort2(10, 12);                                   \
    sort2(11, 14);                                   \
    sort2(5, 11);                                    \
    sort2(10, 16);                                   \
    sort2(12, 18);                                   \
    sort2(3, 12);                                    \
    sort2(9, 16);                                    \
    sort2(11, 20);                                   \
    sort2(14, 18);                                   \
    sort2(5, 16);                                    \
    sort2(6, 12);                                    \
    sort2(7, 14);                                    \
    sort2(11, 17);                                   \
    sort2(7, 16);                                    \
    sort2(11, 24);                                   \
    sort2(12, 17);                                   \
    sort2(7, 11);                                    \
    sort2(16, 24);                                   \
    sort2(6, 11);                                    \
    sort2(12, 16);                                   \
    sort2(11, 12);                                   \
    sort2(12, 16)

/* The block length is a batching choice, not a compiler vector type. */
#define RANK_AUTOVEC_LANES 16

static inline void
RankSort5UINT8Block(UINT8 p[5][2 * RANK_AUTOVEC_LANES], size_t offset) {
    for (size_t lane = 0; lane < RANK_AUTOVEC_LANES; lane++) {
#define RANK_COMPARE_SWAP_SORT5_UINT8(a, b) \
    do {                                    \
        UINT8 lo_ = (a) < (b) ? (a) : (b);  \
        UINT8 hi_ = (a) < (b) ? (b) : (a);  \
        (a) = lo_;                          \
        (b) = hi_;                          \
    } while (0)
        UINT8 p0 = p[0][offset + lane];
        UINT8 p1 = p[1][offset + lane];
        UINT8 p2 = p[2][offset + lane];
        UINT8 p3 = p[3][offset + lane];
        UINT8 p4 = p[4][offset + lane];

        RANK_COMPARE_SWAP_SORT5_UINT8(p0, p3);
        RANK_COMPARE_SWAP_SORT5_UINT8(p1, p4);
        RANK_COMPARE_SWAP_SORT5_UINT8(p0, p2);
        RANK_COMPARE_SWAP_SORT5_UINT8(p1, p3);
        RANK_COMPARE_SWAP_SORT5_UINT8(p0, p1);
        RANK_COMPARE_SWAP_SORT5_UINT8(p2, p4);
        RANK_COMPARE_SWAP_SORT5_UINT8(p1, p2);
        RANK_COMPARE_SWAP_SORT5_UINT8(p3, p4);
        RANK_COMPARE_SWAP_SORT5_UINT8(p2, p3);

        p[0][offset + lane] = p0;
        p[1][offset + lane] = p1;
        p[2][offset + lane] = p2;
        p[3][offset + lane] = p3;
        p[4][offset + lane] = p4;
#undef RANK_COMPARE_SWAP_SORT5_UINT8
    }
}

static inline UINT8
RankMedian5x5SortedColumnsUINT8(
    UINT8 p0,
    UINT8 p1,
    UINT8 p2,
    UINT8 p3,
    UINT8 p4,
    UINT8 p5,
    UINT8 p6,
    UINT8 p7,
    UINT8 p8,
    UINT8 p9,
    UINT8 p10,
    UINT8 p11,
    UINT8 p12,
    UINT8 p13,
    UINT8 p14,
    UINT8 p15,
    UINT8 p16,
    UINT8 p17,
    UINT8 p18,
    UINT8 p19,
    UINT8 p20,
    UINT8 p21,
    UINT8 p22,
    UINT8 p23,
    UINT8 p24
) {
#define RANK_COMPARE_SWAP_SORTED_UINT8(a, b)   \
    do {                                       \
        UINT8 lo_ = p##a < p##b ? p##a : p##b; \
        UINT8 hi_ = p##a < p##b ? p##b : p##a; \
        p##a = lo_;                            \
        p##b = hi_;                            \
    } while (0)

    RANK_MEDIAN_25_SORTED_COLUMNS_NETWORK(RANK_COMPARE_SWAP_SORTED_UINT8);

#undef RANK_COMPARE_SWAP_SORTED_UINT8
    return p12;
}

static void
RankMedian5x5SortedColumnsUINT8Block(
    UINT8 *restrict out, const UINT8 p[5][2 * RANK_AUTOVEC_LANES], size_t lanes
) {
    for (size_t lane = 0; lane < lanes; lane++) {
        out[lane] = RankMedian5x5SortedColumnsUINT8(
            p[0][lane],
            p[1][lane],
            p[2][lane],
            p[3][lane],
            p[4][lane],
            p[0][lane + 1],
            p[1][lane + 1],
            p[2][lane + 1],
            p[3][lane + 1],
            p[4][lane + 1],
            p[0][lane + 2],
            p[1][lane + 2],
            p[2][lane + 2],
            p[3][lane + 2],
            p[4][lane + 2],
            p[0][lane + 3],
            p[1][lane + 3],
            p[2][lane + 3],
            p[3][lane + 3],
            p[4][lane + 3],
            p[0][lane + 4],
            p[1][lane + 4],
            p[2][lane + 4],
            p[3][lane + 4],
            p[4][lane + 4]
        );
    }
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

static void
RankFilter5x5MedianUINT8Row(
    UINT8 *restrict out,
    const UINT8 *row0,
    const UINT8 *row1,
    const UINT8 *row2,
    const UINT8 *row3,
    const UINT8 *row4,
    size_t xsize
) {
    size_t x = 0;

    if (xsize >= 64) {
        UINT8 columns[5][2 * RANK_AUTOVEC_LANES];
        UINT8 result[RANK_AUTOVEC_LANES];
        const UINT8 *rows[5] = {row0, row1, row2, row3, row4};

        for (size_t rank = 0; rank < 5; rank++) {
            memcpy(columns[rank], rows[rank], RANK_AUTOVEC_LANES);
        }
        RankSort5UINT8Block(columns, 0);

        for (; x <= xsize - 28; x += RANK_AUTOVEC_LANES) {
            for (size_t rank = 0; rank < 5; rank++) {
                memcpy(
                    columns[rank] + RANK_AUTOVEC_LANES,
                    rows[rank] + x + RANK_AUTOVEC_LANES,
                    RANK_AUTOVEC_LANES
                );
            }
            RankSort5UINT8Block(columns, RANK_AUTOVEC_LANES);
            RankMedian5x5SortedColumnsUINT8Block(result, columns, RANK_AUTOVEC_LANES);
            memcpy(out + x, result, RANK_AUTOVEC_LANES);
            for (size_t rank = 0; rank < 5; rank++) {
                memcpy(
                    columns[rank],
                    columns[rank] + RANK_AUTOVEC_LANES,
                    RANK_AUTOVEC_LANES
                );
            }
        }
    }

    for (; x < xsize; x++) {
        UINT8 p[25] = {
            row0[x], row0[x + 1], row0[x + 2], row0[x + 3], row0[x + 4],
            row1[x], row1[x + 1], row1[x + 2], row1[x + 3], row1[x + 4],
            row2[x], row2[x + 1], row2[x + 2], row2[x + 3], row2[x + 4],
            row3[x], row3[x + 1], row3[x + 2], row3[x + 3], row3[x + 4],
            row4[x], row4[x + 1], row4[x + 2], row4[x + 3], row4[x + 4],
        };
        out[x] = RankMedian5x5UINT8(p);
    }
}

static void
RankFilter5x5MedianUINT8(Imaging imOut, Imaging im) {
    const size_t xsize = (size_t)imOut->xsize;

    /* imOut is a fresh allocation, so its output row cannot alias im. */
    for (int y = 0; y < imOut->ysize; y++) {
        RankFilter5x5MedianUINT8Row(
            (UINT8 *)imOut->image[y],
            (UINT8 *)im->image[y],
            (UINT8 *)im->image[y + 1],
            (UINT8 *)im->image[y + 2],
            (UINT8 *)im->image[y + 3],
            (UINT8 *)im->image[y + 4],
            xsize
        );
    }
}

Imaging
ImagingRankFilter(Imaging im, int size, int rank) {
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
    } else if (im->image8 && size == 5 && rank == 12) {
        ImagingSectionEnter(&cookie);
        RankFilter5x5MedianUINT8(imOut, im);
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
