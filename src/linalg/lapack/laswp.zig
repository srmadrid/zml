const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const lapack = @import("../lapack.zig");

fn laswp(comptime A: type, order: Order, n: isize, a: [*]A,
                                lda: isize, k1: isize, k2: isize,
                                ipiv: [*]const isize, incx: isize ) isize
{
    var info: isize = 0;
    if(order == LAPACK_COL_MAJOR ) {
        // Call LAPACK function and adjust info
        k_laswp(A, n, a, lda, k1, k2, ipiv, incx );

        if( info < 0 ) {
            info -= 1;
        }
    } else if(order == LAPACK_ROW_MAJOR ) {
        var lda_t: isize = int.max(1, k2);
        var i = k1;
        while(i <= k2) {
            lda_t = int.max( lda_t, ipiv[scast(usize, k1 + ( i - k1 ) * int.abs( incx ) - 1)]);

            i += 1;
        }

        double* a_t = NULL;
        // Check leading dimension(s)
        if( lda < n ) {
            info = -4;
            return info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (double*)LAPACKE_malloc( sizeof(double) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        /* Transpose input matrices */
        API_SUFFIX(LAPACKE_dge_trans)( order, lda_t, n, a, lda, a_t, lda_t );
        /* Call LAPACK function and adjust info */
        LAPACK_dlaswp( &n, a_t, &lda_t, &k1, &k2, ipiv, &incx );
        info = 0;  /* LAPACK call is ok! */
        /* Transpose output matrices */
        API_SUFFIX(LAPACKE_dge_trans)( LAPACK_COL_MAJOR, lda_t, n, a_t, lda_t, a, lda );
        /* Release memory and exit */
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_dlaswp_work", info );
        }
    } else {
        info = -1;
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_dlaswp_work", info );
    }
    return info;
}

fn k_laswp(comptime A: type, n: isize, a: [*]A, lda: isize, k1: isize, k2: isize, ipiv: [*]const isize, incx: isize) void {
    // Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
    // K1 through K2.
    var ix0: isize = undefined;
    var I1: isize = undefined;
    var I2: isize = undefined;
    var inc: isize = undefined;
    if (incx > 0) {
        ix0 = k1;
        I1 = k1;
        I2 = k2;
        inc = 1;
    } else if (incx < 0) {
        ix0 = k1 + (k1 - k2) * incx;
        I1 = k2;
        I2 = k1;
        inc = -1;
    } else {
        return;
    }

    var temp: A = undefined;
    const n32: isize = int.div(n, 32) * 32;
    if (n32 != 0) {
        var j: usize = 0;
        while (j < n32) {
            const ix: isize = ix0;
            var i: isize = I1;
            while (i != I2) {
                const ip: isize = ipiv[scast(usize, ix)];
                if (ip != i) {
                    var k: isize = j;
                    while (k < j + 32) {
                        temp = a[scast(usize, i + k * lda)];
                        a[scast(usize, i + k * lda)] = a[scast(usize, ip + k * lda)];
                        a[scast(usize, ip + k * lda)] = temp;

                        k += 1;
                    }
                }

                i += inc;
                ix += incx;
            }

            j += 32;
        }
    }

    if (n32 != n) {
        const ix: isize = ix0;
        var i: isize = I1;
        while (i != I2) {
            const ip: isize = ipiv[scast(usize, ix)];
            if (ip != i) {
                var k: isize = n32;
                while (k < n) {
                    temp = a[scast(usize, i + k * lda)];
                    a[scast(usize, i + k * lda)] = a[scast(usize, ip + k * lda)];
                    a[scast(usize, ip + k * lda)] = temp;

                    k += 1;
                }
            }

            i += inc;
            ix += incx;
        }
    }
}
