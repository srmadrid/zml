const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Uplo = types.Uplo;

const utils = @import("../utils.zig");

pub fn lahef(
    order: Order,
    uplo: Uplo,
    n: i32,
    nb: i32,
    kb: *i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    w: anytype,
    ldw: i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));
    const W: type = types.Child(@TypeOf(w));

    var info: i32 = 0;

    if (comptime !types.isArbitraryPrecision(A)) {
        // Initialize alpha for use inchoosing the pivot block size.
        const alpha: types.Scalar(A) = (1 + try ops.sqrt(types.scast(types.Scalar(A), 17), ctx)) / 8;

        if (uplo == .upper) {
            // Factorize the trailing columns of a using the upper triangle of a and
            // working backwards, and compute the matrix w = u12*d for use in
            // updating a.

            // k is the main loop index, decreasing from n in steps of 1 or 2.
            // kw is the column of w which corresponds to column k of a.
            var k: i32 = n - 1;
            var kw: i32 = undefined;
            while (true) {
                kw = nb + k - n;

                // Exit from loop
                if ((k <= n - nb and nb < n) or k < 0)
                    break;

                // Copy column k of a to column kw of w and update it
                try blas.copy(
                    k,
                    a + utils.index(order, 0, k, lda),
                    utils.col_ld(order, lda),
                    w + utils.index(order, 0, kw, ldw),
                    utils.col_ld(order, ldw),
                    ctx,
                );

                try ops.set(
                    &w[utils.index(order, k, kw, ldw)],
                    try ops.re(a[utils.index(order, k, k, lda)], ctx),
                    ctx,
                );

                if (k < n - 1) {
                    try blas.gemv(
                        order,
                        .no_trans,
                        k + 1,
                        n - (k + 1),
                        -1,
                        a + utils.index(order, 0, k + 1, lda),
                        lda,
                        w + utils.index(order, k, kw + 1, ldw),
                        utils.row_ld(order, ldw),
                        1,
                        w + utils.index(order, 0, kw, ldw),
                        utils.col_ld(order, ldw),
                        ctx,
                    );

                    try ops.set(
                        &w[utils.index(order, k, kw, ldw)],
                        try ops.re(w[utils.index(order, k, kw, ldw)], ctx),
                        ctx,
                    );
                }

                var kstep: i32 = 1;

                // Determine rows and columns to be interchanged and whether a
                // 1-by-1 or 2-by-2 pivot block will be used.
                const absakk: types.Scalar(W) = ops.abs1(try ops.re(w[utils.index(order, k, kw, ldw)], ctx), ctx) catch unreachable;

                // imax is the row-index of the largest off-diagonal element in
                // column k, and colmax is its absolute value.
                // determine both colmax and imax.
                var imax: i32 = undefined;
                var colmax: types.Scalar(W) = undefined;
                if (k > 0) {
                    imax = types.scast(i32, try blas.iamax(
                        k,
                        w + utils.index(order, 0, kw, ldw),
                        utils.col_ld(order, ldw),
                        ctx,
                    ));

                    colmax = ops.abs1(w[utils.index(order, imax, kw, ldw)], ctx) catch unreachable;
                } else {
                    colmax = constants.zero(types.Scalar(W), ctx) catch unreachable;
                }

                var kp: i32 = undefined;
                if (try ops.eq(try ops.max(absakk, colmax, ctx), 0, ctx)) {
                    // Column k is zero or underflow: set info and continue
                    if (info == 0)
                        info = k + 1;

                    kp = k;

                    try ops.set(
                        &a[utils.index(order, k, k, lda)],
                        try ops.re(a[utils.index(order, k, k, lda)], ctx),
                        ctx,
                    );
                } else {
                    if (try ops.ge(absakk, try ops.mul(alpha, colmax, ctx), ctx)) {
                        // No interchange, use 1-by-1 pivot block
                        kp = k;
                    } else {
                        // Copy column imax to column kw-1 of w and update it
                        try blas.copy(
                            imax,
                            a + utils.index(order, 0, imax, lda),
                            utils.col_ld(order, lda),
                            w + utils.index(order, 0, kw - 1, ldw),
                            utils.col_ld(order, ldw),
                            ctx,
                        );

                        try ops.set(
                            &w[utils.index(order, imax, kw - 1, ldw)],
                            try ops.re(a[utils.index(order, imax, imax, lda)], ctx),
                            ctx,
                        );

                        try blas.copy(
                            k - imax,
                            a + utils.index(order, imax, imax + 1, lda),
                            utils.row_ld(order, lda),
                            w + utils.index(order, imax + 1, kw - 1, ldw),
                            utils.col_ld(order, ldw),
                            ctx,
                        );

                        try lapack.lacgv(
                            k - imax,
                            w + utils.index(order, imax + 1, kw - 1, ldw),
                            utils.col_ld(order, ldw),
                        );

                        if (k < n - 1) {
                            try blas.gemv(
                                order,
                                .no_trans,
                                k + 1,
                                n - (k + 1),
                                -1,
                                a + utils.index(order, 0, k + 1, lda),
                                lda,
                                w + utils.index(order, imax, kw + 1, ldw),
                                utils.row_ld(order, ldw),
                                1,
                                w + utils.index(order, 0, kw - 1, ldw),
                                utils.col_ld(order, ldw),
                                ctx,
                            );

                            try ops.set(
                                &w[utils.index(order, imax, kw - 1, ldw)],
                                try ops.re(w[utils.index(order, imax, kw - 1, ldw)], ctx),
                                ctx,
                            );
                        }

                        // jmax is the column-index of the largest off-diagonal
                        // element in row imax, and rowmax is its absolute value
                        var jmax: i32 = imax + 1 + types.scast(i32, try blas.iamax(
                            k - imax,
                            w + utils.index(order, imax + 1, kw - 1, ldw),
                            utils.col_ld(order, ldw),
                            ctx,
                        ));
                        var rowmax: types.Scalar(W) = ops.abs1(w[utils.index(order, jmax, kw - 1, ldw)], ctx) catch unreachable;
                        if (imax > 0) {
                            jmax = types.scast(i32, try blas.iamax(
                                imax,
                                w + utils.index(order, 0, kw - 1, ldw),
                                utils.col_ld(order, ldw),
                                ctx,
                            ));
                            rowmax = try ops.max(rowmax, try ops.abs1(w[utils.index(order, jmax, kw - 1, ldw)], ctx), ctx);
                        }

                        if (try ops.ge(absakk, try ops.mul(alpha, try ops.mul(colmax, try ops.div(colmax, rowmax, ctx), ctx), ctx), ctx)) {
                            // No interchange, use 1-by-1 pivot block
                            kp = k;
                        } else if (try ops.ge(try ops.abs1(try ops.re(w[utils.index(order, imax, kw - 1, ldw)], ctx), ctx), try ops.mul(alpha, rowmax, ctx), ctx)) {
                            // Interchange rows and columns k and imax, use
                            // 1-by-1 pivot block
                            kp = imax;

                            // Copy column kw-1 of w to column kw of w
                            try blas.copy(
                                k + 1,
                                w + utils.index(order, 0, kw - 1, ldw),
                                utils.col_ld(order, ldw),
                                w + utils.index(order, 0, kw, ldw),
                                utils.col_ld(order, ldw),
                                ctx,
                            );
                        } else {
                            // Interchange rows and columns k-1 and imax, use
                            // 2-by-2 pivot block
                            kp = imax;
                            kstep = 2;
                        }
                    }
                }

                // kk is the column of a where pivoting step stopped
                const kk: i32 = k - kstep + 1;

                // kkw is the column of w which corresponds to column kk of a
                const kkw: i32 = nb + kk - n;

                // Interchange rows and columns kp and kk.
                // Updated column kp is already stored in column kkw of w.
                if (kp != kk) {
                    // Copy non-updated column kk to column kp of submatrix a
                    // at step k. no need to copy element into column k (or k
                    // and k-1 for 2-by-2 pivot) of a, since these columns will
                    // be later overwritten.

                    try ops.set(
                        &a[utils.index(order, kp, kp, lda)],
                        try ops.re(a[utils.index(order, kk, kk, lda)], ctx),
                        ctx,
                    );

                    try blas.copy(
                        kk - kp - 1,
                        a + utils.index(order, kp + 1, kk, lda),
                        utils.col_ld(order, lda),
                        a + utils.index(order, kp, kp + 1, lda),
                        utils.row_ld(order, lda),
                        ctx,
                    );

                    try lapack.lacgv(
                        kk - kp - 1,
                        a + utils.index(order, kp, kp + 1, lda),
                        utils.row_ld(order, lda),
                    );

                    if (kp > 0) {
                        try blas.copy(
                            kp,
                            a + utils.index(order, 0, kk, lda),
                            utils.col_ld(order, lda),
                            a + utils.index(order, 0, kp, lda),
                            utils.col_ld(order, lda),
                            ctx,
                        );
                    }

                    // Interchange rows kk and kp in last k+1 to n columns of a
                    // (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    // later overwritten). interchange rows kk and kp in last
                    // kkw to nb columns of w.
                    if (k < n - 1) {
                        try blas.swap(
                            n - (k + 1),
                            a + utils.index(order, kk, k + 1, lda),
                            utils.row_ld(order, lda),
                            a + utils.index(order, kp, k + 1, lda),
                            utils.row_ld(order, lda),
                            ctx,
                        );
                    }

                    try blas.swap(
                        n - kk,
                        w + utils.index(order, kk, kkw, ldw),
                        utils.row_ld(order, ldw),
                        w + utils.index(order, kp, kkw, ldw),
                        utils.row_ld(order, ldw),
                        ctx,
                    );
                }

                if (kstep == 1) {
                    // 1-by-1 pivot block d(k): column kw of w now holds
                    // w(kw) = u(k)*d(k),
                    // where u(k) is the k-th column of u

                    // Store subdiagonal elements of column u(k) and 1-by-1
                    // block d(k) in column k of a. note: diagonal element
                    // u(k,k) is a unit element and not stored.
                    //   a(k,k) := d(k,k) = w(k,kw)
                    //   a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    try blas.copy(
                        k + 1,
                        w + utils.index(order, 0, kw, ldw),
                        utils.col_ld(order, ldw),
                        a + utils.index(order, 0, k, lda),
                        utils.col_ld(order, lda),
                        ctx,
                    );

                    if (k > 0) {
                        const r1: types.Scalar(A) = try ops.div(1, try ops.re(a[utils.index(order, k, k, lda)], ctx), ctx);
                        try blas.scal(
                            k,
                            r1,
                            a + utils.index(order, 0, k, lda),
                            utils.col_ld(order, lda),
                            ctx,
                        );

                        try lapack.lacgv(
                            k,
                            w + utils.index(order, 0, kw, ldw),
                            utils.col_ld(order, ldw),
                        );
                    }
                } else {
                    // 2-by-2 pivot block d(k): columns kw and kw-1 of w now
                    // hold
                    // ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    // where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    // of u

                    // Store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2 block
                    // d(k-1:k,k-1:k) in columns k-1 and k of a. Note: 2-by-2
                    // diagonal block u(k-1:k,k-1:k) is a unit block and not
                    // stored.
                    //   a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                    //   a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                    //   = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 1) {
                        // Compose the columns of the inverse of 2-by-2 pivot
                        // block d in the following way to reduce the number
                        // of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                        // this inverse

                        // d**(-1) = ( d11 d21 )**(-1) =
                        //           ( d21 d22 )
                        //
                        // = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                        //                        ( (-d21 ) ( d11 ) )
                        //
                        // = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                        //
                        //   * ( ( d22/d21 ) (      -1 ) ) =
                        //     ( (      -1 ) ( d11/d21 ) )
                        //
                        // = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                        //                           ( ( -1  ) ( d22 ) )
                        //
                        // = 1/d21 * t * ( ( d11 ) (  -1 ) )
                        //               ( (  -1 ) ( d22 ) )
                        //
                        // = d21 * ( ( d11 ) (  -1 ) )
                        //         ( (  -1 ) ( d22 ) )
                        var d21: W = w[utils.index(order, k - 1, kw, ldw)];
                        const d11: W = try ops.div(
                            w[utils.index(order, k, kw, ldw)],
                            try ops.conjugate(d21, ctx),
                            ctx,
                        );
                        const d22: W = try ops.div(
                            w[utils.index(order, k - 1, kw - 1, ldw)],
                            d21,
                            ctx,
                        );
                        const t: types.Scalar(W) = try ops.div(
                            1,
                            try ops.sub(
                                try ops.re(try ops.mul(
                                    d11,
                                    d22,
                                    ctx,
                                ), ctx),
                                1,
                                ctx,
                            ),
                            ctx,
                        );
                        d21 = try ops.div(t, d21, ctx);

                        // Apdate elements in columns a(k-1) and a(k) as dot
                        // products of rows of ( w(kw-1) w(kw) ) and columns
                        // of d**(-1)
                        var j: i32 = 0;
                        while (j < k - 1) : (j += 1) {
                            try ops.mul_( // a[j + (k - 1) * lda] = d21 * (d11 * w[j + (kw - 1) * ldw] - w[j + kw * ldw]);
                                &a[utils.index(order, j, k - 1, lda)],
                                d21,
                                try ops.sub(
                                    try ops.mul(
                                        d11,
                                        w[utils.index(order, j, kw - 1, ldw)],
                                        ctx,
                                    ),
                                    w[utils.index(order, j, kw, ldw)],
                                    ctx,
                                ),
                                ctx,
                            );

                            try ops.mul_( // a[j + k * lda] = conj(d21) * (d22 * w[j + kw * ldw] - w[j + (kw - 1) * ldw]);
                                &a[utils.index(order, j, k, lda)],
                                try ops.conjugate(d21, ctx),
                                try ops.sub(
                                    try ops.mul(
                                        d22,
                                        w[utils.index(order, j, kw, ldw)],
                                        ctx,
                                    ),
                                    w[utils.index(order, j, kw - 1, ldw)],
                                    ctx,
                                ),
                                ctx,
                            );
                        }
                    }

                    // Copy d(k) to a
                    try ops.set( // a[k-1 + (k-1)*lda] = w[k-1 + (kw-1)*ldw];
                        &a[utils.index(order, k - 1, k - 1, lda)],
                        w[utils.index(order, k - 1, kw - 1, ldw)],
                        ctx,
                    );
                    try ops.set( // a[k-1 + k*lda] = w[k-1 + kw*ldw];
                        &a[utils.index(order, k - 1, k, lda)],
                        w[utils.index(order, k - 1, kw, ldw)],
                        ctx,
                    );
                    try ops.set( // a[k + k*lda] = w[k + kw*ldw];
                        &a[utils.index(order, k, k, lda)],
                        w[utils.index(order, k, kw, ldw)],
                        ctx,
                    );

                    try lapack.lacgv(
                        k,
                        w + utils.index(order, 0, kw, ldw),
                        utils.col_ld(order, ldw),
                    );

                    try lapack.lacgv(
                        k - 1,
                        w + utils.index(order, 0, kw - 1, ldw),
                        utils.col_ld(order, ldw),
                    );
                }

                // Store details of the interchanges in ipiv
                if (kstep == 1) {
                    ipiv[types.scast(u32, k)] = kp + 1; // convert to 1-based for ipiv
                } else {
                    ipiv[types.scast(u32, k)] = -(kp + 1); // convert to 1-based for ipiv
                    ipiv[types.scast(u32, k - 1)] = -(kp + 1); // convert to 1-based for ipiv
                }

                // Decrease k and return to the start of the main loop
                k -= kstep;
            }

            // Update the upper triangle of a11 (= a(1:k,1:k)) as
            // a11 := a11 - u12*d*u12**t = a11 - u12*w**t

            try blas.gemmtr(
                order,
                .upper,
                .no_trans,
                .trans,
                k + 1,
                n - (k + 1),
                -1,
                a + utils.index(order, 0, k + 1, lda),
                lda,
                w + utils.index(order, 0, kw + 1, ldw),
                ldw,
                1,
                a,
                lda,
                ctx,
            );

            // Put u12 in standard form by partially undoing the interchanges
            // in columns k+1:n looping backwards from k+1 to n
            var j: i32 = k + 1;
            while (j < n) {
                // Undo the interchanges (if any) of rows jj and jp at each step j
                // (here, j is a diagonal index)
                const jj: i32 = j;
                var jp: i32 = ipiv[types.scast(u32, j)];
                if (jp < 0) {
                    jp = -jp;
                    // (here, j is a diagonal index)
                    j += 1;
                }
                // (note: here, j is used to determine row length. length n-j+1
                // of the rows to swap back doesn't include diagonal element)
                j += 1;
                if (jp - 1 != jj and j < n) {
                    try blas.swap(
                        n - j,
                        a + utils.index(order, jp - 1, j, lda),
                        utils.row_ld(order, lda),
                        a + utils.index(order, jj, j, lda),
                        utils.row_ld(order, lda),
                        ctx,
                    );
                }
            }

            // Set kb to the number of columns factorized
            kb.* = n - k - 1;
        } else {
            // Factorize the leading columns of a using the lower triangle
            // of a and working forwards, and compute the matrix w = l21*d
            // for use in updating a22

            // k is the main loop index, increasing from 1 in steps of 1 or 2
            var k: i32 = 0;
            while (true) {
                // Exit from loop
                if ((k >= nb - 1 and nb < n) or k > n - 1)
                    break;

                // Copy column k of a to column k of w and update it
                try ops.set(
                    &w[utils.index(order, k, k, ldw)],
                    try ops.re(a[utils.index(order, k, k, lda)], ctx),
                    ctx,
                );

                if (k < n - 1) {
                    try blas.copy(
                        n - k - 1,
                        a + utils.index(order, k + 1, k, lda),
                        utils.col_ld(order, lda),
                        w + utils.index(order, k + 1, k, ldw),
                        utils.col_ld(order, ldw),
                        ctx,
                    );
                }

                try blas.gemv(
                    order,
                    .no_trans,
                    n - k,
                    k,
                    -1,
                    a + utils.index(order, k, 0, lda),
                    lda,
                    w + utils.index(order, k, 0, ldw),
                    utils.row_ld(order, ldw),
                    1,
                    w + utils.index(order, k, k, ldw),
                    utils.col_ld(order, ldw),
                    ctx,
                );

                try ops.set(
                    &w[utils.index(order, k, k, ldw)],
                    try ops.re(w[utils.index(order, k, k, ldw)], ctx),
                    ctx,
                );

                var kstep: i32 = 1;

                // Determine rows and columns to be interchanged and whether
                // a 1-by-1 or 2-by-2 pivot block will be used
                const absakk: types.Scalar(W) = ops.abs1(try ops.re(w[utils.index(order, k, k, ldw)], ctx), ctx) catch unreachable;

                // imax is the row-index of the largest off-diagonal element in
                // column k, and colmax is its absolute value.
                // determine both colmax and imax.
                var imax: i32 = undefined;
                var colmax: types.Scalar(W) = undefined;
                if (k < n - 1) {
                    imax = k + 1 + types.scast(i32, try blas.iamax(
                        n - k - 1,
                        w + utils.index(order, k + 1, k, ldw),
                        utils.col_ld(order, ldw),
                        ctx,
                    ));
                    colmax = ops.abs1(w[utils.index(order, imax, k, ldw)], ctx) catch unreachable;
                } else {
                    colmax = constants.zero(types.Scalar(W), ctx) catch unreachable;
                }

                var kp: i32 = undefined;
                if (try ops.eq(try ops.max(absakk, colmax, ctx), 0, ctx)) {
                    // Column k is zero or underflow: set info and continue
                    if (info == 0)
                        info = k + 1;

                    kp = k;

                    try ops.set(
                        &a[utils.index(order, k, k, lda)],
                        try ops.re(a[utils.index(order, k, k, lda)], ctx),
                        ctx,
                    );
                } else {
                    if (try ops.ge(absakk, try ops.mul(alpha, colmax, ctx), ctx)) {
                        // No interchange, use 1-by-1 pivot block
                        kp = k;
                    } else {
                        // Copy column imax to column k+1 of w and update it
                        try blas.copy(
                            imax - k,
                            a + utils.index(order, imax, k, lda),
                            utils.row_ld(order, lda),
                            w + utils.index(order, k, k + 1, ldw),
                            utils.col_ld(order, ldw),
                            ctx,
                        );

                        try lapack.lacgv(
                            imax - k,
                            w + utils.index(order, k, k + 1, ldw),
                            utils.col_ld(order, ldw),
                        );

                        try ops.set(
                            &w[utils.index(order, imax, k + 1, ldw)],
                            try ops.re(a[utils.index(order, imax, imax, lda)], ctx),
                            ctx,
                        );

                        if (imax < n - 1) {
                            try blas.copy(
                                n - imax - 1,
                                a + utils.index(order, imax + 1, imax, lda),
                                utils.col_ld(order, lda),
                                w + utils.index(order, imax + 1, k + 1, ldw),
                                utils.col_ld(order, ldw),
                                ctx,
                            );
                        }

                        try blas.gemv(
                            order,
                            .no_trans,
                            n - k,
                            k,
                            -1,
                            a + utils.index(order, k, 0, lda),
                            lda,
                            w + utils.index(order, imax, 0, ldw),
                            utils.row_ld(order, ldw),
                            1,
                            w + utils.index(order, k, k + 1, ldw),
                            utils.col_ld(order, ldw),
                            ctx,
                        );

                        try ops.set(
                            &w[utils.index(order, imax, k + 1, ldw)],
                            try ops.re(w[utils.index(order, imax, k + 1, ldw)], ctx),
                            ctx,
                        );

                        // jmax is the column-index of the largest off-diagonal
                        // element in row imax, and rowmax is its absolute value
                        var jmax: i32 = k + types.scast(i32, try blas.iamax(
                            imax - k,
                            w + utils.index(order, k, k + 1, ldw),
                            utils.col_ld(order, ldw),
                            ctx,
                        ));
                        var rowmax: types.Scalar(W) = ops.abs1(w[utils.index(order, jmax, k + 1, ldw)], ctx) catch unreachable;
                        if (imax < n - 1) {
                            jmax = imax + 1 + types.scast(i32, try blas.iamax(
                                n - imax - 1,
                                w + utils.index(order, imax + 1, k + 1, ldw),
                                utils.col_ld(order, ldw),
                                ctx,
                            ));
                            rowmax = try ops.max(rowmax, try ops.abs1(w[utils.index(order, jmax, k + 1, ldw)], ctx), ctx);
                        }

                        if (try ops.ge(absakk, try ops.mul(alpha, try ops.mul(colmax, try ops.div(colmax, rowmax, ctx), ctx), ctx), ctx)) {
                            // No interchange, use 1-by-1 pivot block
                            kp = k;
                        } else if (try ops.ge(try ops.abs1(try ops.re(w[utils.index(order, imax, k + 1, ldw)], ctx), ctx), try ops.mul(alpha, rowmax, ctx), ctx)) {
                            // Interchange rows and columns k and imax, use 1-by-1
                            // pivot block
                            kp = imax;

                            // Copy column k+1 of w to column k of w
                            try blas.copy(
                                n - k,
                                w + utils.index(order, k, k + 1, ldw),
                                utils.col_ld(order, ldw),
                                w + utils.index(order, k, k, ldw),
                                utils.col_ld(order, ldw),
                                ctx,
                            );
                        } else {
                            // Interchange rows and columns k+1 and imax, use 2-by-2
                            // pivot block
                            kp = imax;
                            kstep = 2;
                        }
                    }

                    // kk is the column of a where pivoting step stopped
                    const kk: i32 = k + kstep - 1;

                    // Interchange rows and columns kp and kk.
                    // updated column kp is already stored in column kk of w.
                    if (kp != kk) {
                        // copy non-updated column kk to column kp of submatrix a
                        // at step k. no need to copy element into column k
                        // (or k and k+1 for 2-by-2 pivot) of a, since these columns
                        // will be later overwritten.
                        try ops.set(
                            &a[utils.index(order, kp, kp, lda)],
                            try ops.re(a[utils.index(order, kk, kk, lda)], ctx),
                            ctx,
                        );

                        try blas.copy(
                            kp - kk - 1,
                            a + utils.index(order, kk + 1, kk, lda),
                            utils.col_ld(order, lda),
                            a + utils.index(order, kp, kk + 1, lda),
                            utils.row_ld(order, lda),
                            ctx,
                        );

                        try lapack.lacgv(
                            kp - kk - 1,
                            a + utils.index(order, kp, kk + 1, lda),
                            utils.row_ld(order, lda),
                        );

                        if (kp < n - 1) {
                            try blas.copy(
                                n - kp - 1,
                                a + utils.index(order, kp + 1, kk, lda),
                                utils.col_ld(order, lda),
                                a + utils.index(order, kp + 1, kp, lda),
                                utils.col_ld(order, lda),
                                ctx,
                            );
                        }

                        // Interchange rows kk and kp in first k-1 columns of a
                        // (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                        // later overwritten). interchange rows kk and kp
                        // in first kk columns of w.
                        if (k > 0) {
                            try blas.swap(
                                k,
                                a + utils.index(order, kk, 0, lda),
                                utils.row_ld(order, lda),
                                a + utils.index(order, kp, 0, lda),
                                utils.row_ld(order, lda),
                                ctx,
                            );
                        }

                        try blas.swap(
                            kk + 1,
                            w + utils.index(order, kk, 0, ldw),
                            utils.row_ld(order, ldw),
                            w + utils.index(order, kp, 0, ldw),
                            utils.row_ld(order, ldw),
                            ctx,
                        );
                    }

                    if (kstep == 1) {
                        // 1-by-1 pivot block d(k): column k of w now holds
                        // w(k) = l(k)*d(k),
                        // where l(k) is the k-th column of l

                        // Store subdiag. elements of column l(k)
                        // and 1-by-1 block d(k) in column k of a.
                        // (note: diagonal element l(k,k) is a unit element
                        // and not stored)
                        //   a(k,k) := d(k,k) = w(k,k)
                        //   a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                        try blas.copy(
                            n - k,
                            w + utils.index(order, k, k, ldw),
                            utils.col_ld(order, ldw),
                            a + utils.index(order, k, k, lda),
                            utils.col_ld(order, lda),
                            ctx,
                        );

                        if (k < n - 1) {
                            const r1: types.Scalar(A) = try ops.div(1, try ops.re(a[utils.index(order, k, k, lda)], ctx), ctx);

                            try blas.scal(
                                n - k - 1,
                                r1,
                                a + utils.index(order, k + 1, k, lda),
                                utils.col_ld(order, lda),
                                ctx,
                            );

                            try lapack.lacgv(
                                n - k - 1,
                                w + utils.index(order, k + 1, k, ldw),
                                utils.col_ld(order, ldw),
                            );
                        }
                    } else {
                        // 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                        // ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                        // where l(k) and l(k+1) are the k-th and (k+1)-th columns of l

                        // Store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                        // block d(k:k+1,k:k+1) in columns k and k+1 of a.
                        // (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                        // block and not stored)
                        //   a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                        //   a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                        //   = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                        if (k < n - 2) {
                            // Compose the columns of the inverse of 2-by-2 pivot
                            // block d in the following way to reduce the number
                            // of flops when we myltiply panel ( w(k) w(k+1) ) by
                            // this inverse

                            // d**(-1) = ( d11 d21 )**(-1) =
                            //           ( d21 d22 )
                            //
                            // = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                            //                        ( (-d21 ) ( d11 ) )
                            //
                            // = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                            //
                            //   * ( ( d22/d21 ) (      -1 ) ) =
                            //     ( (      -1 ) ( d11/d21 ) )
                            //
                            // = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                            //                           ( ( -1  ) ( d22 ) )
                            //
                            // = 1/d21 * t * ( ( d11 ) (  -1 ) )
                            //               ( (  -1 ) ( d22 ) )
                            //
                            // = d21 * ( ( d11 ) (  -1 ) )
                            //         ( (  -1 ) ( d22 ) )
                            var d21: W = w[utils.index(order, k + 1, k, ldw)];
                            const d11: W = try ops.div(
                                w[utils.index(order, k + 1, k + 1, ldw)],
                                d21,
                                ctx,
                            );
                            const d22: W = try ops.div(
                                w[utils.index(order, k, k, ldw)],
                                try ops.conjugate(d21, ctx),
                                ctx,
                            );
                            const t: types.Scalar(W) = try ops.div(
                                1,
                                try ops.sub(
                                    try ops.re(try ops.mul(
                                        d11,
                                        d22,
                                        ctx,
                                    ), ctx),
                                    1,
                                    ctx,
                                ),
                                ctx,
                            );
                            d21 = try ops.div(t, d21, ctx);

                            // Update elements in columns a(k) and a(k+1) as
                            // dot products of rows of ( w(k) w(k+1) ) and columns
                            // of d**(-1)
                            var j: i32 = k + 2;
                            while (j < n) : (j += 1) {
                                try ops.mul_( // a[j + k * lda] = conj(d21) * (d11 * w[j + k * ldw] - w[j + (k + 1) * ldw]);
                                    &a[utils.index(order, j, k, lda)],
                                    try ops.conjugate(d21, ctx),
                                    try ops.sub(
                                        try ops.mul(
                                            d11,
                                            w[utils.index(order, j, k, ldw)],
                                            ctx,
                                        ),
                                        w[utils.index(order, j, k + 1, ldw)],
                                        ctx,
                                    ),
                                    ctx,
                                );

                                try ops.mul_( // a[j + (k + 1) * lda] = d21 * (d22 * w[j + (k + 1) * ldw] - w[j + k * ldw]);
                                    &a[utils.index(order, j, k + 1, lda)],
                                    d21,
                                    try ops.sub(
                                        try ops.mul(
                                            d22,
                                            w[utils.index(order, j, k + 1, ldw)],
                                            ctx,
                                        ),
                                        w[utils.index(order, j, k, ldw)],
                                        ctx,
                                    ),
                                    ctx,
                                );
                            }
                        }

                        // Copy d(k) to a
                        try ops.set( // a[k + k * lda] = w[k + k * ldw];
                            &a[utils.index(order, k, k, lda)],
                            w[utils.index(order, k, k, ldw)],
                            ctx,
                        );
                        try ops.set( // a[k + 1 + k * lda] = w[k + 1 + k * ldw];
                            &a[utils.index(order, k + 1, k, lda)],
                            w[utils.index(order, k + 1, k, ldw)],
                            ctx,
                        );
                        try ops.set( // a[k + 1 + (k + 1) * lda] = w[k + 1 + (k + 1) * ldw];
                            &a[utils.index(order, k + 1, k + 1, lda)],
                            w[utils.index(order, k + 1, k + 1, ldw)],
                            ctx,
                        );

                        try lapack.lacgv(
                            n - k - 1,
                            w + utils.index(order, k + 1, k, ldw),
                            utils.col_ld(order, ldw),
                        );

                        try lapack.lacgv(
                            n - k - 2,
                            w + utils.index(order, k + 2, k + 1, ldw),
                            utils.col_ld(order, ldw),
                        );
                    }
                }

                // store details of the interchanges in ipiv
                if (kstep == 1) {
                    ipiv[types.scast(u32, k)] = kp + 1; // convert to 1-based for ipiv
                } else {
                    ipiv[types.scast(u32, k)] = -(kp + 1); // convert to 1-based for ipiv
                    ipiv[types.scast(u32, k + 1)] = -(kp + 1); // convert to 1-based for ipiv
                }

                // Increase k and return to the start of the main loop
                k += kstep;
            }

            // Update the lower triangle of a22 (= a(k:n,k:n)) as
            // a22 := a22 - l21*d*l21**t = a22 - l21*w**t
            try blas.gemmtr(
                order,
                .lower,
                .no_trans,
                .trans,
                n - k,
                k,
                -1,
                a + utils.index(order, k, 0, lda),
                lda,
                w + utils.index(order, k, 0, ldw),
                ldw,
                1,
                a + utils.index(order, k, k, lda),
                lda,
                ctx,
            );

            // Put l21 in standard form by partially undoing the interchanges
            // of rows in columns 1:k-1 looping backwards from k-1 to 1
            var j: i32 = k - 1;
            while (j > 0) {
                // Undo the interchanges (if any) of rows jj and jp at each step j
                // (here, j is a diagonal index)
                const jj: i32 = j;
                var jp: i32 = ipiv[types.scast(u32, j)];
                if (jp < 0) {
                    jp = -jp;
                    // (here, j is a diagonal index)
                    j -= 1;
                }
                // (note: here, j is used to determine row length. length j
                // of the rows to swap back doesn't include diagonal element)
                j -= 1;
                if (jp - 1 != jj and j >= 0) {
                    try blas.swap(
                        j + 1,
                        a + utils.index(order, jp - 1, 0, lda),
                        utils.row_ld(order, lda),
                        a + utils.index(order, jj, 0, lda),
                        utils.row_ld(order, lda),
                        ctx,
                    );
                }
            }

            // Set kb to the number of columns factorized
            kb.* = k;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.lahef not implemented for arbitrary precision types yet");
    }

    return info;
}
