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

pub fn hetf2(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));

    if (n < 0 or lda < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    if (comptime !types.isArbitraryPrecision(A)) {
        // Initialize alpha for use inchoosing the pivot block size.
        const alpha: types.Scalar(A) = (1 + try ops.sqrt(types.scast(types.Scalar(A), 17), ctx)) / 8;

        if (uplo == .upper) {
            // Factorize a as u*d*u**t using the upper triangle of a
            // k is the main loop index, decreasing from n to 1 in steps of 1 or 2
            var k: i32 = n - 1;
            while (k >= 0) {
                var kstep: i32 = 1;

                // Determine rows and columns to be interchanged and whether
                // a 1-by-1 or 2-by-2 pivot block will be used
                const absakk: types.Scalar(A) = try ops.abs1(try ops.re(a[utils.index(order, k, k, lda)], ctx), ctx);

                // imax is the row-index of the largest off-diagonal element in
                // column k, and colmax is its absolute value.
                // determine both colmax and imax.
                var imax: i32 = 0;
                var colmax: types.Scalar(A) = undefined;
                if (k > 0) {
                    imax = types.scast(i32, try blas.iamax(
                        k,
                        a + utils.index(order, 0, k, lda),
                        utils.col_ld(order, lda),
                        ctx,
                    ));
                    colmax = try ops.abs1(a[utils.index(order, imax, k, lda)], ctx);
                } else {
                    colmax = try constants.zero(types.Scalar(A), ctx);
                }

                var kp: i32 = undefined;
                if (try ops.eq(try ops.max(absakk, colmax, ctx), 0, ctx) or (types.numericType(types.Scalar(A)) == .float and std.math.isNan(absakk))) {
                    // Column k is zero or underflow, or contains a nan:
                    // set info and continue
                    if (info == 0)
                        info = k + 1;

                    kp = k;
                } else {
                    if (try ops.ge(absakk, try ops.mul(alpha, colmax, ctx), ctx)) {
                        // No interchange, use 1-by-1 pivot block
                        kp = k;
                    } else {
                        // jmax is the column-index of the largest off-diagonal
                        // element in row imax, and rowmax is its absolute value
                        var jmax: i32 = imax + 1 + types.scast(i32, try blas.iamax(
                            k - imax,
                            a + utils.index(order, imax, imax + 1, lda),
                            utils.row_ld(order, lda),
                            ctx,
                        ));
                        var rowmax: types.Scalar(A) = try ops.abs1(a[utils.index(order, imax, jmax, lda)], ctx);
                        if (imax > 0) {
                            jmax = types.scast(i32, try blas.iamax(
                                imax,
                                a + utils.index(order, 0, imax, lda),
                                utils.col_ld(order, lda),
                                ctx,
                            ));
                            rowmax = try ops.max(rowmax, try ops.abs1(a[utils.index(order, jmax, imax, lda)], ctx), ctx);
                        }

                        if (try ops.ge(absakk, try ops.mul(alpha, try ops.mul(colmax, try ops.div(colmax, rowmax, ctx), ctx), ctx), ctx)) {
                            // No interchange, use 1-by-1 pivot block
                            kp = k;
                        } else if (try ops.ge(try ops.abs1(a[utils.index(order, imax, imax, lda)], ctx), try ops.mul(alpha, rowmax, ctx), ctx)) {
                            // Interchange rows and columns k and imax, use 1-by-1
                            // pivot block
                            kp = imax;
                        } else {
                            // Interchange rows and columns k-1 and imax, use 2-by-2
                            // pivot block
                            kp = imax;
                            kstep = 2;
                        }
                    }

                    const kk: i32 = k - kstep + 1;
                    if (kp != kk) {
                        // Interchange rows and columns kk and kp in the leading
                        // submatrix a(1:k,1:k)
                        try blas.swap(
                            kp,
                            a + utils.index(order, 0, kk, lda),
                            utils.col_ld(order, lda),
                            a + utils.index(order, 0, kp, lda),
                            utils.col_ld(order, lda),
                            ctx,
                        );

                        var j: i32 = kp + 1;
                        while (j < kk) : (j += 1) {
                            const t: A = try ops.conj(a[utils.index(order, j, kk, lda)], ctx);
                            a[utils.index(order, j, kk, lda)] = try ops.conj(a[utils.index(order, kp, j, lda)], ctx);
                            a[utils.index(order, kp, j, lda)] = t;
                        }

                        try ops.conj_(
                            &a[utils.index(order, kp, kk, lda)],
                            a[utils.index(order, kp, kk, lda)],
                            ctx,
                        );

                        const r1: types.Scalar(A) = try ops.re(a[utils.index(order, kk, kk, lda)], ctx);
                        try ops.set(
                            &a[utils.index(order, kk, kk, lda)],
                            try ops.re(a[utils.index(order, kp, kp, lda)], ctx),
                            ctx,
                        );
                        try ops.set(
                            &a[utils.index(order, kp, kp, lda)],
                            r1,
                            ctx,
                        );

                        if (kstep == 2) {
                            try ops.re_(
                                &a[utils.index(order, k, k, lda)],
                                a[utils.index(order, k, k, lda)],
                                ctx,
                            );
                            const t: A = a[utils.index(order, k - 1, k, lda)];
                            try ops.set(
                                &a[utils.index(order, k - 1, k, lda)],
                                a[utils.index(order, kp, k, lda)],
                                ctx,
                            );
                            try ops.set(
                                &a[utils.index(order, kp, k, lda)],
                                t,
                                ctx,
                            );
                        }
                    } else {
                        try ops.re_(
                            &a[utils.index(order, k, k, lda)],
                            a[utils.index(order, k, k, lda)],
                            ctx,
                        );

                        if (kstep == 2) {
                            try ops.re_(
                                &a[utils.index(order, k - 1, k - 1, lda)],
                                a[utils.index(order, k - 1, k - 1, lda)],
                                ctx,
                            );
                        }
                    }

                    // Update the leading submatrix
                    if (kstep == 1) {
                        // 1-by-1 pivot block d(k): column k now holds
                        // w(k) = u(k)*d(k)
                        // where u(k) is the k-th column of u
                        // perform a rank-1 update of a(1:k-1,1:k-1) as
                        // a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                        const r1: types.Scalar(A) = try ops.div(1, try ops.re(a[utils.index(order, k, k, lda)], ctx), ctx);

                        try blas.her(
                            order,
                            uplo,
                            k,
                            try ops.neg(r1, ctx),
                            a + utils.index(order, 0, k, lda),
                            utils.col_ld(order, lda),
                            a,
                            lda,
                            ctx,
                        );

                        // Store u(k) in column k
                        try blas.scal(
                            k,
                            r1,
                            a + utils.index(order, 0, k, lda),
                            utils.col_ld(order, lda),
                            ctx,
                        );
                    } else {
                        // 2-by-2 pivot block d(k): columns k and k-1 now hold
                        // ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                        // where u(k) and u(k-1) are the k-th and (k-1)-th columns of u
                        // perform a rank-2 update of a(1:k-2,1:k-2) as
                        // a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                        //   = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                        if (k > 1) {
                            var d: types.Scalar(A) = try ops.abs(a[utils.index(order, k - 1, k, lda)], ctx);
                            const d22: types.Scalar(A) = try ops.div(
                                try ops.re(a[utils.index(order, k - 1, k - 1, lda)], ctx),
                                d,
                                ctx,
                            );
                            const d11: types.Scalar(A) = try ops.div(
                                try ops.re(a[utils.index(order, k, k, lda)], ctx),
                                d,
                                ctx,
                            );
                            const t: types.Scalar(A) = try ops.div(
                                1,
                                try ops.sub(
                                    try ops.mul(
                                        d11,
                                        d22,
                                        ctx,
                                    ),
                                    1,
                                    ctx,
                                ),
                                ctx,
                            );
                            const d12: A = try ops.div(
                                a[utils.index(order, k - 1, k, lda)],
                                d,
                                ctx,
                            );
                            d = try ops.div(t, d, ctx);

                            var j: i32 = k - 2;
                            while (j >= 0) : (j -= 1) {
                                const wkm1: A = try ops.mul( // wkm1 = d * (d11 * a[j + (k - 1) * lda] - conj(d12) * a[j + k * lda])
                                    d,
                                    try ops.sub(
                                        try ops.mul(
                                            d11,
                                            a[utils.index(order, j, k - 1, lda)],
                                            ctx,
                                        ),
                                        try ops.mul(
                                            try ops.conj(d12, ctx),
                                            a[utils.index(order, j, k, lda)],
                                            ctx,
                                        ),
                                        ctx,
                                    ),
                                    ctx,
                                );
                                const wk: A = try ops.mul( // wk = d * (d22 * a[j + k * lda] - d12 * a[j + (k - 1) * lda])
                                    d,
                                    try ops.sub(
                                        try ops.mul(
                                            d22,
                                            a[utils.index(order, j, k, lda)],
                                            ctx,
                                        ),
                                        try ops.mul(
                                            d12,
                                            a[utils.index(order, j, k - 1, lda)],
                                            ctx,
                                        ),
                                        ctx,
                                    ),
                                    ctx,
                                );

                                var i: i32 = j;
                                while (i >= 0) : (i -= 1) {
                                    try ops.sub_( // a[i + j * lda] -= a[i + k * lda] * conj(wk) + a[i + (k - 1) * lda] * conj(wkm1)
                                        &a[utils.index(order, i, j, lda)],
                                        a[utils.index(order, i, j, lda)],
                                        try ops.add(
                                            try ops.mul(
                                                a[utils.index(order, i, k, lda)],
                                                try ops.conj(wk, ctx),
                                                ctx,
                                            ),
                                            try ops.mul(
                                                a[utils.index(order, i, k - 1, lda)],
                                                try ops.conj(wkm1, ctx),
                                                ctx,
                                            ),
                                            ctx,
                                        ),
                                        ctx,
                                    );
                                }

                                try ops.set(
                                    &a[utils.index(order, j, k, lda)],
                                    wk,
                                    ctx,
                                );

                                try ops.set(
                                    &a[utils.index(order, j, k - 1, lda)],
                                    wkm1,
                                    ctx,
                                );

                                try ops.set(
                                    &a[utils.index(order, j, j, lda)],
                                    try ops.re(a[utils.index(order, j, j, lda)], ctx),
                                    ctx,
                                );
                            }
                        }
                    }
                }

                // Store details of the interchanges in ipiv
                if (kstep == 1) {
                    ipiv[types.scast(u32, k)] = kp + 1;
                } else {
                    ipiv[types.scast(u32, k)] = -(kp + 1);
                    ipiv[types.scast(u32, k - 1)] = -(kp + 1);
                }

                // Decrease k and return to the start of the main loop
                k -= kstep;
            }
        } else {
            // Factorize a as l*d*l**t using the lower triangle of a
            // k is the main loop index, increasing from 1 to n in steps of 1 or 2
            var k: i32 = 0;
            while (k < n) {
                var kstep: i32 = 1;

                // Determine rows and columns to be interchanged and whether
                // a 1-by-1 or 2-by-2 pivot block will be used
                const absakk: types.Scalar(A) = try ops.abs1(try ops.re(a[utils.index(order, k, k, lda)], ctx), ctx);

                // imax is the row-index of the largest off-diagonal element in
                // column k, and colmax is its absolute value.
                // determine both colmax and imax.
                var imax: i32 = 0;
                var colmax: types.Scalar(A) = undefined;
                if (k < n - 1) {
                    imax = k + 1 + types.scast(i32, try blas.iamax(
                        n - k - 1,
                        a + utils.index(order, k + 1, k, lda),
                        utils.col_ld(order, lda),
                        ctx,
                    ));
                    colmax = try ops.abs1(a[utils.index(order, imax, k, lda)], ctx);
                } else {
                    colmax = try constants.zero(types.Scalar(A), ctx);
                }

                var kp: i32 = undefined;
                if (try ops.eq(try ops.max(absakk, colmax, ctx), 0, ctx) or (types.numericType(types.Scalar(A)) == .float and std.math.isNan(absakk))) {
                    std.debug.print("absakk: {d}, colmax: {d}\n", .{ absakk, colmax });
                    // Column k is zero or underflow, or contains a nan:
                    // set info and continue
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
                        // jmax is the column-index of the largest off-diagonal
                        // element in row imax, and rowmax is its absolute value
                        var jmax: i32 = k + types.scast(i32, try blas.iamax(
                            imax - k,
                            a + utils.index(order, imax, k, lda),
                            utils.row_ld(order, lda),
                            ctx,
                        ));
                        var rowmax: types.Scalar(A) = try ops.abs1(a[utils.index(order, imax, jmax, lda)], ctx);
                        if (imax < n - 1) {
                            jmax = imax + 1 + types.scast(i32, try blas.iamax(
                                n - imax - 1,
                                a + utils.index(order, imax + 1, imax, lda),
                                utils.col_ld(order, lda),
                                ctx,
                            ));
                            rowmax = try ops.max(rowmax, try ops.abs1(a[utils.index(order, jmax, imax, lda)], ctx), ctx);
                        }

                        if (try ops.ge(absakk, try ops.mul(alpha, try ops.mul(colmax, try ops.div(colmax, rowmax, ctx), ctx), ctx), ctx)) {
                            // No interchange, use 1-by-1 pivot block
                            kp = k;
                        } else if (try ops.ge(try ops.abs1(try ops.re(a[utils.index(order, imax, imax, lda)], ctx), ctx), try ops.mul(alpha, rowmax, ctx), ctx)) {
                            // Interchange rows and columns k and imax, use 1-by-1
                            // pivot block
                            kp = imax;
                        } else {
                            // Interchange rows and columns k+1 and imax, use 2-by-2
                            // pivot block
                            kp = imax;
                            kstep = 2;
                        }
                    }

                    const kk: i32 = k + kstep - 1;
                    if (kp != kk) {
                        // interchange rows and columns kk and kp in the trailing
                        // submatrix a(k:n,k:n)
                        if (kp < n - 1) {
                            try blas.swap(
                                n - kp - 1,
                                a + utils.index(order, kp + 1, kk, lda),
                                utils.col_ld(order, lda),
                                a + utils.index(order, kp + 1, kp, lda),
                                utils.col_ld(order, lda),
                                ctx,
                            );
                        }

                        var j: i32 = kk + 1;
                        while (j < kp) : (j += 1) {
                            const t: A = try ops.conj(a[utils.index(order, j, kk, lda)], ctx);
                            a[utils.index(order, j, kk, lda)] = try ops.conj(a[utils.index(order, kp, j, lda)], ctx);
                            a[utils.index(order, kp, j, lda)] = t;
                        }

                        try ops.conj_(
                            &a[utils.index(order, kp, kk, lda)],
                            a[utils.index(order, kp, kk, lda)],
                            ctx,
                        );
                        const r1: types.Scalar(A) = try ops.re(a[utils.index(order, kk, kk, lda)], ctx);
                        try ops.set(
                            &a[utils.index(order, kk, kk, lda)],
                            try ops.re(a[utils.index(order, kp, kp, lda)], ctx),
                            ctx,
                        );
                        try ops.set(
                            &a[utils.index(order, kp, kp, lda)],
                            r1,
                            ctx,
                        );

                        if (kstep == 2) {
                            try ops.re_(
                                &a[utils.index(order, k, k, lda)],
                                a[utils.index(order, k, k, lda)],
                                ctx,
                            );
                            const t: A = a[utils.index(order, k + 1, k, lda)];
                            try ops.set(
                                &a[utils.index(order, k + 1, k, lda)],
                                a[utils.index(order, kp, k, lda)],
                                ctx,
                            );
                            try ops.set(
                                &a[utils.index(order, kp, k, lda)],
                                t,
                                ctx,
                            );
                        }
                    } else {
                        try ops.re_(
                            &a[utils.index(order, k, k, lda)],
                            a[utils.index(order, k, k, lda)],
                            ctx,
                        );

                        if (kstep == 2) {
                            try ops.re_(
                                &a[utils.index(order, k + 1, k + 1, lda)],
                                a[utils.index(order, k + 1, k + 1, lda)],
                                ctx,
                            );
                        }
                    }

                    // Update the trailing submatrix
                    if (kstep == 1) {
                        // 1-by-1 pivot block d(k): column k now holds
                        // w(k) = l(k)*d(k)
                        // where l(k) is the k-th column of l
                        if (k < n - 1) {
                            // Perform a rank-1 update of a(k+1:n,k+1:n) as
                            // a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                            const r1: types.Scalar(A) = try ops.div(1, try ops.re(a[utils.index(order, k, k, lda)], ctx), ctx);

                            try blas.her(
                                order,
                                uplo,
                                n - k - 1,
                                try ops.neg(r1, ctx),
                                a + utils.index(order, k + 1, k, lda),
                                utils.col_ld(order, lda),
                                a + utils.index(order, k + 1, k + 1, lda),
                                lda,
                                ctx,
                            );

                            // Store l(k) in column k
                            try blas.scal(
                                n - k - 1,
                                r1,
                                a + utils.index(order, k + 1, k, lda),
                                utils.col_ld(order, lda),
                                ctx,
                            );
                        }
                    } else {
                        // 2-by-2 pivot block d(k)
                        if (k < n - 2) {
                            // Perform a rank-2 update of a(k+2:n,k+2:n) as
                            // a := a - ( (a(k) a(k+1))*d(k)**(-1) ) * (a(k) a(k+1))**t
                            // where l(k) and l(k+1) are the k-th and (k+1)-th columns of l
                            var d: types.Scalar(A) = try ops.abs(a[utils.index(order, k + 1, k, lda)], ctx);
                            const d11: types.Scalar(A) = try ops.div(
                                try ops.re(a[utils.index(order, k + 1, k + 1, lda)], ctx),
                                d,
                                ctx,
                            );
                            const d22: types.Scalar(A) = try ops.div(
                                try ops.re(a[utils.index(order, k, k, lda)], ctx),
                                d,
                                ctx,
                            );
                            const t: types.Scalar(A) = try ops.div(
                                1,
                                try ops.sub(
                                    try ops.mul(d11, d22, ctx),
                                    1,
                                    ctx,
                                ),
                                ctx,
                            );
                            const d21: A = try ops.div(
                                a[utils.index(order, k + 1, k, lda)],
                                d,
                                ctx,
                            );
                            d = try ops.div(t, d, ctx);

                            var j: i32 = k + 2;
                            while (j < n) : (j += 1) {
                                const wk: A = try ops.mul( // wk = d * (d11 * a[j + k * lda] - d21 * a[j + (k + 1) * lda])
                                    d,
                                    try ops.sub(
                                        try ops.mul(
                                            d11,
                                            a[utils.index(order, j, k, lda)],
                                            ctx,
                                        ),
                                        try ops.mul(
                                            d21,
                                            a[utils.index(order, j, k + 1, lda)],
                                            ctx,
                                        ),
                                        ctx,
                                    ),
                                    ctx,
                                );
                                const wkp1: A = try ops.mul( // wkp1 = d * (d22 * a[j + (k + 1) * lda] - conj(d21) * a[j + k * lda])
                                    d,
                                    try ops.sub(
                                        try ops.mul(
                                            d22,
                                            a[utils.index(order, j, k + 1, lda)],
                                            ctx,
                                        ),
                                        try ops.mul(
                                            try ops.conj(d21, ctx),
                                            a[utils.index(order, j, k, lda)],
                                            ctx,
                                        ),
                                        ctx,
                                    ),
                                    ctx,
                                );

                                var i: i32 = j;
                                while (i < n) : (i += 1) {
                                    try ops.sub_( // a[i + j * lda] -= a[i + k * lda] * conj(wk) + a[i + (k + 1) * lda] * conj(wkp1)
                                        &a[utils.index(order, i, j, lda)],
                                        a[utils.index(order, i, j, lda)],
                                        try ops.add(
                                            try ops.mul(
                                                a[utils.index(order, i, k, lda)],
                                                try ops.conj(wk, ctx),
                                                ctx,
                                            ),
                                            try ops.mul(
                                                a[utils.index(order, i, k + 1, lda)],
                                                try ops.conj(wkp1, ctx),
                                                ctx,
                                            ),
                                            ctx,
                                        ),
                                        ctx,
                                    );
                                }

                                a[utils.index(order, j, k, lda)] = wk;
                                a[utils.index(order, j, k + 1, lda)] = wkp1;
                                try ops.set(
                                    &a[utils.index(order, j, j, lda)],
                                    try ops.re(a[utils.index(order, j, j, lda)], ctx),
                                    ctx,
                                );
                            }
                        }
                    }
                }

                // Store details of the interchanges in ipiv
                if (kstep == 1) {
                    ipiv[types.scast(u32, k)] = kp + 1;
                } else {
                    ipiv[types.scast(u32, k)] = -(kp + 1);
                    ipiv[types.scast(u32, k + 1)] = -(kp + 1);
                }

                // Increase k and return to the start of the main loop
                k += kstep;
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.hetf2 not implemented for arbitrary precision types yet");
    }

    return info;
}
