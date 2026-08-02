const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

pub fn sytf2(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    a: anytype,
    lda: usize,
    ipiv: [*]usize,
) !usize {
    const A: type = @TypeOf(a);

    comptime if (!meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.lapack.sytf2: a must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n");

    if (lda < int.max(1, n))
        return linalg.lapack.Error.InvalidArgument;

    var info = int.highest(usize);

    // Initialize alpha for use in choosing the pivot block size.
    const alpha = numeric.div(
        numeric.add(
            1,
            numeric.sqrt(numeric.cast(meta.Numeric(A), 17)),
        ),
        8,
    );

    if (uplo == .upper) {
        // Factorize A as  U * D * Uᵀ.
        var k: usize = n;
        var kstep: usize = 1;
        while (k > 0) {
            k -|= kstep;
            kstep = 1;

            // Determine rows and columns to be interchanged and whether a 1 × 1
            // or 2 × 2 pivot block will be used.
            const absakk = numeric.abs1(a[utils.index(layout, k, k, lda)]);

            // imax is the row index of the largest off-diagonal element in
            // column k, and colmax is its absolute value.
            const imax: usize = if (k > 0)
                linalg.blas.iamax(
                    k,
                    a + utils.index(layout, 0, k, lda),
                    utils.col_ld(layout, lda),
                ) catch unreachable
            else
                0;

            const colmax = if (k > 0)
                numeric.abs1(a[utils.index(layout, imax, k, lda)])
            else
                numeric.zero(numeric.Abs1(meta.Numeric(A)));

            var kp: usize = undefined;
            if (numeric.eq(numeric.max(absakk, colmax), 0) or numeric.isNan(absakk)) {
                // Column k is zero or underflow, or contains a nan, set info and continue.
                if (info == 0)
                    info = k;

                kp = k;
            } else {
                if (numeric.ge(absakk, numeric.mul(alpha, colmax))) {
                    // No interchange, use 1 × 1 pivot block.
                    kp = k;
                } else {
                    // jmax is the column-index of the largest off-diagonal element
                    // in row imax, and rowmax is its absolute value.
                    var jmax = imax + 1 + linalg.blas.iamax(
                        k - imax,
                        a + utils.index(layout, imax, imax + 1, lda),
                        utils.row_ld(layout, lda),
                    );

                    var rowmax = numeric.abs1(a[utils.index(layout, imax, jmax, lda)]);

                    if (imax > 0) {
                        jmax = linalg.blas.iamax(
                            imax,
                            a + utils.index(layout, 0, imax, lda),
                            utils.col_ld(layout, lda),
                        );

                        numeric.maxInto(
                            &rowmax,
                            rowmax,
                            numeric.abs1(a[utils.index(layout, jmax, imax, lda)]),
                        );
                    }

                    if (numeric.ge(absakk, numeric.mul(alpha, numeric.mul(colmax, numeric.div(colmax, rowmax))))) {
                        // No interchange, use 1 × 1 pivot block
                        kp = k;
                    } else if (numeric.ge(numeric.abs1(a[utils.index(layout, imax, imax, lda)]), numeric.mul(alpha, rowmax))) {
                        // Interchange rows and columns k and imax, use 1 × 1 pivot block.
                        kp = imax;
                    } else {
                        // Interchange rows and columns k - 1 and imax, use 2 × 2 pivot block.
                        kp = imax;
                        kstep = 2;
                    }
                }

                const kk = k - kstep + 1;
                if (kp != kk) {
                    // Interchange rows and columns kk and kp in the leading submatrix a[1..k, 1..k].
                    linalg.blas.swap(
                        kp,
                        a + utils.index(layout, 0, kk, lda),
                        utils.col_ld(layout, lda),
                        a + utils.index(layout, 0, kp, lda),
                        utils.col_ld(layout, lda),
                    ) catch unreachable;

                    linalg.blas.swap(
                        kk - kp - 1,
                        a + utils.index(layout, kp + 1, kk, lda),
                        utils.col_ld(layout, lda),
                        a + utils.index(layout, kp, kp + 1, lda),
                        utils.row_ld(layout, lda),
                    ) catch unreachable;

                    var temp = a[utils.index(layout, kk, kk, lda)];
                    a[utils.index(layout, kk, kk, lda)] = a[utils.index(layout, kp, kp, lda)];
                    a[utils.index(layout, kp, kp, lda)] = temp;

                    if (kstep == 2) {
                        temp = a[utils.index(layout, k - 1, k, lda)];
                        a[utils.index(layout, k - 1, k, lda)] = a[utils.index(layout, kp, k, lda)];
                        a[utils.index(layout, kp, k, lda)] = temp;
                    }
                }

                // Update the leading submatrix.
                if (kstep == 1) {
                    // 1-by-1 pivot block d(k): column k now holds
                    // w(k) = u(k)*d(k)
                    // where u(k) is the k-th column of u
                    // perform a rank-1 update of a(1:k-1,1:k-1) as
                    // a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    const r1: A = try ops.div(1, a[utils.index(order, k, k, lda)], ctx);

                    try blas.syr(
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
                        var d12: A = a[utils.index(order, k - 1, k, lda)];
                        const d22: A = try ops.div(
                            a[utils.index(order, k - 1, k - 1, lda)],
                            d12,
                            ctx,
                        );
                        const d11: A = try ops.div(
                            a[utils.index(order, k, k, lda)],
                            d12,
                            ctx,
                        );
                        const t: A = try ops.div(
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
                        d12 = try ops.div(t, d12, ctx);

                        var j: i32 = k - 2;
                        while (j >= 0) : (j -= 1) {
                            const wkm1: A = try ops.mul( // wkm1 = d12 * (d11 * a[j + (k - 1) * lda] - a[j + k * lda])
                                d12,
                                try ops.sub(
                                    try ops.mul(
                                        d11,
                                        a[utils.index(order, j, k - 1, lda)],
                                        ctx,
                                    ),
                                    a[utils.index(order, j, k, lda)],
                                    ctx,
                                ),
                                ctx,
                            );
                            const wk: A = try ops.mul( // wk = d12 * (d22 * a[j + k * lda] - a[j + (k - 1) * lda])
                                d12,
                                try ops.sub(
                                    try ops.mul(
                                        d22,
                                        a[utils.index(order, j, k, lda)],
                                        ctx,
                                    ),
                                    a[utils.index(order, j, k - 1, lda)],
                                    ctx,
                                ),
                                ctx,
                            );

                            var i: i32 = j;
                            while (i >= 0) : (i -= 1) {
                                try ops.sub_( // a[i + j * lda] -= a[i + k * lda] * wk + a[i + (k - 1) * lda] * wkm1
                                    &a[utils.index(order, i, j, lda)],
                                    a[utils.index(order, i, j, lda)],
                                    try ops.add(
                                        try ops.mul(
                                            a[utils.index(order, i, k, lda)],
                                            wk,
                                            ctx,
                                        ),
                                        try ops.mul(
                                            a[utils.index(order, i, k - 1, lda)],
                                            wkm1,
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
        }
    } else {
        // Factorize a as l*d*l**t using the lower triangle of a
        // k is the main loop index, increasing from 1 to n in steps of 1 or 2
        var k: i32 = 0;
        while (k < n) {
            var kstep: i32 = 1;

            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            const absakk: types.Scalar(A) = try ops.abs1(a[utils.index(order, k, k, lda)], ctx);

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
                    } else if (try ops.ge(try ops.abs1(a[utils.index(order, imax, imax, lda)], ctx), try ops.mul(alpha, rowmax, ctx), ctx)) {
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

                    try blas.swap(
                        kp - kk - 1,
                        a + utils.index(order, kk + 1, kk, lda),
                        utils.col_ld(order, lda),
                        a + utils.index(order, kp, kk + 1, lda),
                        utils.row_ld(order, lda),
                        ctx,
                    );

                    var t: A = a[utils.index(order, kk, kk, lda)];
                    a[utils.index(order, kk, kk, lda)] = a[utils.index(order, kp, kp, lda)];
                    a[utils.index(order, kp, kp, lda)] = t;
                    if (kstep == 2) {
                        t = a[utils.index(order, k + 1, k, lda)];
                        a[utils.index(order, k + 1, k, lda)] = a[utils.index(order, kp, k, lda)];
                        a[utils.index(order, kp, k, lda)] = t;
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
                        const d11: A = try ops.div(1, a[utils.index(order, k, k, lda)], ctx);

                        try blas.syr(
                            order,
                            uplo,
                            n - k - 1,
                            try ops.neg(d11, ctx),
                            a + utils.index(order, k + 1, k, lda),
                            utils.col_ld(order, lda),
                            a + utils.index(order, k + 1, k + 1, lda),
                            lda,
                            ctx,
                        );

                        // Store l(k) in column k
                        try blas.scal(
                            n - k - 1,
                            d11,
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
                        var d21: A = a[utils.index(order, k + 1, k, lda)];
                        const d11: A = try ops.div(a[utils.index(order, k + 1, k + 1, lda)], d21, ctx);
                        const d22: A = try ops.div(a[utils.index(order, k, k, lda)], d21, ctx);
                        const t: A = try ops.div(
                            1,
                            try ops.sub(
                                try ops.mul(d11, d22, ctx),
                                1,
                                ctx,
                            ),
                            ctx,
                        );
                        d21 = try ops.div(t, d21, ctx);

                        var j: i32 = k + 2;
                        while (j < n) : (j += 1) {
                            const wk: A = try ops.mul( // wk = d21 * (d11 * a[j + k * lda] - a[j + (k + 1) * lda])
                                d21,
                                try ops.sub(
                                    try ops.mul(
                                        d11,
                                        a[utils.index(order, j, k, lda)],
                                        ctx,
                                    ),
                                    a[utils.index(order, j, k + 1, lda)],
                                    ctx,
                                ),
                                ctx,
                            );
                            const wkp1: A = try ops.mul( // wkp1 = d21 * (d22 * a[j + (k + 1) * lda] - a[j + k * lda])
                                d21,
                                try ops.sub(
                                    try ops.mul(
                                        d22,
                                        a[utils.index(order, j, k + 1, lda)],
                                        ctx,
                                    ),
                                    a[utils.index(order, j, k, lda)],
                                    ctx,
                                ),
                                ctx,
                            );

                            var i: i32 = j;
                            while (i < n) : (i += 1) {
                                try ops.sub_( // a[i + j * lda] -= a[i + k * lda] * wk + a[i + (k + 1) * lda] * wkp1
                                    &a[utils.index(order, i, j, lda)],
                                    a[utils.index(order, i, j, lda)],
                                    try ops.add(
                                        try ops.mul(
                                            a[utils.index(order, i, k, lda)],
                                            wk,
                                            ctx,
                                        ),
                                        try ops.mul(
                                            a[utils.index(order, i, k + 1, lda)],
                                            wkp1,
                                            ctx,
                                        ),
                                        ctx,
                                    ),
                                    ctx,
                                );
                            }

                            a[utils.index(order, j, k, lda)] = wk;
                            a[utils.index(order, j, k + 1, lda)] = wkp1;
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

    return info;
}
