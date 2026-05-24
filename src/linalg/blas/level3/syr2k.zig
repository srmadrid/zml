const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;
const Uplo = types.Uplo;

/// Performs a symmetric rank-`2k` update.
///
/// The `syr2k` routine performs a rank-`2k` matrix-matrix operation for a
/// symmetric matrix `C` using general matrices `A` and `B`. The operation is
/// defined as:
///
/// ```zig
///     C = alpha * A * B^T + alpha * B * A^T + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^T * B + alpha * B^T * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` symmetric
/// matrix, `A` and `B` are `n`-by-`k` matrices in the first case and `k`-by-`n`
/// matrices in the second case.
///
/// Signature
/// ---------
/// ```zig
/// fn syr2k(order: Order, uplo: Uplo, trans: Transpose, n: i32, k: i32, alpha: Al, a: [*]const A, lda: i32, b: [*]const B, ldb: i32, beta: Be, c: [*]C, ldc: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `C` is used:
/// - If `uplo = upper`, then the upper triangular part of `C` is used.
/// - If `uplo = lower`, then the lower triangular part of `C` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * B^T + alpha * B * A^T + beta * C`.
/// - If `trans = trans`, then `C = alpha * A^T * B + alpha * B^T * A + beta * C`.
///
/// `n` (`i32`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`i32`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`i32`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`        |
/// | `order = row_major` | `ldb * n`           | `ldb * k`        |
///
/// `ldb` (`i32`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `beta`.
///
/// `c` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `ldc * n`.
///
/// `ldc` (`i32`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)` or `max(1, k)`, if `ldb` is less than `max(1, n)` or
/// `max(1, k)`, or if `ldc` is less than `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub fn syr2k(
    order: Layout,
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.syr2k requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.syr2k requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.syr2k requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B))
        @compileError("zml.linalg.blas.syr2k requires b to be a many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.syr2k requires b's child type to numeric, got " ++ @typeName(B));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.syr2k requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.syr2k requires c to be a mutable many-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.syr2k requires c's child type to be numeric, got " ++ @typeName(C));

    comptime if (Al == bool and A == bool and B == bool and Be == bool and C == bool)
        @compileError("zml.linalg.blas.syr2k does not support alpha, a, b, beta and c all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(C))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.syr2k not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and options.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssyr2k(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                } else if (comptime A == f64) {
                    return ci.cblas_dsyr2k(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_csyr2k(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zsyr2k(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return _syr2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

fn _syr2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_syr2k(
            uplo,
            trans,
            n,
            k,
            alpha,
            a,
            lda,
            b,
            ldb,
            beta,
            c,
            ldc,
            ctx,
        );
    } else {
        return k_syr2k(
            uplo.invert(),
            trans.invert(),
            n,
            k,
            alpha,
            b,
            ldb,
            a,
            lda,
            beta,
            c,
            ldc,
            ctx,
        );
    }
}

fn k_syr2k(
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const Be: type = @TypeOf(beta);
    const C: type = types.Child(@TypeOf(c));
    const T1: type = types.Coerce(Al, B);
    const T2: type = types.Coerce(Al, A);
    const T3: type = types.Coerce(A, B);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(B, types.Coerce(Be, C))));

    const nrowa: i32 = if (trans == .no_trans) n else k;

    if (trans == .conj_no_trans or trans == .conj_trans or
        n < 0 or k < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, nrowa) or ldc < int.max(1, n))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or
        ((ops.eq(alpha, 0, ctx) catch unreachable or k == 0) and
            ops.eq(beta, 1, ctx) catch unreachable))
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            if (uplo == .upper) {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = 0;
                        while (i <= j) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = 0;
                        while (i <= j) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            } else {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = j;
                        while (i < n) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = j;
                        while (i < n) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            }

            return;
        }

        if (trans == .no_trans) {
            if (uplo == .upper) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        var i: i32 = 0;
                        while (i <= j) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                        var i: i32 = 0;
                        while (i <= j) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    }

                    var l: i32 = 0;
                    while (l < k) : (l += 1) {
                        if (ops.ne(a[scast(u32, j + l * lda)], 0, ctx) catch unreachable or
                            ops.ne(b[scast(u32, j + l * ldb)], 0, ctx) catch unreachable)
                        {
                            const temp1: T1 = ops.mul( // temp1 = alpha * b[j + l * ldb]
                                alpha,
                                b[scast(u32, j + l * ldb)],
                                ctx,
                            ) catch unreachable;
                            const temp2: T2 = ops.mul( // temp2 = alpha * a[j + l * lda]
                                alpha,
                                a[scast(u32, j + l * lda)],
                                ctx,
                            ) catch unreachable;

                            var i: i32 = 0;
                            while (i <= j) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += a[i + l * lda] * temp1 + b[i + l * ldb] * temp2
                                    &c[scast(u32, i + j * ldc)],
                                    c[scast(u32, i + j * ldc)],
                                    ops.add(
                                        ops.mul(
                                            a[scast(u32, i + l * lda)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            b[scast(u32, i + l * ldb)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    if (ops.eq(beta, 0, ctx) catch unreachable) {
                        var i: i32 = j;
                        while (i < n) : (i += 1) {
                            ops.set( // c[i + j * ldc] = 0
                                &c[scast(u32, i + j * ldc)],
                                0,
                                ctx,
                            ) catch unreachable;
                        }
                    } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                        var i: i32 = j;
                        while (i < n) : (i += 1) {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    }

                    var l: i32 = 0;
                    while (l < k) : (l += 1) {
                        if (ops.ne(a[scast(u32, j + l * lda)], 0, ctx) catch unreachable or
                            ops.ne(b[scast(u32, j + l * ldb)], 0, ctx) catch unreachable)
                        {
                            const temp1: T1 = ops.mul( // temp1 = alpha * b[j + l * ldb]
                                alpha,
                                b[scast(u32, j + l * ldb)],
                                ctx,
                            ) catch unreachable;
                            const temp2: T2 = ops.mul( // temp2 = alpha * a[j + l * lda]
                                alpha,
                                a[scast(u32, j + l * lda)],
                                ctx,
                            ) catch unreachable;

                            var i: i32 = j;
                            while (i < n) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += a[i + l * lda] * temp1 + b[i + l * ldb] * temp2
                                    &c[scast(u32, i + j * ldc)],
                                    c[scast(u32, i + j * ldc)],
                                    ops.add(
                                        ops.mul(
                                            a[scast(u32, i + l * lda)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            b[scast(u32, i + l * ldb)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                }
            }
        } else {
            if (uplo == .upper) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i <= j) : (i += 1) {
                        var temp1: T3 = constants.zero(T3, ctx) catch unreachable;
                        var temp2: T3 = constants.zero(T3, ctx) catch unreachable;

                        var l: i32 = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp1 += a[l + i * lda] * b[l + j * ldb]
                                &temp1,
                                temp1,
                                ops.mul(
                                    a[scast(u32, l + i * lda)],
                                    b[scast(u32, l + j * ldb)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[l + i * ldb] * a[l + j * lda]
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(u32, l + i * ldb)],
                                    a[scast(u32, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.set( // c[i + j * ldc] = alpha * temp1 + alpha * temp2
                                &c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        alpha,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        alpha,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += alpha * temp1 + alpha * temp2
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        alpha,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        alpha,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            } else {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = j;
                    while (i < n) : (i += 1) {
                        var temp1: T3 = constants.zero(T3, ctx) catch unreachable;
                        var temp2: T3 = constants.zero(T3, ctx) catch unreachable;

                        var l: i32 = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp1 += a[l + i * lda] * b[l + j * ldb]
                                &temp1,
                                temp1,
                                ops.mul(
                                    a[scast(u32, l + i * lda)],
                                    b[scast(u32, l + j * ldb)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[l + i * ldb] * a[l + j * lda]
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(u32, l + i * ldb)],
                                    a[scast(u32, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.set( // c[i + j * ldc] = alpha * temp1 + alpha * temp2
                                &c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        alpha,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        alpha,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.mul_( // c[i + j * ldc] *= beta
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += alpha * temp1 + alpha * temp2
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        alpha,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        alpha,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.syr2k not implemented for arbitrary precision types yet");
    }

    return;
}
