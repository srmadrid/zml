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

/// Performs a symmetric rank-`k` update.
///
/// The `syrk` routine performs a rank-`k` matrix-matrix operation for a
/// symmetric matrix `C` using a general matrix `A`. The operation is defined
/// as:
///
/// ```zig
///     C = alpha * A * A^T + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^T * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` symmetric matrix,
/// `A` is an `n`-by-`k` matrix in the first case and a `k`-by-`n` matrix in the
/// second case.
///
/// Signature
/// ---------
/// ```zig
/// fn syrk(order: Order, uplo: Uplo, trans: Transpose, n: i32, k: i32, alpha: Al, a: [*]const A, lda: i32, beta: Be, c: [*]C, ldc: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * A^T + beta * C`.
/// - If `trans = trans`, then `C = alpha * A^T * A + beta * C`.
///
/// `n` (`i32`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
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
/// than `max(1, k)` or `max(1, n)`, or if `ldc` is less than `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub fn syrk(
    order: Layout,
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.syrk requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.syrk requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.syrk requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.syrk requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.syrk requires c to be a mutable many-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.syrk requires c's child type to be numeric, got " ++ @typeName(C));

    comptime if (Al == bool and A == bool and Be == bool and C == bool)
        @compileError("zml.linalg.blas.syrk does not support alpha, a, b, beta and c all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(C))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.syrk not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and options.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssyrk(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), scast(A, beta), c, scast(c_int, ldc));
                } else if (comptime A == f64) {
                    return ci.cblas_dsyrk(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), scast(A, beta), c, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_csyrk(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zsyrk(order.toCUInt(), uplo.toCUInt(), trans.toCUInt(), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return _syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, ctx);
}

fn _syrk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_syrk(
            uplo,
            trans,
            n,
            k,
            alpha,
            a,
            lda,
            beta,
            c,
            ldc,
            ctx,
        );
    } else {
        return k_syrk(
            uplo.invert(),
            trans.invert(),
            n,
            k,
            alpha,
            a,
            lda,
            beta,
            c,
            ldc,
            ctx,
        );
    }
}

fn k_syrk(
    uplo: Uplo,
    trans: Transpose,
    n: i32,
    k: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    beta: anytype,
    c: anytype,
    ldc: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const Be: type = @TypeOf(beta);
    const C: type = types.Child(@TypeOf(c));
    const T1: type = types.Coerce(Al, A);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(Be, C)));

    const nrowa: i32 = if (trans == .no_trans) n else k;

    if (trans == .conj_no_trans or trans == .conj_trans or
        n < 0 or k < 0 or lda < int.max(1, nrowa) or ldc < int.max(1, n))
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
                        if (ops.ne(a[scast(u32, j + l * lda)], 0, ctx) catch unreachable) {
                            const temp: T1 = ops.mul( // temp = alpha * a[j + l * lda]
                                alpha,
                                a[scast(u32, j + l * lda)],
                                ctx,
                            ) catch unreachable;

                            var i: i32 = 0;
                            while (i <= j) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                    &c[scast(u32, i + j * ldc)],
                                    c[scast(u32, i + j * ldc)],
                                    ops.mul(
                                        temp,
                                        a[scast(u32, i + l * lda)],
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
                        if (ops.ne(a[scast(u32, j + l * lda)], 0, ctx) catch unreachable) {
                            const temp: T1 = ops.mul( // temp = alpha * a[j + l * lda]
                                alpha,
                                a[scast(u32, j + l * lda)],
                                ctx,
                            ) catch unreachable;

                            var i: i32 = j;
                            while (i < n) : (i += 1) {
                                ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                    &c[scast(u32, i + j * ldc)],
                                    c[scast(u32, i + j * ldc)],
                                    ops.mul(
                                        temp,
                                        a[scast(u32, i + l * lda)],
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
                        var temp: A = constants.zero(A, ctx) catch unreachable;

                        var l: i32 = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp += a[l + i * lda] * a[l + j * lda]
                                &temp,
                                temp,
                                ops.mul(
                                    a[scast(u32, l + i * lda)],
                                    a[scast(u32, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.set( // c[i + j * ldc] = alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    alpha,
                                    temp,
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

                            ops.add_( // c[i + j * ldc] += alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    alpha,
                                    temp,
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
                        var temp: A = constants.zero(A, ctx) catch unreachable;

                        var l: i32 = 0;
                        while (l < k) : (l += 1) {
                            ops.add_( // temp += a[l + i * lda] * a[l + j * lda]
                                &temp,
                                temp,
                                ops.mul(
                                    a[scast(u32, l + i * lda)],
                                    a[scast(u32, l + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.set( // c[i + j * ldc] = alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    alpha,
                                    temp,
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

                            ops.add_( // c[i + j * ldc] += alpha * temp
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    alpha,
                                    temp,
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
        @compileError("zml.linalg.blas.syrk not implemented for arbitrary precision types yet");
    }

    return;
}
