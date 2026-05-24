const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;
const Side = linalg.Side;
const Uplo = types.Uplo;
const Diag = types.Diag;

/// Solves a triangular matrix equation.
///
/// The `trsm` routine solves one of the following matrix equations:
///
/// ```zig
///     op(A) * X = alpha * B,
/// ```
///
/// or
///
/// ```zig
///     X * op(A) = alpha * B,
/// ```
///
/// where `op(X)` is `X`, `X^T`, `conj(X)`, or `X^H`, `alpha` is a scalar,
/// `A` is a unit, or non-unit, upper or lower triangular matrix, and `B` and
/// `X` are `m`-by-`n` matrices.
///
/// Signature
/// ---------
/// ```zig
/// fn trsm(order: Order, side: Side, uplo: Uplo, transa: Transpose, diag: Diag, m: i32, n: i32, alpha: Al, a: [*]const A, lda: i32, b: [*]B, ldb: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `side` (`Side`): Specifies whether the triangular matrix `A` is on the left
/// or right side of the product:
/// - If `side = left`, then `B = alpha * op(A) * B`.
/// - If `side = right`, then `B = alpha * B * op(A)`.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is upper or lower
/// triangular:
/// - If `uplo = upper`, then `A` is an upper triangular matrix.
/// - If `uplo = lower`, then `A` is a lower triangular matrix.
///
/// `transa` (`Transpose`): Specifies the form of `op(A)`:
/// - If `transa = no_trans`, then `op(A) = A`.
/// - If `transa = trans`, then `op(A) = A^T`.
/// - If `transa = conj_no_trans`, then `op(A) = conj(A)`.
/// - If `transa = conj_trans`, then `op(A) = A^H`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then `A` is a unit triangular matrix.
/// - If `diag = non_unit`, then `A` is not a unit triangular matrix.
///
/// `m` (`i32`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`i32`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * k`, where
/// `k` is `m` if `side = left` and `n` if `side = right`.
///
/// `lda` (`i32`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `ldb * n` if
/// `order = col_major` or `ldb * m` if `order = row_major`.
///
/// `ldb` (`i32`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)` or `max(1, m)`, or if `ldb` is less than `max(1, n)` or
/// `max(1, m)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub fn trsm(
    order: Layout,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: i32,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.trsm requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.trsm requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.trsm requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B) or types.isConstPointer(B))
        @compileError("zml.linalg.blas.trsm requires b to be a mutable many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.trsm requires b's child type to numeric, got " ++ @typeName(B));

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.trsm not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and types.canCoerce(Al, A) and options.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_strsm(order.toCUInt(), side.toCUInt(), uplo.toCUInt(), transa.toCUInt(), diag.toCUInt(), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb));
                } else if (comptime A == f64) {
                    return ci.cblas_dtrsm(order.toCUInt(), side.toCUInt(), uplo.toCUInt(), transa.toCUInt(), diag.toCUInt(), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_ctrsm(order.toCUInt(), side.toCUInt(), uplo.toCUInt(), transa.toCUInt(), diag.toCUInt(), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_ztrsm(order.toCUInt(), side.toCUInt(), uplo.toCUInt(), transa.toCUInt(), diag.toCUInt(), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb));
                }
            },
            else => {},
        }
    }

    return _trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, ctx);
}

fn _trsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: i32,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_trsm(
            side,
            uplo,
            transa,
            diag,
            m,
            n,
            alpha,
            a,
            lda,
            b,
            ldb,
            ctx,
        );
    } else {
        return k_trsm(
            side.invert(),
            uplo.invert(),
            transa,
            diag,
            n,
            m,
            alpha,
            a,
            lda,
            b,
            ldb,
            ctx,
        );
    }
}

fn k_trsm(
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: i32,
    n: i32,
    alpha: anytype,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const T1: type = types.Coerce(Al, types.Coerce(A, B));
    const CC: type = types.Coerce(Al, types.Coerce(A, B));

    const not: bool = transa == .no_trans or transa == .conj_no_trans;
    const noconj: bool = transa == .no_trans or transa == .trans;

    const nrowa: i32 = if (side == .left) m else n;

    if (m < 0 or n < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            var j: i32 = 0;
            while (j < n) : (j += 1) {
                var i: i32 = 0;
                while (i < m) : (i += 1) {
                    ops.set( // b[i + j * ldb] = 0
                        &b[scast(u32, i + j * ldb)],
                        0,
                        ctx,
                    ) catch unreachable;
                }
            }

            return;
        }

        if (side == .left) {
            if (not) {
                if (uplo == .upper) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (noconj) {
                            if (ops.ne(alpha, 1, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.mul_( // b[i + j * ldb] *= alpha
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        alpha,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            var k: i32 = m - 1;
                            while (k >= 0) : (k -= 1) {
                                if (ops.ne(b[scast(u32, k + j * ldb)], 0, ctx) catch unreachable) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= a[k + k * lda]
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            a[scast(u32, k + k * lda)],
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = 0;
                                    while (i < k) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                a[scast(u32, i + k * lda)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                }
                            }
                        } else {
                            if (ops.ne(alpha, 1, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.mul_( // b[i + j * ldb] *= alpha
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        alpha,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            var k: i32 = m - 1;
                            while (k >= 0) : (k -= 1) {
                                if (ops.ne(b[scast(u32, k + j * ldb)], 0, ctx) catch unreachable) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= conj(a[k + k * lda])
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = 0;
                                    while (i < k) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * conj(a[i + k * lda])
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                ops.conj(a[scast(u32, i + k * lda)], ctx) catch unreachable,
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
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(u32, i + j * ldb)],
                                    alpha,
                                    b[scast(u32, i + j * ldb)],
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: i32 = 0;
                        while (k < m) : (k += 1) {
                            if (ops.ne(b[scast(u32, k + j * ldb)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= a[k + k * lda]
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            a[scast(u32, k + k * lda)],
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = k + 1;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                a[scast(u32, i + k * lda)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    if (diag == .non_unit) {
                                        ops.div_( // b[k + j * ldb] /= conj(a[k + k * lda])
                                            &b[scast(u32, k + j * ldb)],
                                            b[scast(u32, k + j * lda)],
                                            ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }

                                    var i: i32 = k + 1;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= b[k + j * ldb] * conj(a[i + k * lda])
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                b[scast(u32, k + j * ldb)],
                                                ops.conj(a[scast(u32, i + k * lda)], ctx) catch unreachable,
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
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
                        while (i < m) : (i += 1) {
                            var temp: T1 = scast(T1, ops.mul( // temp = alpha * b[i + j * ldb]
                                b[scast(u32, i + j * ldb)],
                                alpha,
                                ctx,
                            ) catch unreachable);

                            if (noconj) {
                                var k: i32 = 0;
                                while (k < i) : (k += 1) {
                                    ops.sub_( // temp -= a[k + i * lda] * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, k + i * lda)],
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= a[i + i * lda]
                                        &temp,
                                        temp,
                                        a[scast(u32, i + i * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: i32 = 0;
                                while (k < i) : (k += 1) {
                                    ops.sub_( // temp -= conj(a[k + i * lda]) * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= conj(a[i + i * lda])
                                        &temp,
                                        temp,
                                        ops.conj(a[scast(u32, i + i * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            ops.set( // b[i + j * ldb] = temp
                                &b[scast(u32, i + j * ldb)],
                                temp,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        var i: i32 = m - 1;
                        while (i >= 0) : (i -= 1) {
                            var temp: T1 = scast(T1, ops.mul( // temp = alpha * b[i + j * ldb]
                                alpha,
                                b[scast(u32, i + j * ldb)],
                                ctx,
                            ) catch unreachable);

                            if (noconj) {
                                var k: i32 = i + 1;
                                while (k < m) : (k += 1) {
                                    ops.sub_( // temp -= a[k + i * lda] * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, k + i * lda)],
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= a[i + i * lda]
                                        &temp,
                                        temp,
                                        a[scast(u32, i + i * lda)],
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else {
                                var k: i32 = i + 1;
                                while (k < m) : (k += 1) {
                                    ops.sub_( // temp -= conj(a[k + i * lda]) * b[k + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                            b[scast(u32, k + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (diag == .non_unit) {
                                    ops.div_( // temp /= conj(a[i + i * lda])
                                        &temp,
                                        temp,
                                        ops.conj(a[scast(u32, i + i * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }

                            ops.set( // b[i + j * ldb] = temp
                                &b[scast(u32, i + j * ldb)],
                                temp,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            }
        } else {
            if (not) {
                if (uplo == .upper) {
                    var j: i32 = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: i32 = 0;
                        while (k < j) : (k += 1) {
                            if (ops.ne(a[scast(u32, k + j * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                a[scast(u32, k + j * lda)],
                                                b[scast(u32, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= conj(a[k + j * lda]) * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                ops.conj(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                                b[scast(u32, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                }
                            }
                        }

                        if (diag == .non_unit) {
                            var temp: A = constants.one(A, ctx) catch unreachable;

                            if (noconj) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= temp
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                } else {
                    var j: i32 = n - 1;
                    while (j >= 0) : (j -= 1) {
                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= alpha
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var k: i32 = j + 1;
                        while (k < n) : (k += 1) {
                            if (ops.ne(a[scast(u32, k + j * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                a[scast(u32, k + j * lda)],
                                                b[scast(u32, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                } else {
                                    var i: i32 = 0;
                                    while (i < m) : (i += 1) {
                                        ops.sub_( // b[i + j * ldb] -= conj(a[k + j * lda]) * b[i + k * ldb]
                                            &b[scast(u32, i + j * ldb)],
                                            b[scast(u32, i + j * ldb)],
                                            ops.mul(
                                                ops.conj(a[scast(u32, k + j * lda)], ctx) catch unreachable,
                                                b[scast(u32, i + k * ldb)],
                                                ctx,
                                            ) catch unreachable,
                                            ctx,
                                        ) catch unreachable;
                                    }
                                }
                            }
                        }

                        if (diag == .non_unit) {
                            var temp: A = constants.one(A, ctx) catch unreachable;

                            if (noconj) {
                                ops.div_( // temp /= a[j + j * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, j + j * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp /= conj(a[j + j * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + j * ldb] *= temp
                                    &b[scast(u32, i + j * ldb)],
                                    b[scast(u32, i + j * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                if (uplo == .upper) {
                    var k: i32 = n - 1;
                    while (k >= 0) : (k -= 1) {
                        var temp: A = constants.one(A, ctx) catch unreachable;

                        if (diag == .non_unit) {
                            if (noconj) {
                                ops.div_( // temp = a[k + k * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, k + k * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp = conj(a[k + k * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= temp
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var j: i32 = 0;
                        while (j < k) : (j += 1) {
                            if (ops.ne(a[scast(u32, j + k * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    ops.set( // temp = a[j + k * lda]
                                        &temp,
                                        a[scast(u32, j + k * lda)],
                                        ctx,
                                    ) catch unreachable;
                                } else {
                                    ops.set( // temp = conj(a[j + k * lda])
                                        &temp,
                                        ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.sub_( // b[i + j * ldb] -= temp * b[i + k * ldb]
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        ops.mul(
                                            temp,
                                            b[scast(u32, i + k * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= alpha
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                } else {
                    var k: i32 = 0;
                    while (k < n) : (k += 1) {
                        var temp: A = constants.one(A, ctx) catch unreachable;

                        if (diag == .non_unit) {
                            if (noconj) {
                                ops.div_( // temp = a[k + k * lda]
                                    &temp,
                                    temp,
                                    a[scast(u32, k + k * lda)],
                                    ctx,
                                ) catch unreachable;
                            } else {
                                ops.div_( // temp = conj(a[k + k * lda])
                                    &temp,
                                    temp,
                                    ops.conj(a[scast(u32, k + k * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= temp
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
                                    temp,
                                    ctx,
                                ) catch unreachable;
                            }
                        }

                        var j: i32 = k + 1;
                        while (j < n) : (j += 1) {
                            if (ops.ne(a[scast(u32, j + k * lda)], 0, ctx) catch unreachable) {
                                if (noconj) {
                                    ops.set( // temp = a[j + k * lda]
                                        &temp,
                                        a[scast(u32, j + k * lda)],
                                        ctx,
                                    ) catch unreachable;
                                } else {
                                    ops.set( // temp = conj(a[j + k * lda])
                                        &temp,
                                        ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.sub_( // b[i + j * ldb] -= temp * b[i + k * ldb]
                                        &b[scast(u32, i + j * ldb)],
                                        b[scast(u32, i + j * ldb)],
                                        ops.mul(
                                            temp,
                                            b[scast(u32, i + k * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }

                        if (ops.ne(alpha, 1, ctx) catch unreachable) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                ops.mul_( // b[i + k * ldb] *= alpha
                                    &b[scast(u32, i + k * ldb)],
                                    b[scast(u32, i + k * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;
                            }
                        }
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.trsm not implemented for arbitrary precision types yet");
    }

    return;
}
