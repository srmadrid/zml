const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;
const Side = linalg.Side;
const Uplo = types.Uplo;

/// Computes a matrix-matrix product where one input matrix is Hermitian.
///
/// The `hemm` routines compute a scalar-matrix-matrix product using a Hermitian
/// matrix `A` and a general matrix `B` and add the result to a scalar-matrix
/// product using a general matrix `C`. The operation is defined as:
///
/// ```zig
///     C = alpha * A * B + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * B * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `A` is a Hermitian matrix, `B` and `C`
/// are `m`-by-`n` general matrices.
///
/// Signature
/// ---------
/// ```zig
/// fn hemm(order: Order, side: Side, uplo: Uplo, m: i32, n: i32, alpha: Al, a: [*]const A, lda: i32, b: [*]const B, ldb: i32, beta: Be, c: [*]C, ldc: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `side` (`Side`): Specifies whether the Hermitian matrix `A` appears on the
/// left or right in the operation as follows:
/// - If `side = left`, then the operation is `C = alpha * A * B + beta * C`.
/// - If `side = right`, then the operation is `C = alpha * B * A + beta * C`.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `m` (`i32`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`i32`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * ka`, where
/// `ka` is `m` if `side = left` and `n` if `side = right`.
///
/// `lda` (`i32`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `ldb * n` if
/// `order = col_major` and `ldb * m` if `order = row_major`.
///
/// `ldb` (`i32`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `beta`.
///
/// `c` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `ldc * n` if `order = col_major` and `ldc * m` if `order = row_major`.
///
/// `ldc` (`i32`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, m)` or `max(1, n)`, if `ldb` is less than `max(1, m)` or
/// `max(1, n)`, or if `ldc` is less than `max(1, m)` or `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub fn hemm(
    order: Layout,
    side: Side,
    uplo: Uplo,
    m: i32,
    n: i32,
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
        @compileError("zml.linalg.blas.hemm requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.hemm requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.hemm requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B))
        @compileError("zml.linalg.blas.hemm requires b to be a many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.hemm requires b's child type to numeric, got " ++ @typeName(B));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.hemm requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.hemm requires c to be a mutable many-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.hemm requires c's child type to be numeric, got " ++ @typeName(C));

    comptime if (Al == bool and A == bool and B == bool and Be == bool and C == bool)
        @compileError("zml.linalg.blas.hemm does not support alpha, a, b, beta and c all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(C))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.hemm not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and options.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_chemm(order.toCUInt(), side.toCUInt(), uplo.toCUInt(), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zhemm(order.toCUInt(), side.toCUInt(), uplo.toCUInt(), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return _hemm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

fn _hemm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: i32,
    n: i32,
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
        return k_hemm(
            side,
            uplo,
            m,
            n,
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
        return k_hemm(
            side.invert(),
            uplo.invert(),
            n,
            m,
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
    }
}

fn k_hemm(
    side: Side,
    uplo: Uplo,
    m: i32,
    n: i32,
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
    const T2: type = types.Coerce(A, B);
    const T3: type = types.Coerce(Al, A);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(B, types.Coerce(Be, C))));

    const nrowa: i32 = if (side == .left) m else n;

    if (m < 0 or n < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, m) or ldc < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            if (ops.eq(beta, 0, ctx) catch unreachable) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.set( // c[i + j * ldc] = 0
                            &c[scast(u32, i + j * ldc)],
                            0,
                            ctx,
                        ) catch unreachable;
                    }
                }
            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] *= beta
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            beta,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }

            return;
        }

        if (side == .left) {
            if (uplo == .upper) {
                var j: i32 = 0;
                while (j < n) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        const temp1: T1 = ops.mul( // temp1 = alpha * b[i + j * ldb];
                            alpha,
                            b[scast(u32, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                        var temp2: T2 = constants.zero(T2, ctx) catch unreachable;

                        var k: i32 = 0;
                        while (k < i) : (k += 1) {
                            ops.add_( // c[k + j * ldc] += temp1 * a[k + i * lda];
                                &c[scast(u32, k + j * ldc)],
                                c[scast(u32, k + j * ldc)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, k + i * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[k + j * ldb] * conj(a[k + i * lda]);
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(u32, k + j * ldb)],
                                    ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] = temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    alpha,
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                            ops.mul_( // c[i + j * ldc] *= beta;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                    var i: i32 = m - 1;
                    while (i >= 0) : (i -= 1) {
                        const temp1: T1 = ops.mul( // temp1 = alpha * b[i + j * ldb];
                            alpha,
                            b[scast(u32, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                        var temp2: T2 = constants.zero(T2, ctx) catch unreachable;

                        var k: i32 = i + 1;
                        while (k < m) : (k += 1) {
                            ops.add_( // c[k + j * ldc] += temp1 * a[k + i * lda];
                                &c[scast(u32, k + j * ldc)],
                                c[scast(u32, k + j * ldc)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, k + i * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += b[k + j * ldb] * conj(a[k + i * lda]);
                                &temp2,
                                temp2,
                                ops.mul(
                                    b[scast(u32, k + j * ldb)],
                                    ops.conj(a[scast(u32, k + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        if (ops.eq(beta, 0, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] = temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ops.mul(
                                    alpha,
                                    temp2,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
                            ops.mul_( // c[i + j * ldc] *= beta;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // c[i + j * ldc] += temp1 * re(a[i + i * lda]) + alpha * temp2;
                                &c[scast(u32, i + j * ldc)],
                                c[scast(u32, i + j * ldc)],
                                ops.add(
                                    ops.mul(
                                        temp1,
                                        ops.re(a[scast(u32, i + i * lda)], ctx) catch unreachable,
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
        } else {
            var j: i32 = 0;
            while (j < n) : (j += 1) {
                var temp1: T3 = ops.mul( // temp1 = alpha * re(a[j + j * lda]);
                    alpha,
                    ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                    ctx,
                ) catch unreachable;

                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] = temp1 * b[i + j * ldb];
                            &c[scast(u32, i + j * ldc)],
                            temp1,
                            b[scast(u32, i + j * ldb)],
                            ctx,
                        ) catch unreachable;
                    }
                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + j * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + j * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.mul_( // c[i + j * ldc] *= beta;
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            beta,
                            ctx,
                        ) catch unreachable;

                        ops.add_( // c[i + j * ldc] += temp1 * b[i + j * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + j * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }

                var k: i32 = 0;
                while (k < j) : (k += 1) {
                    if (uplo == .upper) {
                        ops.mul_( // temp1 = alpha * a[k + j * lda];
                            &temp1,
                            alpha,
                            a[scast(u32, k + j * lda)],
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // temp1 = alpha * conj(a[j + k * lda]);
                            &temp1,
                            alpha,
                            ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + k * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + k * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }

                k = j + 1;
                while (k < n) : (k += 1) {
                    if (uplo == .upper) {
                        ops.mul_( // temp1 = alpha * conj(a[j + k * lda]);
                            &temp1,
                            alpha,
                            ops.conj(a[scast(u32, j + k * lda)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.mul_( // temp1 = alpha * a[k + j * lda];
                            &temp1,
                            alpha,
                            a[scast(u32, k + j * lda)],
                            ctx,
                        ) catch unreachable;
                    }

                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // c[i + j * ldc] += temp1 * b[i + k * ldb];
                            &c[scast(u32, i + j * ldc)],
                            c[scast(u32, i + j * ldc)],
                            ops.mul(
                                temp1,
                                b[scast(u32, i + k * ldb)],
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.hemm not implemented for arbitrary precision types yet");
    }

    return;
}
