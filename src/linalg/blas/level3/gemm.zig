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

/// Computes a matrix-matrix product with general matrices.
///
/// The `gemm` routine computes a scalar-matrix-matrix product and adds the
/// result to a scalar-matrix product, with general matrices. The operation is
/// defined as:
///
/// ```zig
///     C = alpha * op(A) * op(B) + beta * C,
/// ```
///
/// where `op(X)` is `X`, `X^T`, `conj(X)`, or `X^H`, `alpha` and `beta` are
/// scalars, `op(A)` is an `m`-by-`k` matrix, `op(B)` is a `k`-by-`n`, and `C`
/// is an `m`-by-`n` matrix.
///
/// Signature
/// ---------
/// ```zig
/// fn gemm(order: Order, transa: Transpose, transb: Transpose, m: i32, n: i32, k: i32, alpha: Al, a: [*]const A, lda: i32, b: [*]const B, ldb: i32, beta: Be, c: [*]C, ldc: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the form of `op(A)`:
/// - If `transa = no_trans`, then `op(A) = A`.
/// - If `transa = trans`, then `op(A) = A^T`.
/// - If `transa = conj_no_trans`, then `op(A) = conj(A)`.
/// - If `transa = conj_trans`, then `op(A) = A^H`.
///
/// `transb` (`Transpose`): Specifies the form of `op(B)`:
/// - If `transb = no_trans`, then `op(B) = B`.
/// - If `transb = trans`, then `op(B) = B^T`.
/// - If `transb = conj_no_trans`, then `op(B) = conj(B)`.
/// - If `transb = conj_trans`, then `op(B) = B^H`.
///
/// `m` (`i32`): Specifies the number of rows of the matrix `op(A)` and of the
/// matrix `C`. Must be greater than or equal to 0.
///
/// `n` (`i32`): Specifies the number of columns of the matrix `op(B)` and the
/// number of columns of the matrix `C`. Must be greater than or equal to 0.
///
/// `k` (`i32`): Specifies the number of columns of the matrix `op(A)` and the
/// number of rows of the matrix `op(B)`. Must be greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `lda * k`                                       | `lda * m`                                 |
/// | `order = row_major` | `lda * m`                                       | `lda * k`                                 |
///
/// `lda` (`i32`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, m)`                                     | `max(1, k)`                               |
/// | `order = row_major` | `max(1, k)`                                     | `max(1, m)`                               |
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `ldb * n`                                       | `ldb * k`                                 |
/// | `order = row_major` | `ldb * k`                                       | `ldb * n`                                 |
///
/// `ldb` (`i32`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, k)`                                     | `max(1, n)`                               |
/// | `order = row_major` | `max(1, n)`                                     | `max(1, k)`                               |
///
/// `beta` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `beta`.
///
/// `c` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `ldc * n` if `order = col_major` or `ldc * m` if `order = row_major`.
///
/// `ldc` (`i32`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)` or `max(1, k)`, if `ldb` is less than `max(1, k)` or
/// `max(1, n)`, or if `ldc` is less than `max(1, m)` or `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub fn gemm(
    order: Layout,
    transa: Transpose,
    transb: Transpose,
    m: i32,
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
        @compileError("zml.linalg.blas.gemm requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.gemm requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.gemm requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B))
        @compileError("zml.linalg.blas.gemm requires b to be a many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.gemm requires b's child type to numeric, got " ++ @typeName(B));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.gemm requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.gemm requires c to be a mutable many-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.gemm requires c's child type to be numeric, got " ++ @typeName(C));

    comptime if (Al == bool and A == bool and B == bool and Be == bool and C == bool)
        @compileError("zml.linalg.blas.gemm does not support alpha, a, b, beta and c all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(C))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.gemm not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and options.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sgemm(order.toCUInt(), transa.toCUInt(), transb.toCUInt(), scast(c_int, m), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                } else if (comptime A == f64) {
                    return ci.cblas_dgemm(order.toCUInt(), transa.toCUInt(), transb.toCUInt(), scast(c_int, m), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_cgemm(order.toCUInt(), transa.toCUInt(), transb.toCUInt(), scast(c_int, m), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zgemm(order.toCUInt(), transa.toCUInt(), transb.toCUInt(), scast(c_int, m), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return _gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

fn _gemm(
    order: Order,
    transa: Transpose,
    transb: Transpose,
    m: i32,
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
        return k_gemm(
            transa,
            transb,
            m,
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
        return k_gemm(
            transb,
            transa,
            n,
            m,
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

fn k_gemm(
    transa: Transpose,
    transb: Transpose,
    m: i32,
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
    const T2: type = types.Coerce(A, B);
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(B, types.Coerce(Be, C))));

    const nota: bool = transa == .no_trans or transa == .conj_no_trans;
    const notb: bool = transb == .no_trans or transb == .conj_no_trans;
    const noconja: bool = transa == .no_trans or transa == .trans;
    const noconjb: bool = transb == .no_trans or transb == .trans;

    const nrowa: i32 = if (nota) m else k;
    const nrowb: i32 = if (notb) k else n;

    if (m < 0 or n < 0 or k < 0 or lda < int.max(1, nrowa) or ldb < int.max(1, nrowb) or ldc < int.max(1, m))
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or
        ((ops.eq(alpha, 0, ctx) catch unreachable or k == 0) and
            ops.eq(beta, 1, ctx) catch unreachable))
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
            } else {
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

        if (notb) {
            if (nota) {
                if (noconjb) {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[l + j * ldb]
                                    alpha,
                                    b[scast(u32, l + j * ldb)],
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
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
                    } else {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[l + j * ldb]
                                    b[scast(u32, l + j * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(u32, i + j * ldc)],
                                        c[scast(u32, i + j * ldc)],
                                        ops.mul(
                                            ops.conj(a[scast(u32, i + l * lda)], ctx) catch unreachable,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[l + j * ldb])
                                    ops.conj(b[scast(u32, l + j * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                        &c[scast(u32, i + j * ldc)],
                                        c[scast(u32, i + j * ldc)],
                                        ops.mul(
                                            a[scast(u32, i + l * lda)],
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
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[l + j * ldb])
                                    ops.conj(b[scast(u32, l + j * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(u32, i + j * ldc)],
                                        c[scast(u32, i + j * ldc)],
                                        ops.mul(
                                            ops.conj(a[scast(u32, i + l * lda)], ctx) catch unreachable,
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
                if (noconjb) {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * b[l + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, l + i * lda)],
                                            b[scast(u32, l + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * b[l + j * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, l + i * lda)], ctx) catch unreachable,
                                            b[scast(u32, l + j * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
                } else {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * conj(b[l + j * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, l + i * lda)],
                                            ops.conj(b[scast(u32, l + j * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * conj(b[l + j * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, l + i * lda)], ctx) catch unreachable,
                                            ops.conj(b[scast(u32, l + j * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
            }
        } else {
            if (nota) {
                if (noconjb) {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[j + l * ldb]
                                    b[scast(u32, j + l * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                        &c[scast(u32, i + j * ldc)],
                                        c[scast(u32, i + j * ldc)],
                                        ops.mul(
                                            a[scast(u32, i + l * lda)],
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
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * b[j + l * ldb]
                                    b[scast(u32, j + l * ldb)],
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(u32, i + j * ldc)],
                                        c[scast(u32, i + j * ldc)],
                                        ops.mul(
                                            ops.conj(a[scast(u32, i + l * lda)], ctx) catch unreachable,
                                            temp,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }
                            }
                        }
                    }
                } else {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[j + l * ldb])
                                    ops.conj(b[scast(u32, j + l * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * a[i + l * lda]
                                        &c[scast(u32, i + j * ldc)],
                                        c[scast(u32, i + j * ldc)],
                                        ops.mul(
                                            a[scast(u32, i + l * lda)],
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
                            if (ops.eq(beta, 0, ctx) catch unreachable) {
                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.set( // c[i + j * ldc] = 0
                                        &c[scast(u32, i + j * ldc)],
                                        0,
                                        ctx,
                                    ) catch unreachable;
                                }
                            } else if (ops.ne(beta, 1, ctx) catch unreachable) {
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

                            var l: i32 = 0;
                            while (l < k) : (l += 1) {
                                const temp: T1 = ops.mul( // temp = alpha * conj(b[j + l * ldb])
                                    ops.conj(b[scast(u32, j + l * ldb)], ctx) catch unreachable,
                                    alpha,
                                    ctx,
                                ) catch unreachable;

                                var i: i32 = 0;
                                while (i < m) : (i += 1) {
                                    ops.add_( // c[i + j * ldc] += temp * conj(a[i + l * lda])
                                        &c[scast(u32, i + j * ldc)],
                                        c[scast(u32, i + j * ldc)],
                                        ops.mul(
                                            ops.conj(a[scast(u32, i + l * lda)], ctx) catch unreachable,
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
                if (noconjb) {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * b[j + l * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, l + i * lda)],
                                            b[scast(u32, j + l * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * b[j + l * ldb]
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, l + i * lda)], ctx) catch unreachable,
                                            b[scast(u32, j + l * ldb)],
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
                } else {
                    if (noconja) {
                        var j: i32 = 0;
                        while (j < n) : (j += 1) {
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += a[l + i * lda] * conj(b[j + l * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            a[scast(u32, l + i * lda)],
                                            ops.conj(b[scast(u32, j + l * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
                            var i: i32 = 0;
                            while (i < m) : (i += 1) {
                                var temp: T2 = constants.zero(T2, ctx) catch unreachable;

                                var l: i32 = 0;
                                while (l < k) : (l += 1) {
                                    ops.add_( // temp += conj(a[l + i * lda]) * conj(b[j + l * ldb])
                                        &temp,
                                        temp,
                                        ops.mul(
                                            ops.conj(a[scast(u32, l + i * lda)], ctx) catch unreachable,
                                            ops.conj(b[scast(u32, j + l * ldb)], ctx) catch unreachable,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable;
                                }

                                if (ops.eq(beta, 0, ctx) catch unreachable) {
                                    ops.mul_( // c[i + j * ldc] = alpha * temp
                                        &c[scast(u32, i + j * ldc)],
                                        alpha,
                                        temp,
                                        ctx,
                                    ) catch unreachable;
                                } else if (ops.eq(beta, 1, ctx) catch unreachable) {
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
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.gemm not implemented for arbitrary precision types yet");
    }

    return;
}
