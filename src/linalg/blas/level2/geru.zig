const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const Order = types.Order;

/// Performs a rank-1 update of a general matrix.
///
/// The `geru` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Signature
/// ---------
/// ```zig
/// fn geru(order: Order, m: i32, n: i32, alpha: Al, x: [*]const X, incx: i32, y: [*]const Y, incy: i32, a: [*]A, lda: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`i32`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`i32`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`i32`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * k`, where `k` is `n` when `order` is `col_major`, or `m` when `order`
/// is `row_major`. On return, contains the result of the operation.
///
/// `lda` (`i32`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `m` or `n` are less than 0, if `lda`
/// is less than `max(1, m)` or `max(1, n)`, or if `incx` or `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub fn geru(
    order: Layout,
    m: i32,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.geru requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.geru requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.geru requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.geru requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.geru requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.geru requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.geru requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and Y == bool and A == bool)
        @compileError("zml.linalg.blas.geru does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.geru not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and options.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_cgeru(order.toCUInt(), scast(c_int, m), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_zgeru(order.toCUInt(), scast(c_int, m), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return _geru(order, m, n, alpha, x, incx, y, incy, a, lda, ctx);
}

fn _geru(
    order: Order,
    m: i32,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_geru(m, n, alpha, x, incx, y, incy, a, lda, ctx);
    } else {
        return k_geru(n, m, alpha, y, incy, x, incx, a, lda, ctx);
    }
}

fn k_geru(
    m: i32,
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C1: type = types.Coerce(Al, Y);
    const A: type = types.Child(@TypeOf(a));
    const CC: type = types.Coerce(Al, types.Coerce(X, types.Coerce(Y, A)));

    if (m < 0 or n < 0 or lda < int.max(1, m) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or ops.eq(alpha, 0, ctx) catch unreachable)
        return;

    var jy: i32 = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (incx == 1) {
            var j: i32 = 0;
            while (j < n) : (j += 1) {
                if (ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable) {
                    const temp: C1 = ops.mul( // temp = alpha * y[jy]
                        alpha,
                        y[scast(u32, jy)],
                        ctx,
                    ) catch unreachable;

                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // a[i + j * lda] += x[i] * temp
                            &a[scast(u32, i + j * lda)],
                            a[scast(u32, i + j * lda)],
                            ops.mul(
                                x[scast(u32, i)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                }

                jy += incy;
            }
        } else {
            const kx: i32 = if (incx < 0) (-m + 1) * incx else 0;

            var j: i32 = 0;
            while (j < n) : (j += 1) {
                if (ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable) {
                    const temp: C1 = ops.mul( // temp = alpha * y[jy]
                        alpha,
                        y[scast(u32, jy)],
                        ctx,
                    ) catch unreachable;

                    var ix: i32 = kx;
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        ops.add_( // a[i + j * lda] += x[ix] * temp
                            &a[scast(u32, i + j * lda)],
                            a[scast(u32, i + j * lda)],
                            ops.mul(
                                x[scast(u32, ix)],
                                temp,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ix += incx;
                    }
                }

                jy += incy;
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.geru not implemented for arbitrary precision types yet");
    }

    return;
}
