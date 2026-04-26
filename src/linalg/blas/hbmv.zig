const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Computes a matrix-vector product with a Hermitian band matrix defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian band matrix with `k` super-diagonals.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.hbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: Al, a: [*]const A, lda: isize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian band matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `beta`. When `beta` is
/// 0, then `y` need not be set on input.
///
/// `y` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` or `k` are less than 0, if `lda`
/// is less than `k + 1`, or if `incx` or `incy` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub fn hbmv(
    order: Layout,
    uplo: Uplo,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.hbmv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.hbmv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.hbmv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.hbmv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.hbmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.hbmv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.hbmv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.hbmv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.hbmv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.hbmv not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and options.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_chbmv(order.toCUInt(), uplo.toCUInt(), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zhbmv(order.toCUInt(), uplo.toCUInt(), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return _hbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, ctx);
}

fn _hbmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    if (order == .col_major) {
        return k_hbmv(
            uplo,
            n,
            k,
            alpha,
            a,
            lda,
            x,
            incx,
            beta,
            y,
            incy,
            true,
            ctx,
        );
    } else {
        return k_hbmv(
            uplo.invert(),
            n,
            k,
            ops.conj(alpha, ctx) catch unreachable,
            a,
            lda,
            x,
            incx,
            ops.conj(beta, ctx) catch unreachable,
            y,
            incy,
            false,
            ctx,
        );
    }
}

fn k_hbmv(
    uplo: Uplo,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    noconj: bool,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const A: type = types.Child(@TypeOf(a));
    const X: type = types.Child(@TypeOf(x));
    const C1: type = types.Coerce(Al, X);
    const C2: type = types.Coerce(A, X);
    const Be: type = @TypeOf(beta);
    const Y: type = types.Child(@TypeOf(y));
    const CC: type = types.Coerce(Al, types.Coerce(A, types.Coerce(X, types.Coerce(Be, Y))));

    if (n < 0 or k < 0 or lda < (k + 1) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or
        (ops.eq(alpha, 0, ctx) catch unreachable and ops.eq(beta, 1, ctx) catch unreachable))
        return;

    var kx: isize = if (incx < 0) (-n + 1) * incx else 0;
    var ky: isize = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        // First form  y = beta * y.
        if (ops.ne(beta, 1, ctx) catch unreachable) {
            if (incy == 1) {
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    for (0..scast(u32, n)) |i| {
                        ops.set( // y[i] = 0
                            &y[i],
                            0,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    if (noconj) {
                        for (0..scast(u32, n)) |i| {
                            ops.mul_( // y[i] *= beta
                                &y[i],
                                y[i],
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    } else {
                        for (0..scast(u32, n)) |i| {
                            ops.mul_( // y[i] = conj(y[i]) * beta
                                &y[i],
                                ops.conj(y[i], ctx) catch unreachable,
                                beta,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                }
            } else {
                var iy: isize = ky;
                if (ops.eq(beta, 0, ctx) catch unreachable) {
                    for (0..scast(u32, n)) |_| {
                        ops.set( // y[iy] = 0
                            &y[scast(u32, iy)],
                            0,
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                } else {
                    if (noconj) {
                        for (0..scast(u32, n)) |_| {
                            ops.mul_( // y[iy] *= beta
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                beta,
                                ctx,
                            ) catch unreachable;

                            iy += incy;
                        }
                    } else {
                        for (0..scast(u32, n)) |_| {
                            ops.mul_( // y[iy] = conj(y[iy]) * beta
                                &y[scast(u32, iy)],
                                ops.conj(y[scast(u32, iy)], ctx) catch unreachable,
                                beta,
                                ctx,
                            ) catch unreachable;

                            iy += incy;
                        }
                    }
                }
            }
        }

        if (ops.eq(alpha, 0, ctx) catch unreachable) {
            if (!noconj) {
                if (incy == 1) {
                    for (0..scast(u32, n)) |i| {
                        ops.conj_( // y[i] = conj(y[i])
                            &y[i],
                            y[i],
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
                    for (0..scast(u32, n)) |_| {
                        ops.conj_( // y[iy] = conj(y[iy])
                            &y[scast(u32, iy)],
                            y[scast(u32, iy)],
                            ctx,
                        ) catch unreachable;

                        iy += incy;
                    }
                }
            }

            return;
        }

        if (uplo == .upper) {
            // Form  y  when upper triangle of A is stored.
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                            alpha,
                            x[scast(u32, j)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        const l: isize = k - j;
                        var i: isize = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[i]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    x[scast(u32, i)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[j] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * x[jx]
                            alpha,
                            x[scast(u32, jx)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        var ix: isize = kx;
                        var iy: isize = ky;
                        const l: isize = k - j;
                        var i: isize = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[ix]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    x[scast(u32, ix)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                            iy += incy;
                        }

                        ops.add_( // y[jy] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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

                        jx += incx;
                        jy += incy;

                        if (j >= k) {
                            kx += incx;
                            ky += incy;
                        }
                    }
                }
            } else {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[j])
                            alpha,
                            ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        const l: isize = k - j;
                        var i: isize = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[i]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conj(x[scast(u32, i)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[j] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[jx])
                            alpha,
                            ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        var ix: isize = kx;
                        var iy: isize = ky;
                        const l: isize = k - j;
                        var i: isize = int.max(0, j - k);
                        while (i < j) : (i += 1) {
                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * conj(x[ix])
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ix += incx;
                            iy += incy;
                        }

                        ops.add_( // y[jy] += temp1 * re(a[k + j * lda]) + alpha * temp2
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.add(
                                ops.mul(
                                    temp1,
                                    ops.re(a[scast(u32, k + j * lda)], ctx) catch unreachable,
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

                        jx += incx;
                        jy += incy;

                        if (j >= k) {
                            kx += incx;
                            ky += incy;
                        }
                    }
                }
            }
        } else {
            // Form  y  when lower triangle of A is stored.
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * x[j]
                            alpha,
                            x[scast(u32, j)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[j] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: isize = -j;
                        var i: isize = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[i]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    x[scast(u32, i)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[j] += alpha * temp2
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.mul(
                                alpha,
                                temp2,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * x[jx]
                            alpha,
                            x[scast(u32, jx)],
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[jy] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: isize = -j;
                        var ix: isize = jx;
                        var iy: isize = jy;
                        var i: isize = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ix += incx;
                            iy += incy;

                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * x[ix]
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    x[scast(u32, ix)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[jy] += alpha * temp2
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.mul(
                                alpha,
                                temp2,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        jx += incx;
                        jy += incy;
                    }
                }
            } else {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[j])
                            alpha,
                            ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[j] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: isize = -j;
                        var i: isize = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ops.add_( // y[i] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, i)],
                                y[scast(u32, i)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * conj(x[i])
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conj(x[scast(u32, i)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[j] += alpha * temp2
                            &y[scast(u32, j)],
                            y[scast(u32, j)],
                            ops.mul(
                                alpha,
                                temp2,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        const temp1: C1 = ops.mul( // temp1 = alpha * conj(x[jx])
                            alpha,
                            ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                            ctx,
                        ) catch unreachable;
                        var temp2: C2 = constants.zero(C2, ctx) catch unreachable;

                        ops.add_( // y[jy] += temp1 * re(a[0 + j * lda])
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.mul(
                                temp1,
                                ops.re(a[scast(u32, 0 + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        const l: isize = -j;
                        var ix: isize = jx;
                        var iy: isize = jy;
                        var i: isize = j + 1;
                        while (i < int.min(n, j + k + 1)) : (i += 1) {
                            ix += incx;
                            iy += incy;

                            ops.add_( // y[iy] += temp1 * a[l + i + j * lda]
                                &y[scast(u32, iy)],
                                y[scast(u32, iy)],
                                ops.mul(
                                    temp1,
                                    a[scast(u32, l + i + j * lda)],
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            ops.add_( // temp2 += conj(a[l + i + j * lda]) * conj(x[ix])
                                &temp2,
                                temp2,
                                ops.mul(
                                    ops.conj(a[scast(u32, l + i + j * lda)], ctx) catch unreachable,
                                    ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
                                    ctx,
                                ) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        ops.add_( // y[jy] += alpha * temp2
                            &y[scast(u32, jy)],
                            y[scast(u32, jy)],
                            ops.mul(
                                alpha,
                                temp2,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        jx += incx;
                        jy += incy;
                    }
                }
            }
        }

        if (!noconj) {
            if (incy == 1) {
                for (0..scast(u32, n)) |i| {
                    ops.conj_( // y[i] = conj(y[i])
                        &y[i],
                        y[i],
                        ctx,
                    ) catch unreachable;
                }
            } else {
                var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
                for (0..scast(u32, n)) |_| {
                    ops.conj_( // y[iy] = conj(y[iy])
                        &y[scast(u32, iy)],
                        y[scast(u32, iy)],
                        ctx,
                    ) catch unreachable;

                    iy += incy;
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.hbmv not implemented for arbitrary precision types yet");
    }

    return;
}
