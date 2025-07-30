const std = @import("std");
const opts = @import("options");

const types = @import("../types.zig");
const scast = types.scast;
const validateContext = types.validateContext;
const Scalar = types.Scalar;
const Child = types.Child;
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const cfloat = @import("../cfloat.zig");
const cf32 = cfloat.cf32;
const cf64 = cfloat.cf64;
const ops = @import("../ops.zig");
const linalg = @import("../linalg.zig");

const ci = @import("../c.zig");

const Order = linalg.Order;
const Transpose = linalg.Transpose;
const Uplo = linalg.Uplo;
const Diag = linalg.Diag;
const Side = linalg.Side;

// Level 1 BLAS

/// Computes the sum of magnitudes of the vector elements.
///
/// The `asum_sub` routine computes the sum of the magnitudes of elements of a
/// real vector, or the sum of magnitudes of the real and imaginary parts of
/// elements of a complex vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n - 1].re) + abs(x[n - 1].im),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// `ret` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to where
/// the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is less than or equal
/// to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn asum_sub(
    n: isize,
    x: anytype,
    incx: isize,
    ret: anytype,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var R: type = @TypeOf(ret);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.asum_sub requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X) or X == bool)
        @compileError("zml.linalg.blas.asum_sub requires x's child type to be a non bool numeric, got " ++ @typeName(X));

    comptime if (!types.isPointer(R) or types.isConstPointer(R))
        @compileError("zml.linalg.blas.asum_sub requires ret to be a mutable one-item pointer, got " ++ @typeName(R));

    R = types.Child(R);

    comptime if (!types.isNumeric(R))
        @compileError("zml.linalg.blas.asum_sub requires ret's child type to be numeric, got " ++ @typeName(R));

    comptime if (types.isArbitraryPrecision(R)) {
        if (types.isArbitraryPrecision(X)) {
            @compileError("zml.linalg.blas.asum_sub not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.asum_sub not implemented for arbitrary precision types yet");
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            @compileError("zml.linalg.blas.asum_sub not implemented for arbitrary precision types yet");
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    if (comptime opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    try ops.set(ret, ci.cblas_sasum(scast(c_int, n), x, scast(c_int, incx)), ctx);
                } else if (comptime X == f64) {
                    try ops.set(ret, ci.cblas_dasum(scast(c_int, n), x, scast(c_int, incx)), ctx);
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    try ops.set(ret, ci.cblas_scasum(scast(c_int, n), x, scast(c_int, incx)), ctx);
                } else if (comptime Scalar(X) == f64) {
                    try ops.set(ret, ci.cblas_dzasum(scast(c_int, n), x, scast(c_int, incx)), ctx);
                }
            },
            else => {},
        }
    }

    return @import("blas/asum_sub.zig").asum_sub(X, n, x, incx, ret, ctx);
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `sasum_sub` routine computes the sum of the magnitudes of elements of a
/// vector:
///
/// ```zig
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n - 1]),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// `ret` (`*f32`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sasum_sub(
    n: isize,
    x: [*]const f32,
    incx: isize,
    ret: *f32,
) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `dasum_sub` routine computes the sum of the magnitudes of elements of a
/// vector:
///
/// ```zig
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n - 1]),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// `ret` (`*f64`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dasum_sub(
    n: isize,
    x: [*]const f64,
    incx: isize,
    ret: *f64,
) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `scasum_sub` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n - 1].re) + abs(x[n - 1].im),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// `ret` (`*f32`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn scasum_sub(
    n: isize,
    x: [*]const cf32,
    incx: isize,
    ret: *f32,
) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `dzasum_sub` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n - 1].re) + abs(x[n - 1].im),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// `ret` (`*f64`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dzasum_sub(
    n: isize,
    x: [*]const cf64,
    incx: isize,
    ret: *f64,
) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `asum` routine computes the sum of the magnitudes of elements of a real
/// vector, or the sum of magnitudes of the real and imaginary parts of elements
/// of a complex vector:
///
/// ```zig
///     abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n - 1].re) + abs(x[n - 1].im),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `Scalar(Child(@TypeOf(x)))`: The sum of magnitudes of real and imaginary
/// parts of all elements of the vector.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is less than or equal
/// to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn asum(
    n: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !Scalar(Child(@TypeOf(x))) {
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.asum requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X) or X == bool)
        @compileError("zml.linalg.blas.asum requires x's child type to be a non bool numeric, got " ++ @typeName(X));

    comptime if (types.isArbitraryPrecision(X)) {
        @compileError("zml.linalg.blas.asum not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_sasum(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime X == f64) {
                    return ci.cblas_dasum(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_scasum(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_dzasum(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/asum.zig").asum(n, x, incx, ctx);
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `sasum` routine computes the sum of the magnitudes of elements of a
/// vector:
///
/// ```zig
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n - 1]),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f32`: The sum of magnitudes of all elements of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sasum(
    n: isize,
    x: [*]const f32,
    incx: isize,
) f32 {
    return asum(n, x, incx, .{}) catch 0;
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `dasum` routine computes the sum of the magnitudes of elements of a
/// vector:
///
/// ```zig
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n - 1]),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f64`: The sum of magnitudes of all elements of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dasum(
    n: isize,
    x: [*]const f64,
    incx: isize,
) f64 {
    return asum(n, x, incx, .{}) catch 0;
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `scasum` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n - 1].re) + abs(x[n - 1].im),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f32`: The sum of magnitudes of the real and imaginary parts of all elements
/// of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn scasum(
    n: isize,
    x: [*]const cf32,
    incx: isize,
) f32 {
    return asum(n, x, incx, .{}) catch 0;
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `dzasum` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n - 1].re) + abs(x[n - 1].im),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f64`: The sum of magnitudes of the real and imaginary parts of all elements
/// of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dzasum(
    n: isize,
    x: [*]const cf64,
    incx: isize,
) f64 {
    return asum(n, x, incx, .{}) catch 0;
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `axpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y,
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with `n`
/// elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `1 + (n - 1) * abs(incy)`. On return contains the updated vector `y`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn axpy(
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.axpy requires alpha to be a numeric type, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.axpy requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);
    const C: type = Coerce(Al, X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.axpy requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.axpy requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.axpy requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and X == bool and Y == bool)
        @compileError("zml.linalg.blas.axpy does not support alpha, x and y all being bool");

    comptime if (types.isArbitraryPrecision(C)) {
        if (types.isArbitraryPrecision(Y)) {
            @compileError("zml.linalg.blas.axpy not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.axpy not implemented for arbitrary precision types yet");
        }
    } else {
        if (types.isArbitraryPrecision(Y)) {
            @compileError("zml.linalg.blas.axpy not implemented for arbitrary precision types yet");
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    if (comptime X == Y and types.canCoerce(Al, X) and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_saxpy(scast(c_int, n), scast(X, alpha), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (comptime X == f64) {
                    return ci.cblas_daxpy(scast(c_int, n), scast(X, alpha), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    const alpha_casted: X = scast(X, alpha);
                    return ci.cblas_caxpy(scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (comptime Scalar(X) == f64) {
                    const alpha_casted: X = scast(X, alpha);
                    return ci.cblas_zaxpy(scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/axpy.zig").axpy(n, alpha, x, incx, y, incy, ctx);
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `saxpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y,
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with `n`
/// elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn saxpy(
    n: isize,
    alpha: f32,
    x: [*]const f32,
    incx: isize,
    y: [*]f32,
    incy: isize,
) void {
    return axpy(n, alpha, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `daxpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y,
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with `n`
/// elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn daxpy(
    n: isize,
    alpha: f64,
    x: [*]const f64,
    incx: isize,
    y: [*]f64,
    incy: isize,
) void {
    return axpy(n, alpha, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `caxpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y,
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with `n`
/// elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn caxpy(
    n: isize,
    alpha: cf32,
    x: [*]const cf32,
    incx: isize,
    y: [*]cf32,
    incy: isize,
) void {
    return axpy(n, alpha, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `zaxpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y,
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with `n`
/// elements.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zaxpy(
    n: isize,
    alpha: cf64,
    x: [*]const cf64,
    incx: isize,
    y: [*]cf64,
    incy: isize,
) void {
    return axpy(n, alpha, x, incx, y, incy, .{}) catch {};
}

/// Copies a vector to another vector.
///
/// The `copy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = x,
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `1 + (n - 1) * abs(incy)`. On return contains the updated vector `y`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn copy(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.copy requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.copy requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.copy requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.copy requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (types.isArbitraryPrecision(X)) {
        if (types.isArbitraryPrecision(Y)) {
            @compileError("zml.linalg.blas.copy not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.copy not implemented for arbitrary precision types yet");
        }
    } else {
        if (types.isArbitraryPrecision(Y)) {
            @compileError("zml.linalg.blas.copy not implemented for arbitrary precision types yet");
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_scopy(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (comptime X == f64) {
                    return ci.cblas_dcopy(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_ccopy(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_zcopy(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/copy.zig").copy(n, x, incx, y, incy, ctx);
}

/// Copies a vector to another vector.
///
/// The `scopy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = x,
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn scopy(
    n: isize,
    x: [*]const f32,
    incx: isize,
    y: [*]f32,
    incy: isize,
) void {
    return copy(n, x, incx, y, incy, .{}) catch {};
}

/// Copies a vector to another vector.
///
/// The `dcopy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = x,
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dcopy(
    n: isize,
    x: [*]const f64,
    incx: isize,
    y: [*]f64,
    incy: isize,
) void {
    return copy(n, x, incx, y, incy, .{}) catch {};
}

/// Copies a vector to another vector.
///
/// The `ccopy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = x,
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ccopy(
    n: isize,
    x: [*]const cf32,
    incx: isize,
    y: [*]cf32,
    incy: isize,
) void {
    return copy(n, x, incx, y, incy, .{}) catch {};
}

/// Copies a vector to another vector.
///
/// The `zcopy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = x,
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zcopy(
    n: isize,
    x: [*]const cf64,
    incx: isize,
    y: [*]cf64,
    incy: isize,
) void {
    return copy(n, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `dot_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to where
/// the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn dot_sub(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ret: anytype,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var R: type = @TypeOf(ret);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.dot_sub requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.dot_sub requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.dot_sub requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);
    const C: type = Coerce(X, Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.dot_sub requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (X == bool and Y == bool)
        @compileError("zml.linalg.blas.dot_sub does not support x and y both being bool");

    comptime if (!types.isPointer(R) or types.isConstPointer(R))
        @compileError("zml.linalg.blas.dot_sub requires ret to be a mutable one-item pointer, got " ++ @typeName(R));

    R = types.Child(R);

    comptime if (!types.isNumeric(R))
        @compileError("zml.linalg.blas.dot_sub requires ret's child type to be numeric, got " ++ @typeName(R));

    comptime if (types.isArbitraryPrecision(R)) {
        if (types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dot_sub not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.dot_sub not implemented for arbitrary precision types yet");
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dot_sub not implemented for arbitrary precision types yet");
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    try ops.set(ret, ci.cblas_sdot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy)), ctx);
                } else if (comptime X == f64) {
                    try ops.set(ret, ci.cblas_ddot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy)), ctx);
                }
            },
            else => {},
        }
    }

    return @import("blas/dot_sub.zig").dot_sub(n, x, incx, y, incy, ret, ctx);
}

/// Computes a vector-vector dot product.
///
/// The `sdot_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*f32`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sdot_sub(
    n: isize,
    x: [*]const f32,
    incx: isize,
    y: [*]const f32,
    incy: isize,
    ret: *f32,
) void {
    return dot_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `ddot_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*f64`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ddot_sub(
    n: isize,
    x: [*]const f64,
    incx: isize,
    y: [*]const f64,
    incy: isize,
    ret: *f64,
) void {
    return dot_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `dot` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer to `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `Scalar(Coerce(Child(@TypeOf(x)), Child(@TypeOf(y))))`: The result of the
/// dot product.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn dot(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !Coerce(Child(@TypeOf(x)), Child(@TypeOf(y))) {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.dot requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.dot requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.dot requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);
    const C: type = Coerce(X, Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.dot requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (X == bool and Y == bool)
        @compileError("zml.linalg.blas.dot does not support x and y both being bool");

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.dot not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_sdot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (comptime X == f64) {
                    return ci.cblas_ddot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/dot.zig").dot(n, x, incx, y, incy, ctx);
}

/// Computes a vector-vector dot product.
///
/// The `sdot` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `f32`: The result of the dot product.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sdot(
    n: isize,
    x: [*]const f32,
    incx: isize,
    y: [*]const f32,
    incy: isize,
) f32 {
    return dot(n, x, incx, y, incy, .{}) catch 0;
}

/// Computes a vector-vector dot product.
///
/// The `ddot` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `f64`: The result of the dot product.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ddot(
    n: isize,
    x: [*]const f64,
    incx: isize,
    y: [*]const f64,
    incy: isize,
) f64 {
    return dot(n, x, incx, y, incy, .{}) catch 0;
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `dotc_sub` routine performs a vector-vector operation defined as:
///
/// ```zig
///     ret = conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n - 1]) * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer to `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to where
/// the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn dotc_sub(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ret: anytype,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var R: type = @TypeOf(ret);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.dotc_sub requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.dotc_sub requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.dotc_sub requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);
    const C: type = Coerce(X, Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.dotc_sub requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (X == bool and Y == bool)
        @compileError("zml.linalg.blas.dotc_sub does not support x and y both being bool");

    comptime if (!types.isPointer(R) or types.isConstPointer(R))
        @compileError("zml.linalg.blas.dotc_sub requires ret to be a mutable one-item pointer, got " ++ @typeName(R));

    R = types.Child(R);

    comptime if (!types.isNumeric(R))
        @compileError("zml.linalg.blas.dotc_sub requires ret's child type to be numeric, got " ++ @typeName(R));

    comptime if (types.isArbitraryPrecision(R)) {
        if (types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dotc_sub not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.dotc_sub not implemented for arbitrary precision types yet");
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dotc_sub not implemented for arbitrary precision types yet");
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_cdotc_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), ret);
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_zdotc_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), ret);
                }
            },
            else => {},
        }
    }

    return @import("blas/dotc_sub.zig").dotc_sub(n, x, incx, y, incy, ret, ctx);
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `cdotc_sub` routine performs a vector-vector operation defined as:
///
/// ```zig
///     ret = conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n - 1]) * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*cf32`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cdotc_sub(
    n: isize,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
    ret: *cf32,
) void {
    return dotc_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `zdotc_sub` routine performs a vector-vector operation defined as:
///
/// ```zig
///     ret = conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n - 1]) * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*cf64`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zdotc_sub(
    n: isize,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
    ret: *cf64,
) void {
    return dotc_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `dotc` routine performs a vector-vector operation defined as:
///
/// ```zig
///     conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n - 1]) * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer to `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `Coerce(Child(@TypeOf(x)), Child(@TypeOf(y)))`: The result of the dot
/// product.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn dotc(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !Coerce(Child(@TypeOf(x)), Child(@TypeOf(y))) {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.dotc requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.dotc requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.dotc requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);
    const C: type = Coerce(X, Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.dotc requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (X == bool and Y == bool)
        @compileError("zml.linalg.blas.dotc does not support x and y both being bool");

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.dotc not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    var temp: cf32 = undefined;
                    ci.cblas_cdotc_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), &temp);
                    return temp;
                } else if (comptime Scalar(X) == f64) {
                    var temp: cf64 = undefined;
                    ci.cblas_zdotc_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), &temp);
                    return temp;
                }
            },
            else => {},
        }
    }

    return @import("blas/dotc.zig").dotc(n, x, incx, y, incy, ctx);
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `cdotc` routine performs a vector-vector operation defined as:
///
/// ```zig
///     conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n - 1]) * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `cf32`: The result of the dot product.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cdotc(
    n: isize,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
) cf32 {
    return dotc(n, x, incx, y, incy, .{}) catch .{ .re = 0, .im = 0 };
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `zdotc` routine performs a vector-vector operation defined as:
///
/// ```zig
///     conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n - 1]) * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `cf64`: The result of the dot product.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zdotc(
    n: isize,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
) cf64 {
    return dotc(n, x, incx, y, incy, .{}) catch .{ .re = 0, .im = 0 };
}

/// Computes a vector-vector dot product.
///
/// The `dotu_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to where
/// the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn dotu_sub(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ret: anytype,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var R: type = @TypeOf(ret);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.dotu_sub requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.dotu_sub requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.dotu_sub requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);
    const C: type = Coerce(X, Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.dotu_sub requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (X == bool and Y == bool)
        @compileError("zml.linalg.blas.dotu_sub does not support x and y both being bool");

    comptime if (!types.isPointer(R) or types.isConstPointer(R))
        @compileError("zml.linalg.blas.dotu_sub requires ret to be a mutable one-item pointer, got " ++ @typeName(R));

    R = types.Child(R);

    comptime if (!types.isNumeric(R))
        @compileError("zml.linalg.blas.dotu_sub requires ret's child type to be numeric, got " ++ @typeName(R));

    comptime if (types.isArbitraryPrecision(R)) {
        if (types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dotu_sub not implemented for arbitrary precision types yet");
        } else {
            @compileError("zml.linalg.blas.dotu_sub not implemented for arbitrary precision types yet");
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            @compileError("zml.linalg.blas.dotu_sub not implemented for arbitrary precision types yet");
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_cdotu_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), ret);
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_zdotu_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), ret);
                }
            },
            else => {},
        }
    }

    return @import("blas/dotu_sub.zig").dotu_sub(n, x, incx, y, incy, ret, ctx);
}

/// Computes a vector-vector dot product.
///
/// The `cdotu_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*cf32`): Pointer to where the result will be stored.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cdotu_sub(
    n: isize,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
    ret: *cf32,
) void {
    return dotu_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `zdotu_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*cf64`): Pointer to the result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zdotu_sub(
    n: isize,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
    ret: *cf64,
) void {
    return dotu_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `dotu` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `Scalar(Coerce(Child(@TypeOf(x)), Child(@TypeOf(y))))`: The result of the
/// dot product.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn dotu(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !Coerce(Child(@TypeOf(x)), Child(@TypeOf(y))) {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.dotu requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.dotu requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.dotu requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);
    const C: type = Coerce(X, Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.dotu requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (X == bool and Y == bool)
        @compileError("zml.linalg.blas.dotu does not support x and y both being bool");

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.dotu not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    var temp: cf32 = undefined;
                    ci.cblas_cdotu_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), &temp);
                    return temp;
                } else if (comptime Scalar(X) == f64) {
                    var temp: cf64 = undefined;
                    ci.cblas_zdotu_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), &temp);
                    return temp;
                }
            },
            else => {},
        }
    }

    return @import("blas/dotu.zig").dotu(n, x, incx, y, incy, ctx);
}

/// Computes a vector-vector dot product.
///
/// The `cdotu` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `cf32`: The result of the dot product.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cdotu(
    n: isize,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
) cf32 {
    return dotu(n, x, incx, y, incy, .{}) catch .{ .re = 0, .im = 0 };
}

/// Computes a vector-vector dot product.
///
/// The `zdotu` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `cf64`: The result of the dot product.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zdotu(
    n: isize,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
) cf64 {
    return dotu(n, x, incx, y, incy, .{}) catch .{ .re = 0, .im = 0 };
}

/// Computes the Euclidean norm of a vector.
///
/// The `nrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// where `x` is a vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer to, `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `EnsureFloat(Scalar(Child(@TypeOf(x))))`: The Euclidean norm of the vector.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn nrm2(
    n: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !EnsureFloat(Scalar(Child(@TypeOf(x)))) {
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.nrm2 requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.nrm2 requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (X == bool)
        @compileError("zml.linalg.blas.nrm2 does not support x being bool");

    comptime if (types.isArbitraryPrecision(X)) {
        @compileError("zml.linalg.blas.nrm2 not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_snrm2(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime X == f64) {
                    return ci.cblas_dnrm2(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_scnrm2(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_dznrm2(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/nrm2.zig").nrm2(n, x, incx, ctx);
}

/// Computes the Euclidean norm of a vector.
///
/// The `snrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// where `x` is a vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f32`: The Euclidean norm of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn snrm2(
    n: isize,
    x: [*]const f32,
    incx: isize,
) f32 {
    return nrm2(n, x, incx, .{}) catch 0;
}

/// Computes the Euclidean norm of a vector.
///
/// The `dnrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// where `x` is a vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f64`: The Euclidean norm of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dnrm2(
    n: isize,
    x: [*]const f64,
    incx: isize,
) f64 {
    return nrm2(n, x, incx, .{}) catch 0;
}

/// Computes the Euclidean norm of a vector.
///
/// The `scnrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// where `x` is a vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f32`: The Euclidean norm of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn scnrm2(
    n: isize,
    x: [*]const cf32,
    incx: isize,
) f32 {
    return nrm2(n, x, incx, .{}) catch 0;
}

/// Computes the Euclidean norm of a vector.
///
/// The `dznrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// where `x` is a vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `f64`: The Euclidean norm of the vector.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dznrm2(
    n: isize,
    x: [*]const cf64,
    incx: isize,
) f64 {
    return nrm2(n, x, incx, .{}) catch 0;
}

/// Performs rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is r
/// eplaced as follows:
///
/// ```zig
///     x[i] = c * x[i] + s * y[i]
///     y[i] = c * y[i] - s * x[i]
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `c` (`bool`, `int`, `float`, `integer`, `rational`, `real` or `expression`):
/// The cosine of the rotation angle.
///
/// `s` (`bool`, `int`, `float`, `integer`, `rational`, `real` or `expression`):
/// The sine of the rotation angle.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn rot(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    c: anytype,
    s: anytype,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    const C: type = @TypeOf(c);
    const S: type = @TypeOf(s);

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.rot requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.rot requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.rot requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.rot requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.rot requires c to be numeric, got " ++ @typeName(C));

    comptime if (types.isComplex(C))
        @compileError("zml.linalg.blas.rot does not support c being complex, got " ++ @typeName(C));

    comptime if (!types.isNumeric(S))
        @compileError("zml.linalg.blas.rot requires s to be numeric, got " ++ @typeName(S));

    comptime if (types.isComplex(S))
        @compileError("zml.linalg.blas.rot does not support s being complex, got " ++ @typeName(S));

    comptime if (X == bool and Y == bool and C == bool and S == bool)
        @compileError("zml.linalg.blas.rot does not support x, y, c and s all being bool");

    comptime if (types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(C) or
        types.isArbitraryPrecision(S))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.rot not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and X == C and X == S and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_srot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), c, s);
                } else if (comptime X == f64) {
                    return ci.cblas_drot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), c, s);
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_csrot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), c, s);
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_zdrot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), c, s);
                }
            },
            else => {},
        }
    }

    return @import("blas/rot.zig").rot(n, x, incx, y, incy, c, s, ctx);
}

/// Performs rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is r
/// eplaced as follows:
///
/// ```zig
///     x[i] = c * x[i] + s * y[i]
///     y[i] = c * y[i] - s * x[i]
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `c` (`f32`): The cosine of the rotation angle.
///
/// `s` (`f32`) The sine of the rotation angle.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn srot(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, c: f32, s: f32) void {
    return rot(n, x, incx, y, incy, c, s, .{}) catch {};
}

/// Performs rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is r
/// eplaced as follows:
///
/// ```zig
///     x[i] = c * x[i] + s * y[i]
///     y[i] = c * y[i] - s * x[i]
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `c` (`f64`): The cosine of the rotation angle.
///
/// `s` (`f64`) The sine of the rotation angle.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn drot(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, c: f64, s: f64) void {
    return rot(n, x, incx, y, incy, c, s, .{}) catch {};
}

/// Performs rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is r
/// eplaced as follows:
///
/// ```zig
///     x[i] = c * x[i] + s * y[i]
///     y[i] = c * y[i] - s * x[i]
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `c` (`f32`): The cosine of the rotation angle.
///
/// `s` (`f32`) The sine of the rotation angle.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn csrot(n: isize, x: [*]cf32, incx: isize, y: [*]cf32, incy: isize, c: f32, s: f32) void {
    return rot(n, x, incx, y, incy, c, s, .{}) catch {};
}

/// Performs rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is r
/// eplaced as follows:
///
/// ```zig
///     x[i] = c * x[i] + s * y[i]
///     y[i] = c * y[i] - s * x[i]
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `c` (`f64`): The cosine of the rotation angle.
///
/// `s` (`f64`) The sine of the rotation angle.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zdrot(n: isize, x: [*]cf64, incx: isize, y: [*]cf64, incy: isize, c: f64, s: f64) void {
    return rot(n, x, incx, y, incy, c, s, .{}) catch {};
}

/// Computes the parameters for a Givens rotation.
///
/// Given the Cartesian coordinates `(a, b)` of a point, this routine returns
/// the parameters `c`, `s`, `r`, and `z` associated with the Givens rotation.
/// The parameters `c` and `s` define a unitary matrix such that:
///
/// ```zig
///     [ c s ] [ a ]   [ r ]
///     [-s c ] [ b ] = [ 0 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `a` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Provides the
/// `x`-coordinate of the point `p`. On return, it contains the parameter `r`
/// associated with the Givens rotation.
///
/// `b` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Provides the
/// `y`-coordinate of the point `p`. On return, it contains the parameter `z`
/// associated with the Givens rotation.
///
/// `c` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): On return, it
/// contains the parameter `c` associated with the Givens rotation.
///
/// `s` (mutable one-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): On return, it
/// contains the parameter `s` associated with the Givens rotation.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn rotg(
    a: anytype,
    b: anytype,
    c: anytype,
    s: anytype,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);
    comptime var C: type = @TypeOf(c);
    comptime var S: type = @TypeOf(s);

    comptime if (!types.isPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.rotg requires a to be a mutable one-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.rotg requires a's child type to be numeric, got " ++ @typeName(A));

    comptime if (!types.isPointer(B) or types.isConstPointer(B))
        @compileError("zml.linalg.blas.rotg requires b to be a mutable one-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.rotg requires b's child type to be numeric, got " ++ @typeName(B));

    comptime if (!types.isPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.rotg requires c to be a mutable one-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.rotg requires c to be numeric, got " ++ @typeName(C));

    comptime if (types.isComplex(C))
        @compileError("zml.linalg.blas.rotg does not support c being complex, got " ++ @typeName(C));

    comptime if (!types.isPointer(S) or types.isConstPointer(S))
        @compileError("zml.linalg.blas.rotg requires s to be a mutable one-item pointer, got " ++ @typeName(S));

    S = types.Child(S);

    comptime if (!types.isNumeric(S))
        @compileError("zml.linalg.blas.rotg requires s to be numeric, got " ++ @typeName(S));

    comptime if (A == bool and B == bool and C == bool and S == bool)
        @compileError("zml.linalg.blas.rotg does not support a, b, c and s all being bool");

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B) or
        types.isArbitraryPrecision(C) or
        types.isArbitraryPrecision(S))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.rotg not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and A == S and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_srotg(a, b, c, s);
                } else if (comptime A == f64) {
                    return ci.cblas_drotg(a, b, c, s);
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_crotg(a, b, c, s);
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_zrotg(a, b, c, s);
                }
            },
            else => {},
        }
    }

    return @import("blas/rotg.zig").rotg(a, b, c, s, ctx);
}

/// Computes the parameters for a Givens rotation.
///
/// Given the Cartesian coordinates `(a, b)` of a point, this routine returns
/// the parameters `c`, `s`, `r`, and `z` associated with the Givens rotation.
/// The parameters `c` and `s` define a unitary matrix such that:
///
/// ```zig
///     [ c s ] [ a ]   [ r ]
///     [-s c ] [ b ] = [ 0 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `a` (`*f32`): Provides the `x`-coordinate of the point `p`. On return, it
/// contains the parameter `r` associated with the Givens rotation.
///
/// `b` (`*f32`): Provides the `y`-coordinate of the point `p`. On return, it
/// contains the parameter `z` associated with the Givens rotation.
///
/// `c` (`*f32`): On return, it contains the parameter `c` associated with the
/// Givens rotation.
///
/// `s` (`*f32`): On return, it contains the parameter `s` associated with the
/// Givens rotation.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
    return rotg(a, b, c, s, .{}) catch {};
}

/// Computes the parameters for a Givens rotation.
///
/// Given the Cartesian coordinates `(a, b)` of a point, this routine returns
/// the parameters `c`, `s`, `r`, and `z` associated with the Givens rotation.
/// The parameters `c` and `s` define a unitary matrix such that:
///
/// ```zig
///     [ c s ] [ a ]   [ r ]
///     [-s c ] [ b ] = [ 0 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `a` (`*f64`): Provides the `x`-coordinate of the point `p`. On return, it
/// contains the parameter `r` associated with the Givens rotation.
///
/// `b` (`*f64`): Provides the `y`-coordinate of the point `p`. On return, it
/// contains the parameter `z` associated with the Givens rotation.
///
/// `c` (`*f64`): On return, it contains the parameter `c` associated with the
/// Givens rotation.
///
/// `s` (`*f64`): On return, it contains the parameter `s` associated with the
/// Givens rotation.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
    return rotg(a, b, c, s, .{}) catch {};
}

/// Computes the parameters for a Givens rotation.
///
/// Given the Cartesian coordinates `(a, b)` of a point, this routine returns
/// the parameters `c`, `s`, `r`, and `z` associated with the Givens rotation.
/// The parameters `c` and `s` define a unitary matrix such that:
///
/// ```zig
///     [ c s ] [ a ]   [ r ]
///     [-s c ] [ b ] = [ 0 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `a` (`*cf32`): Provides the `x`-coordinate of the point `p`. On return, it
/// contains the parameter `r` associated with the Givens rotation.
///
/// `b` (`*cf32`): Provides the `y`-coordinate of the point `p`. On return, it
/// contains the parameter `z` associated with the Givens rotation.
///
/// `c` (`*f32`): On return, it contains the parameter `c` associated with the
/// Givens rotation.
///
/// `s` (`*cf32`): On return, it contains the parameter `s` associated with the
/// Givens rotation.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn crotg(a: *cf32, b: *cf32, c: *f32, s: *cf32) void {
    return rotg(a, b, c, s, .{}) catch {};
}

/// Computes the parameters for a Givens rotation.
///
/// Given the Cartesian coordinates `(a, b)` of a point, this routine returns
/// the parameters `c`, `s`, `r`, and `z` associated with the Givens rotation.
/// The parameters `c` and `s` define a unitary matrix such that:
///
/// ```zig
///     [ c s ] [ a ]   [ r ]
///     [-s c ] [ b ] = [ 0 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `a` (`*cf64`): Provides the `x`-coordinate of the point `p`. On return, it
/// contains the parameter `r` associated with the Givens rotation.
///
/// `b` (`*cf64`): Provides the `y`-coordinate of the point `p`. On return, it
/// contains the parameter `z` associated with the Givens rotation.
///
/// `c` (`*f64`): On return, it contains the parameter `c` associated with the
/// Givens rotation.
///
/// `s` (`*cf64`): On return, it contains the parameter `s` associated with the
/// Givens rotation.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zrotg(a: *cf64, b: *cf64, c: *f64, s: *cf64) void {
    return rotg(a, b, c, s, .{}) catch {};
}

/// Performs modified Givens rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is
/// replaced as follows:
///
/// ```zig
///     [ x[i] ]     [ x[i] ]
///     [ y[i] ] = H [ y[i] ]
/// ```
///
/// for `i = 1` to `n`, where `H` is a modified Givens transformation matrix
/// whose values are stored in the `param[1]` through `param[4]` array. See
/// discussion on the `param` argument.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (mutable many-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`. On return, every element `x[i]` is replaced
/// by `h11 * x[i] + h12 * y[i]`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`. On return, every element `y[i]` is replaced
/// by `h21 * x[i] + h22 * y[i]`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `param` (many-item pointer to `bool`, `int`, `float`, `integer`, `rational`,
/// `real` or `expression`): Array, size 5. The elements of the `param` array are:
///
/// - param[0] contains a switch, flag.
/// - param[1-4] contain `h11`, `h21`, `h12`, and `h22`, respectively, the
/// components of the array `H`.
///
/// Depending on the values of flag, the components of `H` are set as follows:
///
/// - `flag = -1`:
///
/// ```zig
///         [ h11 h12 ]
///     H = [ h21 h22 ]
/// ```
///
/// - `flag = 0`:
///
/// ```zig
///          [   1 h12 ]
///     H =  [ h21   1 ]
/// ```
///
/// - `flag = 1`:
///
/// ```zig
///          [ h11   1 ]
///     H =  [  -1 h22 ]
/// ```
///
/// - `flag = 2`:
///
/// ```zig
///          [ 1 0 ]
///     H =  [ 0 1 ]
/// ```
///
/// In the last three cases, the matrix entries of 1, -1, and 0 are assumed
/// based on the value of flag and are not required to be set in the `param`
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn rotm(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    param: anytype,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var P: type = @TypeOf(param);

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.rotm requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.rotm requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (types.isComplex(X))
        @compileError("zml.linalg.blas.rotm does not support x being complex, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.rot requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.rotm requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (types.isComplex(Y))
        @compileError("zml.linalg.blas.rotm does not support y being complex, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(P))
        @compileError("zml.linalg.blas.rotm requires param to be a many-item pointer, got " ++ @typeName(P));

    P = types.Child(P);

    comptime if (!types.isNumeric(P))
        @compileError("zml.linalg.blas.rotm requires param's child type to be numeric, got " ++ @typeName(P));

    comptime if (types.isComplex(P))
        @compileError("zml.linalg.blas.rotm does not support param being complex, got " ++ @typeName(P));

    comptime if (X == bool and Y == bool and P == bool)
        @compileError("zml.linalg.blas.rotm does not support x, y and param all being bool");

    comptime if (types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(P))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.rotm not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and X == P and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_srotm(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), param);
                } else if (comptime X == f64) {
                    return ci.cblas_drotm(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), param);
                }
            },
            else => {},
        }
    }

    return @import("blas/rotm.zig").rotm(n, x, incx, y, incy, param, ctx);
}

/// Performs modified Givens rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is
/// replaced as follows:
///
/// ```zig
///     [ x[i] ]     [ x[i] ]
///     [ y[i] ] = H [ y[i] ]
/// ```
///
/// for `i = 1` to `n`, where `H` is a modified Givens transformation matrix
/// whose values are stored in the `param[1]` through `param[4]` array. See
/// discussion on the `param` argument.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`. On return,
/// every element `x[i]` is replaced by `h11 * x[i] + h12 * y[i]`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`. On return,
/// every element `y[i]` is replaced by `h21 * x[i] + h22 * y[i]`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `param` (`[*]const f32`): Array, size 5. The elements of the `param` array
/// are:
///
/// - param[0] contains a switch, flag.
/// - param[1-4] contain `h11`, `h21`, `h12`, and `h22`, respectively, the
/// components of the array `H`.
///
/// Depending on the values of flag, the components of `H` are set as follows:
///
/// - `flag = -1`:
///
/// ```zig
///         [ h11 h12 ]
///     H = [ h21 h22 ]
/// ```
///
/// - `flag = 0`:
///
/// ```zig
///          [   1 h12 ]
///     H =  [ h21   1 ]
/// ```
///
/// - `flag = 1`:
///
/// ```zig
///          [ h11   1 ]
///     H =  [  -1 h22 ]
/// ```
///
/// - `flag = 2`:
///
/// ```zig
///          [ 1 0 ]
///     H =  [ 0 1 ]
/// ```
///
/// In the last three cases, the matrix entries of 1, -1, and 0 are assumed
/// based on the value of flag and are not required to be set in the `param`
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn srotm(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, param: [*]const f32) void {
    return rotm(n, x, incx, y, incy, param, .{}) catch {};
}

/// Performs modified Givens rotation of points in the plane.
///
/// Given two vectors `x` and `y`, each vector element of these vectors is
/// replaced as follows:
///
/// ```zig
///     [ x[i] ]     [ x[i] ]
///     [ y[i] ] = H [ y[i] ]
/// ```
///
/// for `i = 1` to `n`, where `H` is a modified Givens transformation matrix
/// whose values are stored in the `param[1]` through `param[4]` array. See
/// discussion on the `param` argument.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`. On return,
/// every element `x[i]` is replaced by `h11 * x[i] + h12 * y[i]`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`. On return,
/// every element `y[i]` is replaced by `h21 * x[i] + h22 * y[i]`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `param` (`[*]const f64`): Array, size 5. The elements of the `param` array
/// are:
///
/// - param[0] contains a switch, flag.
/// - param[1-4] contain `h11`, `h21`, `h12`, and `h22`, respectively, the
/// components of the array `H`.
///
/// Depending on the values of flag, the components of `H` are set as follows:
///
/// - `flag = -1`:
///
/// ```zig
///         [ h11 h12 ]
///     H = [ h21 h22 ]
/// ```
///
/// - `flag = 0`:
///
/// ```zig
///          [   1 h12 ]
///     H =  [ h21   1 ]
/// ```
///
/// - `flag = 1`:
///
/// ```zig
///          [ h11   1 ]
///     H =  [  -1 h22 ]
/// ```
///
/// - `flag = 2`:
///
/// ```zig
///          [ 1 0 ]
///     H =  [ 0 1 ]
/// ```
///
/// In the last three cases, the matrix entries of 1, -1, and 0 are assumed
/// based on the value of flag and are not required to be set in the `param`
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn drotm(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, param: [*]const f64) void {
    return rotm(n, x, incx, y, incy, param, .{}) catch {};
}

/// Computes the parameters for a modified Givens rotation.
///
/// Given the Cartesian coordinates `(x1, y1)` of an input vector, this  routine
/// computes the components of a modified Givens transformation matrix `H` that
/// zeros the `y`-component of the resulting vector:
///
/// ```zig
///     [ x1 ]     [ x1 √d1 ]
///     [  0 ] = H [ y1 √d2 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `d1` (mutable one-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Provides the scaling factor for the
/// `x`-coordinate of the input vector. On return it provides the first diagonal
/// element of the updated matrix.
///
/// `d2` (mutable one-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Provides the scaling factor for the
/// `y`-coordinate of the input vector. On return it provides the second diagonal
/// element of the updated matrix.
///
/// `x1` (mutable one-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Provides the `x`-coordinate of the
/// input vector. On return it provides the `x`-coordinate of the rotated vector
/// before scaling.
///
/// `y1` (`bool`, `int`, `float`, `integer`, `rational`, `real` or
/// `expression`): Provides the `y`-coordinate of the input vector.
///
/// `param` (mutable many-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Array, size 5. On return the elements
/// of the `param` array are:
///
/// - param[0] contains a switch, flag.
/// - param[1-4] contain `h11`, `h21`, `h12`, and `h22`, respectively, the
/// components of the array `H`.
///
/// Depending on the values of flag, the components of `H` are set as follows:
///
/// - `flag = -1`:
///
/// ```zig
///         [ h11 h12 ]
///     H = [ h21 h22 ]
/// ```
///
/// - `flag = 0`:
///
/// ```zig
///          [   1 h12 ]
///     H =  [ h21   1 ]
/// ```
///
/// - `flag = 1`:
///
/// ```zig
///          [ h11   1 ]
///     H =  [  -1 h22 ]
/// ```
///
/// - `flag = 2`:
///
/// ```zig
///          [ 1 0 ]
///     H =  [ 0 1 ]
/// ```
///
/// In the last three cases, the matrix entries of 1, -1, and 0 are assumed
/// based on the value of flag and are not required to be set in the `param`
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn rotmg(
    d1: anytype,
    d2: anytype,
    x1: anytype,
    y1: anytype,
    param: anytype,
    ctx: anytype,
) !void {
    comptime var D1: type = @TypeOf(d1);
    comptime var D2: type = @TypeOf(d2);
    comptime var X1: type = @TypeOf(x1);
    const Y1: type = @TypeOf(y1);
    comptime var P: type = @TypeOf(param);

    comptime if (!types.isPointer(D1) or types.isConstPointer(D1))
        @compileError("zml.linalg.blas.rotmg requires d1 to be a mutable one-item pointer, got " ++ @typeName(D1));

    D1 = types.Child(D1);

    comptime if (!types.isNumeric(D1))
        @compileError("zml.linalg.blas.rotmg requires d1's child type to be numeric, got " ++ @typeName(D1));

    comptime if (types.isComplex(D1))
        @compileError("zml.linalg.blas.rotmg does not support d1 being complex, got " ++ @typeName(D1));

    comptime if (!types.isPointer(D2) or types.isConstPointer(D2))
        @compileError("zml.linalg.blas.rotmg requires d2 to be a mutable one-item pointer, got " ++ @typeName(D2));

    D2 = types.Child(D2);

    comptime if (!types.isNumeric(D2))
        @compileError("zml.linalg.blas.rotmg requires d2's child type to be numeric, got " ++ @typeName(D2));

    comptime if (types.isComplex(D2))
        @compileError("zml.linalg.blas.rotmg does not support d2 being complex, got " ++ @typeName(D2));

    comptime if (!types.isPointer(X1) or types.isConstPointer(X1))
        @compileError("zml.linalg.blas.rotmg requires x1 to be a mutable one-item pointer, got " ++ @typeName(X1));

    X1 = types.Child(X1);

    comptime if (!types.isNumeric(X1))
        @compileError("zml.linalg.blas.rotmg requires x1's child type to be numeric, got " ++ @typeName(X1));

    comptime if (types.isComplex(X1))
        @compileError("zml.linalg.blas.rotmg does not support x1 being complex, got " ++ @typeName(X1));

    comptime if (!types.isNumeric(Y1))
        @compileError("zml.linalg.blas.rotmg requires y1 to be numeric, got " ++ @typeName(Y1));

    comptime if (!types.isManyPointer(P) or types.isConstPointer(P))
        @compileError("zml.linalg.blas.rotmg requires param to be a mutable many-item pointer, got " ++ @typeName(P));

    P = types.Child(P);

    comptime if (!types.isNumeric(P))
        @compileError("zml.linalg.blas.rotmg requires param's child type to be numeric, got " ++ @typeName(P));

    comptime if (types.isComplex(P))
        @compileError("zml.linalg.blas.rotmg does not support param being complex, got " ++ @typeName(P));

    comptime if (D1 == bool and D2 == bool and X1 == bool and Y1 == bool and P == bool)
        @compileError("zml.linalg.blas.rotmg does not support d1, d2, x1, y1 and param all being bool");

    comptime if (types.isArbitraryPrecision(D1) or
        types.isArbitraryPrecision(D2) or
        types.isArbitraryPrecision(X1) or
        types.isArbitraryPrecision(Y1) or
        types.isArbitraryPrecision(P))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.rotmg not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime D1 == D2 and D1 == X1 and D1 == Y1 and D1 == P and opts.link_cblas != null) {
        switch (comptime types.numericType(D1)) {
            .float => {
                if (comptime D1 == f32) {
                    return ci.cblas_srotmg(d1, d2, x1, y1, param);
                } else if (comptime D1 == f64) {
                    return ci.cblas_drotmg(d1, d2, x1, y1, param);
                }
            },
            else => {},
        }
    }

    return @import("blas/rotmg.zig").rotmg(d1, d2, x1, y1, param, ctx);
}

/// Computes the parameters for a modified Givens rotation.
///
/// Given the Cartesian coordinates `(x1, y1)` of an input vector, this  routine
/// computes the components of a modified Givens transformation matrix `H` that
/// zeros the `y`-component of the resulting vector:
///
/// ```zig
///     [ x1 ]     [ x1 √d1 ]
///     [  0 ] = H [ y1 √d2 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `d1` (`*f32`): Provides the scaling factor for the `x`-coordinate of the
/// input vector. On return it provides the first diagonal element of the
/// updated matrix.
///
/// `d2` (`*f32`): Provides the scaling factor for the `y`-coordinate of the
/// input vector. On return it provides the second diagonal element of the
/// updated matrix.
///
/// `x1` (`*f32`): Provides the `x`-coordinate of the input vector. On return it
/// provides the `x`-coordinate of the rotated vector before scaling.
///
/// `y1` (`f32`): Provides the `y`-coordinate of the input vector.
///
/// `param` (`[*]f32`): Array, size 5. On return the elements of the `param`
/// array are:
///
/// - param[0] contains a switch, flag.
/// - param[1-4] contain `h11`, `h21`, `h12`, and `h22`, respectively, the
/// components of the array `H`.
///
/// Depending on the values of flag, the components of `H` are set as follows:
///
/// - `flag = -1`:
///
/// ```zig
///         [ h11 h12 ]
///     H = [ h21 h22 ]
/// ```
///
/// - `flag = 0`:
///
/// ```zig
///          [   1 h12 ]
///     H =  [ h21   1 ]
/// ```
///
/// - `flag = 1`:
///
/// ```zig
///          [ h11   1 ]
///     H =  [  -1 h22 ]
/// ```
///
/// - `flag = 2`:
///
/// ```zig
///          [ 1 0 ]
///     H =  [ 0 1 ]
/// ```
///
/// In the last three cases, the matrix entries of 1, -1, and 0 are assumed
/// based on the value of flag and are not required to be set in the `param`
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*]f32) void {
    return rotmg(d1, d2, x1, y1, param, .{}) catch {};
}

/// Computes the parameters for a modified Givens rotation.
///
/// Given the Cartesian coordinates `(x1, y1)` of an input vector, this  routine
/// computes the components of a modified Givens transformation matrix `H` that
/// zeros the `y`-component of the resulting vector:
///
/// ```zig
///     [ x1 ]     [ x1 √d1 ]
///     [  0 ] = H [ y1 √d2 ]
/// ```
///
/// The parameter `z` is defined such that if `|a| > |b|`, `z` is `s`; otherwise
/// if `c` is not 0 `z` is `1/c`; otherwise `z` is `1`.
///
/// Parameters
/// ----------
/// `d1` (`*f64`): Provides the scaling factor for the `x`-coordinate of the
/// input vector. On return it provides the first diagonal element of the
/// updated matrix.
///
/// `d2` (`*f64`): Provides the scaling factor for the `y`-coordinate of the
/// input vector. On return it provides the second diagonal element of the
/// updated matrix.
///
/// `x1` (`*f64`): Provides the `x`-coordinate of the input vector. On return it
/// provides the `x`-coordinate of the rotated vector before scaling.
///
/// `y1` (`f64`): Provides the `y`-coordinate of the input vector.
///
/// `param` (`[*]f64`): Array, size 5. On return the elements of the `param`
/// array are:
///
/// - param[0] contains a switch, flag.
/// - param[1-4] contain `h11`, `h21`, `h12`, and `h22`, respectively, the
/// components of the array `H`.
///
/// Depending on the values of flag, the components of `H` are set as follows:
///
/// - `flag = -1`:
///
/// ```zig
///         [ h11 h12 ]
///     H = [ h21 h22 ]
/// ```
///
/// - `flag = 0`:
///
/// ```zig
///          [   1 h12 ]
///     H =  [ h21   1 ]
/// ```
///
/// - `flag = 1`:
///
/// ```zig
///          [ h11   1 ]
///     H =  [  -1 h22 ]
/// ```
///
/// - `flag = 2`:
///
/// ```zig
///          [ 1 0 ]
///     H =  [ 0 1 ]
/// ```
///
/// In the last three cases, the matrix entries of 1, -1, and 0 are assumed
/// based on the value of flag and are not required to be set in the `param`
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`, `b`, `c` and `s`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*]f64) void {
    return rotmg(d1, d2, x1, y1, param, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x,
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `1 + (n - 1) * abs(incx)`. On return contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn scal(
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.scal requires alpha to be a numeric type, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.scal requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.scal requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (Al == bool and X == bool)
        @compileError("zml.linalg.blas.scal does not support alpha and x both being bool");

    comptime if (types.isArbitraryPrecision(Al) or types.isArbitraryPrecision(X)) {
        // When implemented, expand if
        @compileError("zml.linalg.blas.scal not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime types.canCoerce(Al, X) and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_sscal(scast(c_int, n), scast(X, alpha), x, scast(c_int, incx));
                } else if (comptime X == f64) {
                    return ci.cblas_dscal(scast(c_int, n), scast(X, alpha), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime types.isComplex(Al)) {
                    if (comptime Scalar(X) == f32) {
                        const alpha_casted: X = scast(X, alpha);
                        return ci.cblas_cscal(scast(c_int, n), &alpha_casted, x, scast(c_int, incx));
                    } else if (comptime Scalar(X) == f64) {
                        const alpha_casted: X = scast(X, alpha);
                        return ci.cblas_zscal(scast(c_int, n), &alpha_casted, x, scast(c_int, incx));
                    }
                } else {
                    if (comptime Scalar(X) == f32) {
                        return ci.cblas_csscal(scast(c_int, n), scast(Scalar(X), alpha), x, scast(c_int, incx));
                    } else if (comptime Scalar(X) == f64) {
                        return ci.cblas_zdscal(scast(c_int, n), scast(Scalar(X), alpha), x, scast(c_int, incx));
                    }
                }
            },
            else => {},
        }
    }

    return @import("blas/scal.zig").scal(n, alpha, x, incx, ctx);
}

/// Computes the product of a vector by a scalar.
///
/// The `sscal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x,
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sscal(
    n: isize,
    alpha: f32,
    x: [*]f32,
    incx: isize,
) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `dscal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x,
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dscal(
    n: isize,
    alpha: f64,
    x: [*]f64,
    incx: isize,
) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `cscal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x,
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cscal(
    n: isize,
    alpha: cf32,
    x: [*]cf32,
    incx: isize,
) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `zscal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x,
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zscal(
    n: isize,
    alpha: cf64,
    x: [*]cf64,
    incx: isize,
) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `csscal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x,
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn csscal(
    n: isize,
    alpha: f32,
    x: [*]cf32,
    incx: isize,
) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `zdscal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x,
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zdscal(n: isize, alpha: f64, x: [*]cf64, incx: isize) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `swap` routine returns vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `1 + (n - 1) * abs(incx)`. On return contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `1 + (n - 1) * abs(incy)`. On return contains the updated vector `y`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn swap(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.swap requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.swap requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.swap requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.swap requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.swap not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_sswap(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (comptime X == f64) {
                    return ci.cblas_dswap(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_cswap(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_zswap(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/swap.zig").swap(n, x, incx, y, incy, ctx);
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `sswap` routine returns vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`. On return
/// contains the updated vector `y`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sswap(
    n: isize,
    x: [*]f32,
    incx: isize,
    y: [*]f32,
    incy: isize,
) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `dswap` routine returns vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`. On return
/// contains the updated vector `y`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dswap(
    n: isize,
    x: [*]f64,
    incx: isize,
    y: [*]f64,
    incy: isize,
) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `cswap` routine returns vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incy)`. On return
/// contains the updated vector `y`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cswap(
    n: isize,
    x: [*]cf32,
    incx: isize,
    y: [*]cf32,
    incy: isize,
) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `zswap` routine returns vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`. On return
/// contains the updated vector `x`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incy)`. On return
/// contains the updated vector `y`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// Returns
/// -------
/// `void`: The result is stored in `x` and `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zswap(
    n: isize,
    x: [*]cf64,
    incx: isize,
    y: [*]cf64,
    incy: isize,
) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Finds the index of the element with maximum absolute value.
///
/// Given a vector `x`, the `iamax` routine returns the position of the vector
/// element `x[i]` that has the largest absolute value for real vectors, or the
/// largest sum `|x[i].re| + |x[i].im|` for complex vectors.
///
/// If more than one vector element is found with the same largest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the maximum absolute value in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is less than or equal
/// to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn iamax(
    n: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !usize {
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.iamax requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X) or X == bool)
        @compileError("zml.linalg.blas.iamax requires x's child type to be a non bool numeric, got " ++ @typeName(X));

    comptime if (types.isArbitraryPrecision(X)) {
        // When implemented, expand if
        // Might need but only when arbitrary p complex
        @compileError("zml.linalg.blas.iamax not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_isamax(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime X == f64) {
                    return ci.cblas_idamax(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_icamax(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_izamax(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/iamax.zig").iamax(n, x, incx, ctx);
}

/// Finds the index of the element with maximum absolute value.
///
/// Given a vector `x`, the `isamax` routine returns the position of the vector
/// element `x[i]` that has the largest absolute value.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same largest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the maximum absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn isamax(
    n: isize,
    x: [*]const f32,
    incx: isize,
) usize {
    return iamax(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with maximum absolute value.
///
/// Given a vector `x`, the `idamax` routine returns the position of the vector
/// element `x[i]` that has the largest absolute value.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same largest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the maximum absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn idamax(
    n: isize,
    x: [*]const f64,
    incx: isize,
) usize {
    return iamax(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with maximum absolute value.
///
/// Given a vector `x`, the `icamax` routine returns the position of the vector
/// element `x[i]` that has the largest sum `|x[i].re| + |x[i].im|`.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same largest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the maximum absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn icamax(
    n: isize,
    x: [*]const cf32,
    incx: isize,
) usize {
    return iamax(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with maximum absolute value.
///
/// Given a vector `x`, the `izamax` routine returns the position of the vector
/// element `x[i]` that has the largest sum `|x[i].re| + |x[i].im|`.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same largest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the maximum absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn izamax(
    n: isize,
    x: [*]const cf64,
    incx: isize,
) usize {
    return iamax(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `iamin` routine returns the position of the vector
/// element `x[i]` that has the smallest absolute value for real vectors, or the
/// smallest sum `|x[i].re| + |x[i].im|` for complex vectors.
///
/// If more than one vector element is found with the same smallest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the smallest absolute value in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is less than or equal
/// to 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn iamin(
    n: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !usize {
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.iamin requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X) or X == bool)
        @compileError("zml.linalg.blas.iamin requires x's child type to be a non bool numeric, got " ++ @typeName(X));

    comptime if (types.isArbitraryPrecision(X)) {
        // When implemented, expand if
        // Might need but only when arbitrary p complex
        @compileError("zml.linalg.blas.iamin not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (comptime X == f32) {
                    return ci.cblas_isamin(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime X == f64) {
                    return ci.cblas_idamin(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(X) == f32) {
                    return ci.cblas_icamin(scast(c_int, n), x, scast(c_int, incx));
                } else if (comptime Scalar(X) == f64) {
                    return ci.cblas_izamin(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/iamin.zig").iamin(n, x, incx, ctx);
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `isamin` routine returns the position of the vector
/// element `x[i]` that has the smallest absolute value.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same smallest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the smallest absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn isamin(
    n: isize,
    x: [*]const f32,
    incx: isize,
) usize {
    return iamin(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `idamin` routine returns the position of the vector
/// element `x[i]` that has the smallest absolute value.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same smallest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the smallest absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn idamin(
    n: isize,
    x: [*]const f64,
    incx: isize,
) usize {
    return iamin(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `icamin` routine returns the position of the vector
/// element `x[i]` that has the smallest sum `|x[i].re| + |x[i].im|`.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same smallest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the smallest absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn icamin(
    n: isize,
    x: [*]const cf32,
    incx: isize,
) usize {
    return iamin(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `izamin` routine returns the position of the vector
/// element `x[i]` that has the smallest sum `|x[i].re| + |x[i].im|`.
///
/// If either `n` or `incx` are not positive, the routine returns 0.
///
/// If more than one vector element is found with the same smallest absolute
/// value, the index of the first one encountered is returned.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// Returns
/// -------
/// `usize`: The index of the element with the smallest absolute value in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn izamin(
    n: isize,
    x: [*]const cf64,
    incx: isize,
) usize {
    return iamin(n, x, incx, .{}) catch 0;
}

// Level 2 BLAS

/// Computes a matrix-vector product with a general band matrix.
///
/// The `gbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` band matrix with `kl` sub-diagonals and `ku` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `kl` (`isize`): Specifies the number of sub-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `ku` (`isize`): Specifies the number of super-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `kl + ku + 1`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)` when `transa` is `no_transpose` or
/// `conj_no_transpose`, or `1 + (m - 1) * abs(incx)` otherwise.
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
/// `1 + (m - 1) * abs(incy)` when `transa` is `no_transpose` or
/// `conj_no_transpose`, or `1 + (n - 1) * abs(incy)` otherwise. On return,
/// contains the result of the operation.
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
/// `linalg.blas.Error.InvalidArgument`: If `m`, `n`, `kl` or `ku` are less than
/// 0, if `lda` is less than `kl + ku + 1`, or if `incx` or `incy` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn gbmv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
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
        @compileError("zml.linalg.blas.gbmv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.gbmv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.gbmv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.gbmv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.gbmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.gbmv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.gbmv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.gbmv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.gbmv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.gbmv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sgbmv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                } else if (comptime A == f64) {
                    return ci.cblas_dgbmv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_cgbmv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zgbmv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/gbmv.zig").gbmv(order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product with a general band matrix.
///
/// The `sgbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` band matrix with `kl` sub-diagonals and `ku` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `kl` (`isize`): Specifies the number of sub-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `ku` (`isize`): Specifies the number of super-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `kl + ku + 1`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f32`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]f32`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sgbmv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    x: [*]const f32,
    incx: isize,
    beta: f32,
    y: [*]f32,
    incy: isize,
) void {
    return gbmv(order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product with a general band matrix.
///
/// The `dgbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` band matrix with `kl` sub-diagonals and `ku` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `kl` (`isize`): Specifies the number of sub-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `ku` (`isize`): Specifies the number of super-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `kl + ku + 1`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f64`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]f64`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dgbmv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    x: [*]const f64,
    incx: isize,
    beta: f64,
    y: [*]f64,
    incy: isize,
) void {
    return gbmv(order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product with a general band matrix.
///
/// The `cgbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` band matrix with `kl` sub-diagonals and `ku` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `kl` (`isize`): Specifies the number of sub-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `ku` (`isize`): Specifies the number of super-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `kl + ku + 1`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf32`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cgbmv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    x: [*]const cf32,
    incx: isize,
    beta: cf32,
    y: [*]cf32,
    incy: isize,
) void {
    return gbmv(order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product with a general band matrix.
///
/// The `zgbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` band matrix, with `kl` sub-diagonals and `ku` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `kl` (`isize`): Specifies the number of sub-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `ku` (`isize`): Specifies the number of super-diagonals of the matrix `A`.
/// Must be greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `kl + ku + 1`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf64`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zgbmv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    kl: isize,
    ku: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    x: [*]const cf64,
    incx: isize,
    beta: cf64,
    y: [*]cf64,
    incy: isize,
) void {
    return gbmv(order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a general matrix.
///
/// The `gemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * k`, where
/// `k` is `n` when `order` is `col_major`, or `m` when `order` is `row_major`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)` when `transa` is `no_transpose` or
/// `conj_no_transpose`, or `1 + (m - 1) * abs(incx)` otherwise.
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
/// `1 + (m - 1) * abs(incy)` when `transa` is `no_transpose` or
/// `conj_no_transpose`, or `1 + (n - 1) * abs(incy)` otherwise. On return,
/// contains the result of the operation.
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
/// `linalg.blas.Error.InvalidArgument`: If `m` or `n` are less than 0, if `lda`
/// is less than `max(1, m)` or `max(1, n)`, or if `incx` or `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn gemv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
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
        @compileError("zml.linalg.blas.gemv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.gemv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.gemv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.gemv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.gemv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.gemv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.gemv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.gemv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.gemv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.gemv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sgemv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                } else if (comptime A == f64) {
                    return ci.cblas_dgemv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_cgemv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zgemv(@intFromEnum(order), @intFromEnum(transa), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/gemv.zig").gemv(order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product using a general matrix.
///
/// The `sgemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f32`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]f32`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sgemv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    x: [*]const f32,
    incx: isize,
    beta: f32,
    y: [*]f32,
    incy: isize,
) void {
    return gemv(order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a general matrix.
///
/// The `dgemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f64`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]f64`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dgemv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    x: [*]const f64,
    incx: isize,
    beta: f64,
    y: [*]f64,
    incy: isize,
) void {
    return gemv(order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a general matrix.
///
/// The `cgemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf32`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cgemv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    x: [*]const cf32,
    incx: isize,
    beta: cf32,
    y: [*]cf32,
    incy: isize,
) void {
    return gemv(order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a general matrix.
///
/// The `zgemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^T * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
///     y = alpha * A^H * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are vectors, `A` is an
/// `m`-by-`n` matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the operation to be performed on `A`:
/// - `no_transpose`: `y = alpha * A * x + beta * y`
/// - `transpose`: `y = alpha * A^T * x + beta * y`
/// - `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
/// - `conj_transpose`: `y = alpha * A^H * x + beta * y`
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (m - 1) * abs(incx)` otherwise.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf64`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (m - 1) * abs(incy)` when
/// `transa` is `no_transpose` or `conj_no_transpose`, or
/// `1 + (n - 1) * abs(incy)` otherwise. On return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zgemv(
    order: Order,
    transa: Transpose,
    m: isize,
    n: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    x: [*]const cf64,
    incx: isize,
    beta: cf64,
    y: [*]cf64,
    incy: isize,
) void {
    return gemv(order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Performs a rank-1 update of a general matrix.
///
/// The `ger` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * k`, where `k` is `n` when `order` is `col_major`, or `m` when `order`
/// is `row_major`. On return, contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
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
pub inline fn ger(
    order: Order,
    m: isize,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.ger requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.ger requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.ger requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.ger requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.ger requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.ger requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.ger requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and Y == bool and A == bool)
        @compileError("zml.linalg.blas.ger does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.ger not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sger(@intFromEnum(order), scast(c_int, m), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                } else if (comptime A == f64) {
                    return ci.cblas_dger(@intFromEnum(order), scast(c_int, m), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/ger.zig").ger(order, m, n, alpha, x, incx, y, incy, a, lda, ctx);
}

/// Performs a rank-1 update of a general matrix.
///
/// The `sger` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`. On return,
/// contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sger(
    order: Order,
    m: isize,
    n: isize,
    alpha: f32,
    x: [*]const f32,
    incx: isize,
    y: [*]const f32,
    incy: isize,
    a: [*]f32,
    lda: isize,
) void {
    return ger(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Performs a rank-1 update of a general matrix.
///
/// The `dger` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`. On return,
/// contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dger(
    order: Order,
    m: isize,
    n: isize,
    alpha: f64,
    x: [*]const f64,
    incx: isize,
    y: [*]const f64,
    incy: isize,
    a: [*]f64,
    lda: isize,
) void {
    return ger(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Performs a rank-1 update (conjugated) of a general matrix.
///
/// The `gerc` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(y^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * k`, where `k` is `n` when `order` is `col_major`, or `m` when `order`
/// is `row_major`. On return, contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
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
pub inline fn gerc(
    order: Order,
    m: isize,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.gerc requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.gerc requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.gerc requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.gerc requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.gerc requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.gerc requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.gerc requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and Y == bool and A == bool)
        @compileError("zml.linalg.blas.gerc does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.gerc not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_cgerc(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_zgerc(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/gerc.zig").gerc(order, m, n, alpha, x, incx, y, incy, a, lda, ctx);
}

/// Performs a rank-1 update (conjugated) of a general matrix.
///
/// The `cgerc` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(y^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`. On return,
/// contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cgerc(
    order: Order,
    m: isize,
    n: isize,
    alpha: cf32,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
    a: [*]cf32,
    lda: isize,
) void {
    return gerc(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Performs a rank-1 update (conjugated) of a general matrix.
///
/// The `zgerc` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(y^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`. On return,
/// contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zgerc(
    order: Order,
    m: isize,
    n: isize,
    alpha: cf64,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
    a: [*]cf64,
    lda: isize,
) void {
    return gerc(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

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
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * k`, where `k` is `n` when `order` is `col_major`, or `m` when `order`
/// is `row_major`. On return, contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
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
pub inline fn geru(
    order: Order,
    m: isize,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
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
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_cgeru(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_zgeru(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/geru.zig").geru(order, m, n, alpha, x, incx, y, incy, a, lda, ctx);
}

/// Performs a rank-1 update of a general matrix.
///
/// The `cgeru` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`. On return,
/// contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cgeru(
    order: Order,
    m: isize,
    n: isize,
    alpha: cf32,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
    a: [*]cf32,
    lda: isize,
) void {
    return geru(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Performs a rank-1 update of a general matrix.
///
/// The `zgeru` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `A`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (m - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * k`, where `k` is `n` when
/// `order` is `col_major`, or `m` when `order` is `row_major`. On return,
/// contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` when
/// `order` is `col_major`, or `max(1, n)` when `order` is `row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zgeru(
    order: Order,
    m: isize,
    n: isize,
    alpha: cf64,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
    a: [*]cf64,
    lda: isize,
) void {
    return geru(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Computes a matrix-vector product using a Hermitian band matrix.
///
/// The `hbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian band matrix with `k` super-diagonals.
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
pub inline fn hbmv(
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
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_chbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zhbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hbmv.zig").hbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product using a Hermitian band matrix.
///
/// The `chbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian band matrix, with `k` super-diagonals.
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
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf32`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn chbmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    k: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    x: [*]const cf32,
    incx: isize,
    beta: cf32,
    y: [*]cf32,
    incy: isize,
) void {
    return hbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a Hermitian band matrix.
///
/// The `zhbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian band matrix, with `k` super-diagonals.
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
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf64`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zhbmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    k: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    x: [*]const cf64,
    incx: isize,
    beta: cf64,
    y: [*]cf64,
    incy: isize,
) void {
    return hbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a Hermitian matrix.
///
/// The `hemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
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
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` or `incy` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn hemv(
    order: Order,
    uplo: Uplo,
    n: isize,
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
        @compileError("zml.linalg.blas.hemv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.hemv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.hemv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.hemv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.hemv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.hemv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.hemv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.hemv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.hemv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.hemv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_chemv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zhemv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hemv.zig").hemv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product using a Hermitian matrix.
///
/// The `chemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf32`): Specifies the scalar `beta`. When `beta` is
/// 0, then `y` need not be set on input.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn chemv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    x: [*]const cf32,
    incx: isize,
    beta: cf32,
    y: [*]cf32,
    incy: isize,
) void {
    return hemv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a Hermitian matrix.
///
/// The `zhemv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf64`): Specifies the scalar `beta`. When `beta` is
/// 0, then `y` need not be set on input.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zhemv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    x: [*]const cf64,
    incx: isize,
    beta: cf64,
    y: [*]cf64,
    incy: isize,
) void {
    return hemv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Performs a rank-1 update of a Hermitian matrix.
///
/// The `her` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(x^T) + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `integer`, `rational`, `real` or
/// `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * n`. On return, contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` or `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn her(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.her requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (types.isComplex(Al))
        @compileError("zml.linalg.blas.her does not support complex alpha, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.her requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.her requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.her requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.her requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and A == bool)
        @compileError("zml.linalg.blas.her does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.her not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and types.canCoerce(Al, Scalar(A)) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_cher(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(Scalar(A), alpha), x, scast(c_int, incx), a, scast(c_int, lda));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_zher(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(Scalar(A), alpha), x, scast(c_int, incx), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/her.zig").her(order, uplo, n, alpha, x, incx, a, lda, ctx);
}

/// Performs a rank-1 update of a Hermitian matrix.
///
/// The `cher` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(x^T) + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * n`. On return, contains the
/// result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cher(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    x: [*]const cf32,
    incx: isize,
    a: [*]cf32,
    lda: isize,
) void {
    return her(order, uplo, n, alpha, x, incx, a, lda, .{}) catch {};
}

/// Performs a rank-1 update of a Hermitian matrix.
///
/// The `zher` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(x^T) + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * n`. On return, contains the
/// result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zher(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    x: [*]const cf64,
    incx: isize,
    a: [*]cf64,
    lda: isize,
) void {
    return her(order, uplo, n, alpha, x, incx, a, lda, .{}) catch {};
}

/// Performs a rank-2 update of a Hermitian matrix.
///
/// The `her2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conjg(y^T) + conjg(alpha) * y * conjg(x^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * n`. On return, contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` or `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn her2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.her2 requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.her2 requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.her2 requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.her2 requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.her2 requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.her2 requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.her2 requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and Y == bool and A == bool)
        @compileError("zml.linalg.blas.her2 does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.her2 not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(Al)) {
            .cfloat => {
                if (comptime Scalar(Al) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_cher2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                } else if (comptime Scalar(Al) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_zher2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/her2.zig").her2(order, uplo, n, alpha, x, incx, y, incy, a, lda, ctx);
}

/// Performs a rank-2 update of a Hermitian matrix.
///
/// The `cher2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conjg(y^T) + conjg(alpha) * y * conjg(x^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * n`. On return, contains the
/// result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cher2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf32,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
    a: [*]cf32,
    lda: isize,
) void {
    return her2(order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Performs a rank-2 update of a Hermitian matrix.
///
/// The `zher2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conjg(y^T) + conjg(alpha) * y * conjg(x^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * n`. On return, contains the
/// result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zher2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf64,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
    a: [*]cf64,
    lda: isize,
) void {
    return her2(order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Computes a matrix-vector product using a Hermitian packed matrix.
///
/// The `hpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `ap` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`.
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
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` or
/// `incy` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn hpmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    ap: anytype,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.hpmv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.hpmv requires ap to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.hpmv requires ap's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.hpmv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.hpmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.hpmv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.hpmv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.hpmv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.hpmv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.hpmv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_chpmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, ap, x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zhpmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, ap, x, scast(c_int, incx), &beta_casted, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hpmv.zig").hpmv(order, uplo, n, alpha, ap, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product using a Hermitian packed matrix.
///
/// The `chpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `ap` (`[*]const cf32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf32`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn chpmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf32,
    ap: [*]const cf32,
    x: [*]const cf32,
    incx: isize,
    beta: cf32,
    y: [*]cf32,
    incy: isize,
) void {
    return hpmv(order, uplo, n, alpha, ap, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a Hermitian packed matrix.
///
/// The `zhpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `ap` (`[*]const cf64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`cf64`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zhpmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf64,
    ap: [*]const cf64,
    x: [*]const cf64,
    incx: isize,
    beta: cf64,
    y: [*]cf64,
    incy: isize,
) void {
    return hpmv(order, uplo, n, alpha, ap, x, incx, beta, y, incy, .{}) catch {};
}

/// Performs a rank-1 update of a Hermitian packed matrix.
///
/// The `hpr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(x^T) + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
///`uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `integer`, `rational`, `real` or
/// `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `ap` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` or
/// `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn hpr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    ap: anytype,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var A: type = @TypeOf(ap);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.hpr requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (types.isComplex(Al))
        @compileError("zml.linalg.blas.hpr does not support complex alpha, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.hpr requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.hpr requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.hpr requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.hpr requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and A == bool)
        @compileError("zml.linalg.blas.hpr does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.hpr not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_chpr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(Scalar(A), alpha), x, scast(c_int, incx), ap);
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_zhpr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(Scalar(A), alpha), x, scast(c_int, incx), ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/hpr.zig").hpr(order, uplo, n, alpha, x, incx, ap, ctx);
}

/// Performs a rank-1 update of a Hermitian packed matrix.
///
/// The `chpr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(x^T) + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
///`uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `ap` (`[*]cf32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn chpr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    x: [*]const cf32,
    incx: isize,
    ap: [*]cf32,
) void {
    return hpr(order, uplo, n, alpha, x, incx, ap, .{}) catch {};
}

/// Performs a rank-1 update of a Hermitian packed matrix.
///
/// The `zhpr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conj(x^T) + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
///`uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `ap` (`[*]cf64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zhpr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    x: [*]const cf64,
    incx: isize,
    ap: [*]cf64,
) void {
    return hpr(order, uplo, n, alpha, x, incx, ap, .{}) catch {};
}

/// Performs a rank-2 update of a Hermitian packed matrix.
///
/// The `hpr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conjg(y^T) + conjg(alpha) * y * conjg(x^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `ap` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`. On return, contains the result of the operation.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` or
/// `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn hpr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ap: anytype,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(ap);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.hpr2 requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.hpr2 requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.hpr2 requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.hpr2 requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.hpr2 requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.hpr2 requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.hpr2 requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and Y == bool and A == bool)
        @compileError("zml.linalg.blas.hpr2 does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.hpr2 not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_chpr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), ap);
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_zhpr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha_casted, x, scast(c_int, incx), y, scast(c_int, incy), ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/hpr2.zig").hpr2(order, uplo, n, alpha, x, incx, y, incy, ap, ctx);
}

/// Performs a rank-2 update of a Hermitian packed matrix.
///
/// The `chpr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conjg(y^T) + conjg(alpha) * y * conjg(x^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `ap` (`[*]cf32`): Array, size at least `(n * (n + 1)) / 2`. On return,
/// contains the result of the operation.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn chpr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf32,
    x: [*]const cf32,
    incx: isize,
    y: [*]const cf32,
    incy: isize,
    ap: [*]cf32,
) void {
    return hpr2(order, uplo, n, alpha, x, incx, y, incy, ap, .{}) catch {};
}

/// Performs a rank-2 update of a Hermitian packed matrix.
///
/// The `zhpr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * conjg(y^T) + conjg(alpha) * y * conjg(x^T) + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const cf64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `ap` (`[*]cf364`): Array, size at least `(n * (n + 1)) / 2`. On return,
/// contains the result of the operation.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zhpr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: cf64,
    x: [*]const cf64,
    incx: isize,
    y: [*]const cf64,
    incy: isize,
    ap: [*]cf64,
) void {
    return hpr2(order, uplo, n, alpha, x, incx, y, incy, ap, .{}) catch {};
}

/// Computes a matrix-vector product with a symmetric band matrix.
///
/// The `sbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetryc band matrix, with `k` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric band matrix `A` is used:
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
pub inline fn sbmv(
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
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.sbmv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.sbmv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.sbmv requires aT, 's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.sbmv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.sbmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.sbmv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.sbmv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.sbmv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.sbmv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.sbmv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                } else if (comptime A == f64) {
                    return ci.cblas_dsbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/sbmv.zig").sbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product with a symmetric band matrix.
///
/// The `ssbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetryc band matrix, with `k` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric band matrix `A` is used:
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
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f32`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`. On return,
/// contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ssbmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    k: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    x: [*]const f32,
    incx: isize,
    beta: f32,
    y: [*]f32,
    incy: isize,
) void {
    return sbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product with a symmetric band matrix.
///
/// The `dsbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetryc band matrix, with `k` super-diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric band matrix `A` is used:
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
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f64`): Specifies the scalar `beta`. When `beta` is 0, then `y` need
/// not be set on input.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`. On return,
/// contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dsbmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    k: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    x: [*]const f64,
    incx: isize,
    beta: f64,
    y: [*]f64,
    incy: isize,
) void {
    return sbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product with a symmetric packed matrix.
///
/// The `spmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `ap` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`.
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
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` or
/// `incy` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn spmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    ap: anytype,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.spmv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.spmv requires ap to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.spmv requires ap's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.spmv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.spmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.spmv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.spmv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.spmv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.spmv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.spmv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sspmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), ap, x, scast(c_int, incx), beta, y, scast(c_int, incy));
                } else if (comptime A == f64) {
                    return ci.cblas_dspmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), ap, x, scast(c_int, incx), beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/spmv.zig").spmv(order, uplo, n, alpha, ap, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product using a symmetric packed matrix.
///
/// The `sspmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `ap` (`[*]const f32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f32`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sspmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    ap: [*]const f32,
    x: [*]const f32,
    incx: isize,
    beta: f32,
    y: [*]f32,
    incy: isize,
) void {
    return spmv(order, uplo, n, alpha, ap, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a symmetric packed matrix.
///
/// The `dspmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `ap` (`[*]const f64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f64`): Specifies the scalar `beta`. When `beta` is 0, then `y`
/// need not be set on input.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the
/// operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dspmv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    ap: [*]const f64,
    x: [*]const f64,
    incx: isize,
    beta: f64,
    y: [*]f64,
    incy: isize,
) void {
    return spmv(order, uplo, n, alpha, ap, x, incx, beta, y, incy, .{}) catch {};
}

/// Performs a rank-1 update of a symmetric packed matrix.
///
/// The `spr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * x^T + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
///`uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `ap` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` or
/// `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn spr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    ap: anytype,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var A: type = @TypeOf(ap);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.spr requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.spr requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.spr requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.spr requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.spr requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and A == bool)
        @compileError("zml.linalg.blas.spr does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.spr not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sspr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), ap);
                } else if (comptime A == f64) {
                    return ci.cblas_dspr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/spr.zig").spr(order, uplo, n, alpha, x, incx, ap, ctx);
}

/// Performs a rank-1 update of a symmetric packed matrix.
///
/// The `sspr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * x^T + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
///`uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `ap` (`[*]f32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sspr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    x: [*]const f32,
    incx: isize,
    ap: [*]f32,
) void {
    return spr(order, uplo, n, alpha, x, incx, ap, .{}) catch {};
}

/// Performs a rank-1 update of a symmetric packed matrix.
///
/// The `dspr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * x^T + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
///`uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `ap` (`[*]f64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dspr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    x: [*]const f64,
    incx: isize,
    ap: [*]f64,
) void {
    return spr(order, uplo, n, alpha, x, incx, ap, .{}) catch {};
}

/// Performs a rank-2 update of a symmetric packed matrix.
///
/// The `spr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + alpha * y * x^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `ap` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`. On return, contains the result of the operation.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` or
/// `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn spr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ap: anytype,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(ap);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.spr2 requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.spr2 requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.spr2 requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.spr2 requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.spr2 requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.spr2 requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.spr2 requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and Y == bool and A == bool)
        @compileError("zml.linalg.blas.spr2 does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.spr2 not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sspr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), y, scast(c_int, incy), ap);
                } else if (comptime A == f64) {
                    return ci.cblas_dspr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), y, scast(c_int, incy), ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/spr2.zig").spr2(order, uplo, n, alpha, x, incx, y, incy, ap, ctx);
}

/// Performs a rank-2 update of a symmetric packed matrix.
///
/// The `sspr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + alpha * y * x^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `ap` (`[*]f32`): Array, size at least `(n * (n + 1)) / 2`. On return,
/// contains the result of the operation.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sspr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    x: [*]const f32,
    incx: isize,
    y: [*]const f32,
    incy: isize,
    ap: [*]f32,
) void {
    return spr2(order, uplo, n, alpha, x, incx, y, incy, ap, .{}) catch {};
}

/// Performs a rank-2 update of a symmetric packed matrix.
///
/// The `dspr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + alpha * y * x^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` symmetric matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// matrix `A` is supplied in the packed array `ap`:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// supplied in `ap`.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// supplied in `ap`.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `ap` (`[*]f64`): Array, size at least `(n * (n + 1)) / 2`. On return,
/// contains the result of the operation.
///
/// Returns
/// -------
/// `void`: The result is stored in `ap`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dspr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    x: [*]const f64,
    incx: isize,
    y: [*]const f64,
    incy: isize,
    ap: [*]f64,
) void {
    return spr2(order, uplo, n, alpha, x, incx, y, incy, ap, .{}) catch {};
}

/// Computes a matrix-vector product using a symmetric matrix.
///
/// The `symv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
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
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` or `incy` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn symv(
    order: Order,
    uplo: Uplo,
    n: isize,
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
        @compileError("zml.linalg.blas.symv requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.symv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.symv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.symv requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.symv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.symv requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.symv requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.symv requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (Al == bool and A == bool and X == bool and Be == bool and Y == bool)
        @compileError("zml.linalg.blas.symv does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(Y))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.symv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssymv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                } else if (comptime A == f64) {
                    return ci.cblas_dsymv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), x, scast(c_int, incx), scast(A, beta), y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/symv.zig").symv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, ctx);
}

/// Computes a matrix-vector product using a symmetric matrix.
///
/// The `ssymv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f32`): Specifies the scalar `beta`. When `beta` is
/// 0, then `y` need not be set on input.
///
/// `y` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ssymv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    x: [*]const f32,
    incx: isize,
    beta: f32,
    y: [*]f32,
    incy: isize,
) void {
    return symv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Computes a matrix-vector product using a symmetric matrix.
///
/// The `dsymv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are scalars, `x` and `y` are `n`-element vectors,
/// `A` is an `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `beta` (`f64`): Specifies the scalar `beta`. When `beta` is
/// 0, then `y` need not be set on input.
///
/// `y` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incy)`. On
/// return, contains the result of the operation.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `y`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dsymv(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    x: [*]const f64,
    incx: isize,
    beta: f64,
    y: [*]f64,
    incy: isize,
) void {
    return symv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch {};
}

/// Performs a rank-1 update of a symmetric matrix.
///
/// The `syr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * x^T + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `integer`, `rational`, `real` or
/// `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` or `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn syr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.syr requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.syr requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.syr requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.syr requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.syr requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and A == bool)
        @compileError("zml.linalg.blas.syr does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.syr not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssyr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), a, scast(c_int, lda));
                } else if (comptime A == f64) {
                    return ci.cblas_dsyr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr.zig").syr(order, uplo, n, alpha, x, incx, a, lda, ctx);
}

/// Performs a rank-1 update of a symmetric matrix.
///
/// The `ssyr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * x^T + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ssyr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    x: [*]const f32,
    incx: isize,
    a: [*]f32,
    lda: isize,
) void {
    return syr(order, uplo, n, alpha, x, incx, a, lda, .{}) catch {};
}

/// Performs a rank-1 update of a symmetric matrix.
///
/// The `dsyr` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * x^T + A,
/// ```
///
/// where `alpha` is a real scalar, `x` is an `n`-element vector, and `A` is an
/// `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dsyr(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    x: [*]const f64,
    incx: isize,
    a: [*]f64,
    lda: isize,
) void {
    return syr(order, uplo, n, alpha, x, incx, a, lda, .{}) catch {};
}

/// Performs a rank-2 update of a symmetric matrix.
///
/// The `syr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + alpha * y * x^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * n`. On return, contains the result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` or `incy` are 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn syr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.syr2 requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(X))
        @compileError("zml.linalg.blas.syr2 requires x to be a many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.syr2 requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (!types.isManyPointer(Y))
        @compileError("zml.linalg.blas.syr2 requires y to be a many-item pointer, got " ++ @typeName(Y));

    Y = types.Child(Y);

    comptime if (!types.isNumeric(Y))
        @compileError("zml.linalg.blas.syr2 requires y's child type to be numeric, got " ++ @typeName(Y));

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.blas.syr2 requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.syr2 requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (Al == bool and X == bool and Y == bool and A == bool)
        @compileError("zml.linalg.blas.syr2 does not support alpha, a, x, beta and y all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(A))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.syr2 not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and A == Y and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssyr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                } else if (comptime A == f64) {
                    return ci.cblas_dsyr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(A, alpha), x, scast(c_int, incx), y, scast(c_int, incy), a, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr2.zig").syr2(order, uplo, n, alpha, x, incx, y, incy, a, lda, ctx);
}

/// Performs a rank-2 update of a symmetric matrix.
///
/// The `ssyr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + alpha * y * x^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const f32`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * n`. On return, contains the
/// result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ssyr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f32,
    x: [*]const f32,
    incx: isize,
    y: [*]const f32,
    incy: isize,
    a: [*]f32,
    lda: isize,
) void {
    return syr2(order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Performs a rank-2 update of a symmetric matrix.
///
/// The `dsyr2` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     A = alpha * x * y^T + alpha * y * x^T + A,
/// ```
///
/// where `alpha` is a scalar, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` symmetric matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of the matrix `A` is
/// used.
/// - If `uplo = lower`, then the lower triangular part of the matrix `A` is
/// used.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// `y` (`[*]const f64`): Array, size at least `1 + (n - 1) * abs(incy)`.
///
/// `incy` (`isize`): Specifies the increment for indexing vector `y`. Must be
/// different from 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * n`. On return, contains the
/// result of the operation.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dsyr2(
    order: Order,
    uplo: Uplo,
    n: isize,
    alpha: f64,
    x: [*]const f64,
    incx: isize,
    y: [*]const f64,
    incy: isize,
    a: [*]f64,
    lda: isize,
) void {
    return syr2(order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular band matrix.
///
/// The `tbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular band matrix, with `k + 1` diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` or `k` are less than 0, if `lda`
/// is less than `k + 1`, or if `incx` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn tbmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.tbmv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.tbmv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.tbmv requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.tbmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (A == bool and X == bool)
        @compileError("zml.linalg.blas.tbmv does not support a and x both being bool");

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.tbmv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_stbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime A == f64) {
                    return ci.cblas_dtbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_ctbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_ztbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tbmv.zig").tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, ctx);
}

/// Computes a matrix-vector product using a triangular band matrix.
///
/// The `stbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular band matrix, with `k + 1` diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn stbmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const f32,
    lda: isize,
    x: [*]f32,
    incx: isize,
) void {
    return tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular band matrix.
///
/// The `dtbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular band matrix, with `k + 1` diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtbmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const f64,
    lda: isize,
    x: [*]f64,
    incx: isize,
) void {
    return tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular band matrix.
///
/// The `ctbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular band matrix, with `k + 1` diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctbmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const cf32,
    lda: isize,
    x: [*]cf32,
    incx: isize,
) void {
    return tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular band matrix.
///
/// The `ztbmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular band matrix, with `k + 1` diagonals.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztbmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const cf64,
    lda: isize,
    x: [*]cf64,
    incx: isize,
) void {
    return tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// band matrix.
///
/// The `tbsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular band matrix, with `k + 1` diagonals. The routine
/// does not test for singularity or near-singularity, such tests must be
/// performed before calling this routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`. Before entry, the incremented array `x` must
/// contain the n-element right-hand side vector `b`. On return, it contains the
/// solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` or `k` are less than 0, if `lda`
/// is less than `k + 1`, or if `incx` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn tbsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.tbsv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.tbsv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.tbsv requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.tbsv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (A == bool and X == bool)
        @compileError("zml.linalg.blas.tbsv does not support a and x both being bool");

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.tbsv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_stbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime A == f64) {
                    return ci.cblas_dtbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_ctbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_ztbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tbsv.zig").tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, ctx);
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// band matrix.
///
/// The `stbsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular band matrix, with `k + 1` diagonals. The routine
/// does not test for singularity or near-singularity, such tests must be
/// performed before calling this routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`. Before
/// entry, the incremented array `x` must contain the n-element right-hand side
/// vector `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn stbsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const f32,
    lda: isize,
    x: [*]f32,
    incx: isize,
) void {
    return tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// band matrix.
///
/// The `dtbsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular band matrix, with `k + 1` diagonals. The routine
/// does not test for singularity or near-singularity, such tests must be
/// performed before calling this routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`. Before
/// entry, the incremented array `x` must contain the n-element right-hand side
/// vector `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtbsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const f64,
    lda: isize,
    x: [*]f64,
    incx: isize,
) void {
    return tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// band matrix.
///
/// The `ctbsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular band matrix, with `k + 1` diagonals. The routine
/// does not test for singularity or near-singularity, such tests must be
/// performed before calling this routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`. Before
/// entry, the incremented array `x` must contain the n-element right-hand side vector `b`. On return, it contains the
/// solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctbsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const cf32,
    lda: isize,
    x: [*]cf32,
    incx: isize,
) void {
    return tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// band matrix.
///
/// The `ztbsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular band matrix, with `k + 1` diagonals. The routine
/// does not test for singularity or near-singularity, such tests must be
/// performed before calling this routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): Specifies the number of super-diagonals or sub-diagonals of
/// the matrix `A`. Must be greater than or equal to 0.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `k + 1`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`. Before
/// entry, the incremented array `x` must contain the n-element right-hand side vector `b`. On return, it contains the
/// solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztbsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]const cf64,
    lda: isize,
    x: [*]cf64,
    incx: isize,
) void {
    return tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular packed matrix.
///
/// The `tpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`.
///
/// `x` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` is
/// 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn tpmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: anytype,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.tpmv requires ap to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.tpmv requires ap's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.tpmv requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.tpmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (A == bool and X == bool)
        @compileError("zml.linalg.blas.tpmv does not support a and x both being bool");

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.tpmv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_stpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                } else if (comptime A == f64) {
                    return ci.cblas_dtpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_ctpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_ztpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tpmv.zig").tpmv(order, uplo, transa, diag, n, ap, x, incx, ctx);
}

/// Computes a matrix-vector product using a triangular packed matrix.
///
/// The `stpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const f32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn stpmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const f32,
    x: [*]f32,
    incx: isize,
) void {
    return tpmv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular packed matrix.
///
/// The `dtpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const f64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtpmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const f64,
    x: [*]f64,
    incx: isize,
) void {
    return tpmv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular packed matrix.
///
/// The `ctpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const cf32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctpmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const cf32,
    x: [*]cf32,
    incx: isize,
) void {
    return tpmv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular packed matrix.
///
/// The `ztpmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const cf64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztpmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const cf64,
    x: [*]cf64,
    incx: isize,
) void {
    return tpmv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// packed matrix.
///
/// The `tpsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// where `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(n * (n + 1)) / 2`.
///
/// `x` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`. On entry, the incremented array `x` must
/// contain the n-element right-hand side vector `b`. On return, it contains
/// the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, or if `incx` is
/// 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn tpsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: anytype,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.tpsv requires ap to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.tpsv requires ap's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.tpsv requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.tpsv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (A == bool and X == bool)
        @compileError("zml.linalg.blas.tpsv does not support a and x both being bool");

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.tpsv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_stpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                } else if (comptime A == f64) {
                    return ci.cblas_dtpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_ctpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_ztpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), ap, x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tpsv.zig").tpsv(order, uplo, transa, diag, n, ap, x, incx, ctx);
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// packed matrix.
///
/// The `stpsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// where `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const f32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn stpsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const f32,
    x: [*]f32,
    incx: isize,
) void {
    return tpsv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// packed matrix.
///
/// The `dtpsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// where `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const f64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtpsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const f64,
    x: [*]f64,
    incx: isize,
) void {
    return tpsv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// packed matrix.
///
/// The `ctpsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// where `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const cf32`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctpsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const cf32,
    x: [*]cf32,
    incx: isize,
) void {
    return tpsv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// packed matrix.
///
/// The `ztpsv` routine solves one of the following systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// where `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix, supplied in packed form.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `ap` (`[*]const cf64`): Array, size at least `(n * (n + 1)) / 2`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztpsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    ap: [*]const cf64,
    x: [*]cf64,
    incx: isize,
) void {
    return tpsv(order, uplo, transa, diag, n, ap, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular matrix.
///
/// The `trmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn trmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.trmv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.trmv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.trmv requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.trmv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (A == bool and X == bool)
        @compileError("zml.linalg.blas.trmv does not support a and x both being bool");

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.trmv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_strmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime A == f64) {
                    return ci.cblas_dtrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_ctrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_ztrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/trmv.zig").trmv(order, uplo, transa, diag, n, a, lda, x, incx, ctx);
}

/// Computes a matrix-vector product using a triangular matrix.
///
/// The `strmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn strmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const f32,
    lda: isize,
    x: [*]f32,
    incx: isize,
) void {
    return trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular matrix.
///
/// The `dtrmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtrmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const f64,
    lda: isize,
    x: [*]f64,
    incx: isize,
) void {
    return trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular matrix.
///
/// The `ctrmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctrmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const cf32,
    lda: isize,
    x: [*]cf32,
    incx: isize,
) void {
    return trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

/// Computes a matrix-vector product using a triangular matrix.
///
/// The `ztrmv` routine performs a matrix-vector operation defined as:
///
/// ```zig
///     x = A * x,
/// ```
///
/// or
///
/// ```zig
///     x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^T * x,
/// ```
///
/// or
///
/// ```zig
///     x = A^H * x,
/// ```
///
/// `x` is an `n`-element vector, `A` is an `n`-by-`n` unit, or non-unit, upper
/// or lower triangular matrix.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the operation to be performed:
/// - If `transa = no_trans`, then the operation is `x = A * x`.
/// - If `transa = trans`, then the operation is `x = A^T * x`.
/// - If `transa = conj_no_trans`, then the operation is `x = conj(A) * x`.
/// - If `transa = conj_trans`, then the operation is `x = A^H * x`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztrmv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const cf64,
    lda: isize,
    x: [*]cf64,
    incx: isize,
) void {
    return trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// matrix.
///
/// The `trsv` routine solves one of the systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular matrix. The routine does not test for singularity
/// or near-singularity, such tests must be performed before calling this
/// routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`. On entry, the incremented array `x` must
/// contain the n-element right-hand side vector `b`. On return, it contains
/// the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)`, or if `incx` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn trsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: anytype,
    lda: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.trsv requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.trsv requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.trsv requires x to be a mutable many-item pointer, got " ++ @typeName(X));

    X = types.Child(X);

    comptime if (!types.isNumeric(X))
        @compileError("zml.linalg.blas.trsv requires x's child type to be numeric, got " ++ @typeName(X));

    comptime if (A == bool and X == bool)
        @compileError("zml.linalg.blas.trsv does not support a and x both being bool");

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(X))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.trsv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == X and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_strsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime A == f64) {
                    return ci.cblas_dtrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_ctrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_ztrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, n), a, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/trsv.zig").trsv(order, uplo, transa, diag, n, a, lda, x, incx, ctx);
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// matrix.
///
/// The `strsv` routine solves one of the systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular matrix. The routine does not test for singularity
/// or near-singularity, such tests must be performed before calling this
/// routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]f32`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn strsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const f32,
    lda: isize,
    x: [*]f32,
    incx: isize,
) void {
    return trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// matrix.
///
/// The `dtrsv` routine solves one of the systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular matrix. The routine does not test for singularity
/// or near-singularity, such tests must be performed before calling this
/// routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]f64`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtrsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const f64,
    lda: isize,
    x: [*]f64,
    incx: isize,
) void {
    return trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// matrix.
///
/// The `ctrsv` routine solves one of the systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular matrix. The routine does not test for singularity
/// or near-singularity, such tests must be performed before calling this
/// routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]cf32`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctrsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const cf32,
    lda: isize,
    x: [*]cf32,
    incx: isize,
) void {
    return trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

/// Solves a system of linear equations whose coefficients are in a triangular
/// matrix.
///
/// The `ztrsv` routine solves one of the systems of equations:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or non-unit,
/// upper or lower triangular matrix. The routine does not test for singularity
/// or near-singularity, such tests must be performed before calling this
/// routine.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
/// triangular matrix:
/// - If `uplo = upper`, then the matrix is upper triangular.
/// - If `uplo = lower`, then the matrix is lower triangular.
///
/// `transa` (`Transpose`): Specifies the system of equations to be solved:
/// - If `transa = no_trans`, then the system is `A * x = b`.
/// - If `transa = trans`, then the system is `A^T * x = b`.
/// - If `transa = conj_no_trans`, then the system is `conj(A) * x = b`.
/// - If `transa = conj_trans`, then the system is `A^H * x = b`.
///
/// `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular:
/// - If `diag = unit`, then the matrix is unit triangular.
/// - If `diag = non_unit`, then the matrix is non-unit triangular.
///
/// `n` (`isize`): Specifies the order of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// `x` (`[*]cf64`): Array, size at least `1 + (n - 1) * abs(incx)`. On entry,
/// the incremented array `x` must contain the n-element right-hand side vector
/// `b`. On return, it contains the solution vector `x`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `x`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztrsv(
    order: Order,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    n: isize,
    a: [*]const cf64,
    lda: isize,
    x: [*]cf64,
    incx: isize,
) void {
    return trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch {};
}

// Level 3 BLAS

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
/// `m` (`isize`): Specifies the number of rows of the matrix `op(A)` and of the
/// matrix `C`. Must be greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `op(B)` and the
/// number of columns of the matrix `C`. Must be greater than or equal to 0.
///
/// `k` (`isize`): Specifies the number of columns of the matrix `op(A)` and the
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
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
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
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
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
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
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
pub inline fn gemm(
    order: Order,
    transa: Transpose,
    transb: Transpose,
    m: isize,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
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
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_sgemm(@intFromEnum(order), @intFromEnum(transa), @intFromEnum(transb), scast(c_int, m), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                } else if (comptime A == f64) {
                    return ci.cblas_dgemm(@intFromEnum(order), @intFromEnum(transa), @intFromEnum(transb), scast(c_int, m), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_cgemm(@intFromEnum(order), @intFromEnum(transa), @intFromEnum(transb), scast(c_int, m), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zgemm(@intFromEnum(order), @intFromEnum(transa), @intFromEnum(transb), scast(c_int, m), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/gemm.zig").gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

/// Computes a matrix-matrix product with general matrices.
///
/// The `sgemm` routine computes a scalar-matrix-matrix product and adds the
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
/// `m` (`isize`): Specifies the number of rows of the matrix `op(A)` and of the
/// matrix `C`. Must be greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `op(B)` and the
/// number of columns of the matrix `C`. Must be greater than or equal to 0.
///
/// `k` (`isize`): Specifies the number of columns of the matrix `op(A)` and the
/// number of rows of the matrix `op(B)`. Must be greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `lda * k`                                       | `lda * m`                                 |
/// | `order = row_major` | `lda * m`                                       | `lda * k`                                 |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, m)`                                     | `max(1, k)`                               |
/// | `order = row_major` | `max(1, k)`                                     | `max(1, m)`                               |
///
/// `b` (`[*]const f32`): Array, size at least:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `ldb * n`                                       | `ldb * k`                                 |
/// | `order = row_major` | `ldb * k`                                       | `ldb * n`                                 |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, k)`                                     | `max(1, n)`                               |
/// | `order = row_major` | `max(1, n)`                                     | `max(1, k)`                               |
///
/// `beta` (`f32`): Specifies the scalar `beta`.
///
/// `c` (`[*]f32`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn sgemm(
    order: Order,
    transa: Transpose,
    transb: Transpose,
    m: isize,
    n: isize,
    k: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    b: [*]const f32,
    ldb: isize,
    beta: f32,
    c: [*]f32,
    ldc: isize,
) void {
    return gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product with general matrices.
///
/// The `dgemm` routine computes a scalar-matrix-matrix product and adds the
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
/// `m` (`isize`): Specifies the number of rows of the matrix `op(A)` and of the
/// matrix `C`. Must be greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `op(B)` and the
/// number of columns of the matrix `C`. Must be greater than or equal to 0.
///
/// `k` (`isize`): Specifies the number of columns of the matrix `op(A)` and the
/// number of rows of the matrix `op(B)`. Must be greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `lda * k`                                       | `lda * m`                                 |
/// | `order = row_major` | `lda * m`                                       | `lda * k`                                 |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, m)`                                     | `max(1, k)`                               |
/// | `order = row_major` | `max(1, k)`                                     | `max(1, m)`                               |
///
/// `b` (`[*]const f64`): Array, size at least:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `ldb * n`                                       | `ldb * k`                                 |
/// | `order = row_major` | `ldb * k`                                       | `ldb * n`                                 |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, k)`                                     | `max(1, n)`                               |
/// | `order = row_major` | `max(1, n)`                                     | `max(1, k)`                               |
///
/// `beta` (`f64`): Specifies the scalar `beta`.
///
/// `c` (`[*]f64`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dgemm(
    order: Order,
    transa: Transpose,
    transb: Transpose,
    m: isize,
    n: isize,
    k: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    b: [*]const f64,
    ldb: isize,
    beta: f64,
    c: [*]f64,
    ldc: isize,
) void {
    return gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product with general matrices.
///
/// The `cgemm` routine computes a scalar-matrix-matrix product and adds the
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
/// `m` (`isize`): Specifies the number of rows of the matrix `op(A)` and of the
/// matrix `C`. Must be greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `op(B)` and the
/// number of columns of the matrix `C`. Must be greater than or equal to 0.
///
/// `k` (`isize`): Specifies the number of columns of the matrix `op(A)` and the
/// number of rows of the matrix `op(B)`. Must be greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `lda * k`                                       | `lda * m`                                 |
/// | `order = row_major` | `lda * m`                                       | `lda * k`                                 |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, m)`                                     | `max(1, k)`                               |
/// | `order = row_major` | `max(1, k)`                                     | `max(1, m)`                               |
///
/// `b` (`[*]const cf32`): Array, size at least:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `ldb * n`                                       | `ldb * k`                                 |
/// | `order = row_major` | `ldb * k`                                       | `ldb * n`                                 |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, k)`                                     | `max(1, n)`                               |
/// | `order = row_major` | `max(1, n)`                                     | `max(1, k)`                               |
///
/// `beta` (`cf32`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf32`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cgemm(
    order: Order,
    transa: Transpose,
    transb: Transpose,
    m: isize,
    n: isize,
    k: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    b: [*]const cf32,
    ldb: isize,
    beta: cf32,
    c: [*]cf32,
    ldc: isize,
) void {
    return gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product with general matrices.
///
/// The `zgemm` routine computes a scalar-matrix-matrix product and adds the
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
/// `m` (`isize`): Specifies the number of rows of the matrix `op(A)` and of the
/// matrix `C`. Must be greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `op(B)` and the
/// number of columns of the matrix `C`. Must be greater than or equal to 0.
///
/// `k` (`isize`): Specifies the number of columns of the matrix `op(A)` and the
/// number of rows of the matrix `op(B)`. Must be greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `lda * k`                                       | `lda * m`                                 |
/// | `order = row_major` | `lda * m`                                       | `lda * k`                                 |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program:
/// |                     | `transa = no_trans` or `transa = conj_no_trans` | `transa = trans` or `transa = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, m)`                                     | `max(1, k)`                               |
/// | `order = row_major` | `max(1, k)`                                     | `max(1, m)`                               |
///
/// `b` (`[*]const cf64`): Array, size at least:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `ldb * n`                                       | `ldb * k`                                 |
/// | `order = row_major` | `ldb * k`                                       | `ldb * n`                                 |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program:
/// |                     | `transb = no_trans` or `transb = conj_no_trans` | `transb = trans` or `transb = conj_trans` |
/// |---------------------|-------------------------------------------------|-------------------------------------------|
/// | `order = col_major` | `max(1, k)`                                     | `max(1, n)`                               |
/// | `order = row_major` | `max(1, n)`                                     | `max(1, k)`                               |
///
/// `beta` (`cf64`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf64`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zgemm(
    order: Order,
    transa: Transpose,
    transb: Transpose,
    m: isize,
    n: isize,
    k: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    b: [*]const cf64,
    ldb: isize,
    beta: cf64,
    c: [*]cf64,
    ldc: isize,
) void {
    return gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

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
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * ka`, where
/// `ka` is `m` if `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `ldb * n` if
/// `order = col_major` and `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
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
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
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
pub inline fn hemm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
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
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_chemm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zhemm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/hemm.zig").hemm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

/// Computes a matrix-matrix product where one input matrix is Hermitian.
///
/// The `chemm` routines compute a scalar-matrix-matrix product using a Hermitian
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
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * ka`, where `ka` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (`[*]const cf32`): Array, size at least `ldb * n` if `order = col_major`
/// and `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`cf32`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf32`): Array, size at least `ldc * n` if `order = col_major` and
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn chemm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    b: [*]const cf32,
    ldb: isize,
    beta: cf32,
    c: [*]cf32,
    ldc: isize,
) void {
    return hemm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is Hermitian.
///
/// The `zhemm` routines compute a scalar-matrix-matrix product using a Hermitian
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
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * ka`, where `ka` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (`[*]const cf64`): Array, size at least `ldb * n` if `order = col_major`
/// and `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`cf64`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf64`): Array, size at least `ldc * n` if `order = col_major` and
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zhemm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    b: [*]const cf64,
    ldb: isize,
    beta: cf64,
    c: [*]cf64,
    ldc: isize,
) void {
    return hemm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Performs a Hermitian rank-`k` update.
///
/// The `herk` routines perform a rank-`k` matrix-matrix operation using a
/// general matrix `A` and a Hermitian matrix `C`. The operation is defined as:
///
/// ```zig
///     C = alpha * A * A^H + beta * C
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^H * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` Hermitian matrix,
/// `A` is an `n`-by-`k` matrix in the first case and a `k`-by-`n` matrix in the
/// second case.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * A^H + beta * C`.
/// - If `trans = conj_trans`, then `C = alpha * A^H * A + beta * C`.
///
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `integer`, `rational`, `real` or
/// `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`             |
/// | `order = row_major` | `lda * n`           | `lda * k`             |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `beta` (`bool`, `int`, `float`, `integer`, `rational`, `real`, or
/// `expression`): Specifies the scalar `beta`.
///
/// `c` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than 0, if `lda` is less
/// than `max(1, n)` or `max(1, k)`, or if `ldc` is less than `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn herk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.herk requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (types.isComplex(Al))
        @compileError("zml.linalg.blas.herk does not support complex alpha, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.herk requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.herk requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.herk requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (types.isComplex(Be))
        @compileError("zml.linalg.blas.herk does not support complex beta, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.herk requires c to be a mutable many-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.herk requires c's child type to be numeric, got " ++ @typeName(C));

    comptime if (Al == bool and A == bool and Be == bool and C == bool)
        @compileError("zml.linalg.blas.herk does not support alpha, a, b, beta and c all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(C))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.herk not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == C and types.canCoerce(Al, Scalar(A)) and types.canCoerce(Be, Scalar(A)) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return ci.cblas_cherk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), scast(Scalar(A), alpha), a, scast(c_int, lda), scast(Scalar(A), beta), c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    return ci.cblas_zherk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), scast(Scalar(A), alpha), a, scast(c_int, lda), scast(Scalar(A), beta), c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/herk.zig").herk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, ctx);
}

/// Performs a Hermitian rank-`k` update.
///
/// The `cherk` routines perform a rank-`k` matrix-matrix operation using a
/// general matrix `A` and a Hermitian matrix `C`. The operation is defined as:
///
/// ```zig
///     C = alpha * A * A^H + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^H * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` Hermitian matrix,
/// `A` is an `n`-by-`k` matrix in the first case and a `k`-by-`n` matrix in the
/// second case.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * A^H + beta * C`.
/// - If `trans = conj_trans`, then `C = alpha * A^H * A + beta * C`.
///
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`             |
/// | `order = row_major` | `lda * n`           | `lda * k`             |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `beta` (`f32`): Specifies the scalar `beta`.
///
/// `c` ([*]cf32`): Array, size at least
/// `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cherk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: f32,
    a: [*]const cf32,
    lda: isize,
    beta: f32,
    c: [*]cf32,
    ldc: isize,
) void {
    return herk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch {};
}

/// Performs a Hermitian rank-`k` update.
///
/// The `zherk` routines perform a rank-`k` matrix-matrix operation using a
/// general matrix `A` and a Hermitian matrix `C`. The operation is defined as:
///
/// ```zig
///     C = alpha * A * A^H + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^H * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` Hermitian matrix,
/// `A` is an `n`-by-`k` matrix in the first case and a `k`-by-`n` matrix in the
/// second case.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * A^H + beta * C`.
/// - If `trans = conj_trans`, then `C = alpha * A^H * A + beta * C`.
///
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`             |
/// | `order = row_major` | `lda * n`           | `lda * k`             |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `beta` (`f64`): Specifies the scalar `beta`.
///
/// `c` ([*]cf64`): Array, size at least
/// `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zherk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: f64,
    a: [*]const cf64,
    lda: isize,
    beta: f64,
    c: [*]cf64,
    ldc: isize,
) void {
    return herk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch {};
}

/// Performs a Hermitian rank-`2k` update.
///
/// The `her2k` routine performs a rank-`2k` matrix-matrix operation using
/// general matrices `A` and `B` and a Hermitian matrix `C`. The operation is
/// defined as:
///
/// ```zig
///     C = alpha * A * B^H + conj(alpha) * B * A^H + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^H * B + conj(alpha) * B^H * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` Hermitian
/// matrix, `A` and `B` are `n`-by-`k` matrices in the first case and `k`-by-`n`
/// matrices in the second case.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `C` is used:
/// - If `uplo = upper`, then the upper triangular part of `C` is used.
/// - If `uplo = lower`, then the lower triangular part of `C` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * B^H + conj(alpha) * B * A^H + beta * C`.
/// - If `trans = conj_trans`, then `C = alpha * A^H * B + conj(alpha) * B^H * A + beta * C`.
///
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = conj_trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`             |
/// | `order = row_major` | `lda * n`           | `lda * k`             |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`             |
/// | `order = row_major` | `ldb * n`           | `ldb * k`             |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `beta` (`bool`, `int`, `float`, `integer`, `rational`, `real` or
/// `expression`): Specifies the scalar `beta`.
///
/// `c` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
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
pub inline fn her2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.her2k requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.her2k requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.her2k requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B))
        @compileError("zml.linalg.blas.her2k requires b to be a many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.her2k requires b's child type to numeric, got " ++ @typeName(B));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.her2k requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (types.isComplex(Be))
        @compileError("zml.linalg.blas.her2k does not support complex beta, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.her2k requires c to be a mutable many-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.her2k requires c's child type to be numeric, got " ++ @typeName(C));

    comptime if (Al == bool and A == bool and B == bool and Be == bool and C == bool)
        @compileError("zml.linalg.blas.her2k does not support alpha, a, b, beta and c all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(C))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.her2k not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, Scalar(A)) and opts.link_cblas != null) {
        switch (comptime types.numericType(Al)) {
            .cfloat => {
                if (comptime Scalar(Al) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_cher2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), scast(Scalar(A), beta), c, scast(c_int, ldc));
                } else if (comptime Scalar(Al) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_zher2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), scast(Scalar(A), beta), c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/her2k.zig").her2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

/// Performs a Hermitian rank-`2k` update.
///
/// The `cher2k` routine performs a rank-`2k` matrix-matrix operation using
/// general matrices `A` and `B` and a Hermitian matrix `C`. The operation is
/// defined as:
///
/// ```zig
///     C = alpha * A * B^H + conj(alpha) * B * A^H + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^H * B + conj(alpha) * B^H * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` Hermitian
/// matrix, `A` and `B` are `n`-by-`k` matrices in the first case and `k`-by-`n`
/// matrices in the second case.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `C` is used:
/// - If `uplo = upper`, then the upper triangular part of `C` is used.
/// - If `uplo = lower`, then the lower triangular part of `C` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * B^H + conj(alpha) * B * A^H + beta * C`.
/// - If `trans = conj_trans`, then `C = alpha * A^H * B + conj(alpha) * B^H * A + beta * C`.
///
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = conj_trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`             |
/// | `order = row_major` | `lda * n`           | `lda * k`             |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `b` (`[*]const cf32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`             |
/// | `order = row_major` | `ldb * n`           | `ldb * k`             |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `beta` (`f32`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf32`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn cher2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    b: [*]const cf32,
    ldb: isize,
    beta: f32,
    c: [*]cf32,
    ldc: isize,
) void {
    return her2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Performs a Hermitian rank-`2k` update.
///
/// The `zher2k` routine performs a rank-`2k` matrix-matrix operation using
/// general matrices `A` and `B` and a Hermitian matrix `C`. The operation is
/// defined as:
///
/// ```zig
///     C = alpha * A * B^H + conj(alpha) * B * A^H + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^H * B + conj(alpha) * B^H * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` Hermitian
/// matrix, `A` and `B` are `n`-by-`k` matrices in the first case and `k`-by-`n`
/// matrices in the second case.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// Hermitian matrix `C` is used:
/// - If `uplo = upper`, then the upper triangular part of `C` is used.
/// - If `uplo = lower`, then the lower triangular part of `C` is used.
///
/// `trans` (`Transpose`): Specifies the operation:
/// - If `trans = no_trans`, then `C = alpha * A * B^H + conj(alpha) * B * A^H + beta * C`.
/// - If `trans = conj_trans`, then `C = alpha * A^H * B + conj(alpha) * B^H * A + beta * C`.
///
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = conj_trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`             |
/// | `order = row_major` | `lda * n`           | `lda * k`             |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `b` (`[*]const cf64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`             |
/// | `order = row_major` | `ldb * n`           | `ldb * k`             |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = conj_trans` |
/// |---------------------|---------------------|-----------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`           |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`           |
///
/// `beta` (`f64`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf64`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zher2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    b: [*]const cf64,
    ldb: isize,
    beta: f64,
    c: [*]cf64,
    ldc: isize,
) void {
    return her2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is symmetric.
///
/// The `symm` routine computes a scalar-matrix-matrix product with one
/// symmetric matrix and adds the result to a scalar-matrix product. The
/// operation is defined as:
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
/// where `alpha` and `beta` are scalars, `A` is a symmetric matrix, `B` and `C`
/// are `m`-by-`n` matrices.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `side` (`Side`): Specifies whether the symmetric matrix `A` appears on the
/// left or right in the operation:
/// - If `side = left`, then `C = alpha * A * B + beta * C`.
/// - If `side = right`, then `C = alpha * B * A + beta * C`.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * ka`, where
/// `ka` is `m` if `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `ldb * n` if
/// `order = col_major` or `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `beta`.
///
/// `c` (mutable many-item pointer to `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `ldc * n` if `order = col_major` or `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
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
/// than `max(1, n)` or `max(1, m)`, if `ldb` is less than `max(1, m)` or
/// `max(1, n)`, or if `ldc` is less than `max(1, n)` or `max(1, m)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn symm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);
    const Be: type = @TypeOf(beta);
    comptime var C: type = @TypeOf(c);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.symm requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.symm requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.symm requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B))
        @compileError("zml.linalg.blas.symm requires b to be a many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.symm requires b's child type to numeric, got " ++ @typeName(B));

    comptime if (!types.isNumeric(Be))
        @compileError("zml.linalg.blas.symm requires beta to be numeric, got " ++ @typeName(Be));

    comptime if (!types.isManyPointer(C) or types.isConstPointer(C))
        @compileError("zml.linalg.blas.symm requires c to be a mutable many-item pointer, got " ++ @typeName(C));

    C = types.Child(C);

    comptime if (!types.isNumeric(C))
        @compileError("zml.linalg.blas.symm requires c's child type to be numeric, got " ++ @typeName(C));

    comptime if (Al == bool and A == bool and B == bool and Be == bool and C == bool)
        @compileError("zml.linalg.blas.symm does not support alpha, a, b, beta and c all being bool");

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B) or
        types.isArbitraryPrecision(Be) or
        types.isArbitraryPrecision(C))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.symm not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                } else if (comptime A == f64) {
                    return ci.cblas_dsymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_csymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zsymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/symm.zig").symm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

/// Computes a matrix-matrix product where one input matrix is symmetric.
///
/// The `ssymm` routine computes a scalar-matrix-matrix product with one
/// symmetric matrix and adds the result to a scalar-matrix product. The
/// operation is defined as:
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
/// where `alpha` and `beta` are scalars, `A` is a symmetric matrix, `B` and `C`
/// are `m`-by-`n` matrices.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `side` (`Side`): Specifies whether the symmetric matrix `A` appears on the
/// left or right in the operation:
/// - If `side = left`, then `C = alpha * A * B + beta * C`.
/// - If `side = right`, then `C = alpha * B * A + beta * C`.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least `lda * ka`, where `ka` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (`[*]const f32`): Array, size at least `ldb * n` if `order = col_major`
/// or `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`f32`): Specifies the scalar `beta`.
///
/// `c` (`[*]f32`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ssymm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    b: [*]const f32,
    ldb: isize,
    beta: f32,
    c: [*]f32,
    ldc: isize,
) void {
    return symm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is symmetric.
///
/// The `dsymm` routine computes a scalar-matrix-matrix product with one
/// symmetric matrix and adds the result to a scalar-matrix product. The
/// operation is defined as:
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
/// where `alpha` and `beta` are scalars, `A` is a symmetric matrix, `B` and `C`
/// are `m`-by-`n` matrices.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `side` (`Side`): Specifies whether the symmetric matrix `A` appears on the
/// left or right in the operation:
/// - If `side = left`, then `C = alpha * A * B + beta * C`.
/// - If `side = right`, then `C = alpha * B * A + beta * C`.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least `lda * ka`, where `ka` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (`[*]const f64`): Array, size at least `ldb * n` if `order = col_major`
/// or `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`f64`): Specifies the scalar `beta`.
///
/// `c` (`[*]f64`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dsymm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    b: [*]const f64,
    ldb: isize,
    beta: f64,
    c: [*]f64,
    ldc: isize,
) void {
    return symm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is symmetric.
///
/// The `csymm` routine computes a scalar-matrix-matrix product with one
/// symmetric matrix and adds the result to a scalar-matrix product. The
/// operation is defined as:
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
/// where `alpha` and `beta` are scalars, `A` is a symmetric matrix, `B` and `C`
/// are `m`-by-`n` matrices.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `side` (`Side`): Specifies whether the symmetric matrix `A` appears on the
/// left or right in the operation:
/// - If `side = left`, then `C = alpha * A * B + beta * C`.
/// - If `side = right`, then `C = alpha * B * A + beta * C`.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * ka`, where `ka` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (`[*]const cf32`): Array, size at least `ldb * n` if `order = col_major`
/// or `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`cf32`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf32`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn csymm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    b: [*]const cf32,
    ldb: isize,
    beta: cf32,
    c: [*]cf32,
    ldc: isize,
) void {
    return symm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is symmetric.
///
/// The `zsymm` routine computes a scalar-matrix-matrix product with one
/// symmetric matrix and adds the result to a scalar-matrix product. The
/// operation is defined as:
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
/// where `alpha` and `beta` are scalars, `A` is a symmetric matrix, `B` and `C`
/// are `m`-by-`n` matrices.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `side` (`Side`): Specifies whether the symmetric matrix `A` appears on the
/// left or right in the operation:
/// - If `side = left`, then `C = alpha * A * B + beta * C`.
/// - If `side = right`, then `C = alpha * B * A + beta * C`.
///
/// `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of the
/// symmetric matrix `A` is used:
/// - If `uplo = upper`, then the upper triangular part of `A` is used.
/// - If `uplo = lower`, then the lower triangular part of `A` is used.
///
/// `m` (`isize`): Specifies the number of rows of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `C`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * ka`, where `ka` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` and `max(1, n)` if `side = right`.
///
/// `b` (`[*]const cf64`): Array, size at least `ldb * n` if `order = col_major`
/// or `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// `beta` (`cf64`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf64`): Array, size at least `ldc * n` if `order = col_major` or
/// `ldc * m` if `order = row_major`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` and `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zsymm(
    order: Order,
    side: Side,
    uplo: Uplo,
    m: isize,
    n: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    b: [*]const cf64,
    ldb: isize,
    beta: cf64,
    c: [*]cf64,
    ldc: isize,
) void {
    return symm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
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
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
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
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
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
pub inline fn syrk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
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
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), scast(A, beta), c, scast(c_int, ldc));
                } else if (comptime A == f64) {
                    return ci.cblas_dsyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), scast(A, beta), c, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_csyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zsyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/syrk.zig").syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, ctx);
}

/// Performs a symmetric rank-`k` update.
///
/// The `ssyrk` routine performs a rank-`k` matrix-matrix operation for a
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`f32`): Specifies the scalar `beta`.
///
/// `c` (`[*]f32`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ssyrk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    beta: f32,
    c: [*]f32,
    ldc: isize,
) void {
    return syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch {};
}

/// Performs a symmetric rank-`k` update.
///
/// The `dsyrk` routine performs a rank-`k` matrix-matrix operation for a
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`f64`): Specifies the scalar `beta`.
///
/// `c` (`[*]f64`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dsyrk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    beta: f64,
    c: [*]f64,
    ldc: isize,
) void {
    return syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch {};
}

/// Performs a symmetric rank-`k` update.
///
/// The `csyrk` routine performs a rank-`k` matrix-matrix operation for a
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`cf32`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf32`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn csyrk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    beta: cf32,
    c: [*]cf32,
    ldc: isize,
) void {
    return syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch {};
}

/// Performs a symmetric rank-`k` update.
///
/// The `zsyrk` routine performs a rank-`k` matrix-matrix operation for a
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`cf64`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf64`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zsyrk(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    beta: cf64,
    c: [*]cf64,
    ldc: isize,
) void {
    return syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch {};
}

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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
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
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
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
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
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
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
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
pub inline fn syr2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    beta: anytype,
    c: anytype,
    ldc: isize,
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
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and types.canCoerce(Al, A) and types.canCoerce(Be, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_ssyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                } else if (comptime A == f64) {
                    return ci.cblas_dsyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb), scast(A, beta), c, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_csyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    const beta_casted: A = scast(A, beta);
                    return ci.cblas_zsyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb), &beta_casted, c, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr2k.zig").syr2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, ctx);
}

/// Performs a symmetric rank-`2k` update.
///
/// The `ssyr2k` routine performs a rank-`2k` matrix-matrix operation for a
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `b` (`[*]const f32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`        |
/// | `order = row_major` | `ldb * n`           | `ldb * k`        |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`f32`): Specifies the scalar `beta`.
///
/// `c` (`[*]f32`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ssyr2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    b: [*]const f32,
    ldb: isize,
    beta: f32,
    c: [*]f32,
    ldc: isize,
) void {
    return syr2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Performs a symmetric rank-`2k` update.
///
/// The `dsyr2k` routine performs a rank-`2k` matrix-matrix operation for a
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `b` (`[*]const f64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`        |
/// | `order = row_major` | `ldb * n`           | `ldb * k`        |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`f64`): Specifies the scalar `beta`.
///
/// `c` (`[*]f64`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dsyr2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    b: [*]const f64,
    ldb: isize,
    beta: f64,
    c: [*]f64,
    ldc: isize,
) void {
    return syr2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Performs a symmetric rank-`2k` update.
///
/// The `csyr2k` routine performs a rank-`2k` matrix-matrix operation for a
/// symmetric matrix `C` using general matrices `A` and `B`. The operation is
/// defined as:
///
/// ```zig
///     C = alpha * A * B^T + conj(alpha) * B * A^T + beta * C,
/// ```
///
/// or
///
/// ```zig
///     C = alpha * A^T * B + conj(alpha) * B^T * A + beta * C,
/// ```
///
/// where `alpha` and `beta` are scalars, `C` is an `n`-by-`n` symmetric
/// matrix, `A` and `B` are `n`-by-`k` matrices in the first case and `k`-by-`n`
/// matrices in the second case.
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `b` (`[*]const cf32`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`        |
/// | `order = row_major` | `ldb * n`           | `ldb * k`        |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`cf32`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf32`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn csyr2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    b: [*]const cf32,
    ldb: isize,
    beta: cf32,
    c: [*]cf32,
    ldc: isize,
) void {
    return syr2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Performs a symmetric rank-`2k` update.
///
/// The `zsyr2k` routine performs a rank-`2k` matrix-matrix operation for a
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
/// `n` (`isize`): Specifies the order of the matrix `C`. Must be greater than
/// or equal to 0.
///
/// `k` (`isize`): With `trans = no_trans`, specifies the number of columns of
/// the matrices `A` and `B`. With `trans = trans`, specifies the number
/// of rows of the matrices `A` and `B`. Must be greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `lda * k`           | `lda * n`        |
/// | `order = row_major` | `lda * n`           | `lda * k`        |
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `b` (`[*]const cf64`): Array, size at least:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `ldb * k`           | `ldb * n`        |
/// | `order = row_major` | `ldb * n`           | `ldb * k`        |
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to:
/// |                     | `transa = no_trans` | `transa = trans` |
/// |---------------------|---------------------|------------------|
/// | `order = col_major` | `max(1, n)`         | `max(1, k)`      |
/// | `order = row_major` | `max(1, k)`         | `max(1, n)`      |
///
/// `beta` (`cf64`): Specifies the scalar `beta`.
///
/// `c` (`[*]cf64`): Array, size at least `ldc * n`.
///
/// `ldc` (`isize`): Specifies the leading dimension of `c` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `void`: The result is stored in `c`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn zsyr2k(
    order: Order,
    uplo: Uplo,
    trans: Transpose,
    n: isize,
    k: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    b: [*]const cf64,
    ldb: isize,
    beta: cf64,
    c: [*]cf64,
    ldc: isize,
) void {
    return syr2k(order, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is triangular.
///
/// The `trmm` routines compute a scalar-matrix-matrix product with one
/// triangular matrix . The operation is defined as:
///
/// ```zig
///     B = alpha * op(A) * B,
/// ```
///
/// or
///
/// ```zig
///     B = alpha * B * op(A),
/// ```
///
/// where `op(X)` is `X`, `X^T`, `conj(X)`, or `X^H`, `alpha` is a scalar,
/// `A` is a unit, or non-unit, upper or lower triangular matrix, and `B` is an
/// `m`-by-`n` matrix.
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * k`, where
/// `k` is `m` if `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `ldb * n` if
/// `order = col_major` or `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
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
/// than `max(1, m)` or `max(1, n)`, or if `ldb` is less than `max(1, m)` or
/// `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available. In that case, no errors will be
/// raised even if the arguments are invalid.
pub inline fn trmm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);

    comptime if (!types.isNumeric(Al))
        @compileError("zml.linalg.blas.trmm requires alpha to be numeric, got " ++ @typeName(Al));

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.blas.trmm requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.blas.trmm requires a's child type to numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B) or types.isConstPointer(B))
        @compileError("zml.linalg.blas.trmm requires b to be a mutable many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.blas.trmm requires b's child type to numeric, got " ++ @typeName(B));

    comptime if (types.isArbitraryPrecision(Al) or
        types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B))
    {
        // When implemented, expand if
        @compileError("zml.linalg.blas.trmm not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_strmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb));
                } else if (comptime A == f64) {
                    return ci.cblas_dtrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_ctrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_ztrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb));
                }
            },
            else => {},
        }
    }

    return @import("blas/trmm.zig").trmm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, ctx);
}

/// Computes a matrix-matrix product where one input matrix is triangular.
///
/// The `strmm` routines compute a scalar-matrix-matrix product with one
/// triangular matrix . The operation is defined as:
///
/// ```zig
///     B = alpha * op(A) * B,
/// ```
///
/// or
///
/// ```zig
///     B = alpha * B * op(A),
/// ```
///
/// where `op(X)` is `X`, `X^T`, `conj(X)`, or `X^H`, `alpha` is a scalar,
/// `A` is a unit, or non-unit, upper or lower triangular matrix, and `B` is an
/// `m`-by-`n` matrix.
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`[*]f32`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn strmm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    b: [*]f32,
    ldb: isize,
) void {
    return trmm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is triangular.
///
/// The `dtrmm` routines compute a scalar-matrix-matrix product with one
/// triangular matrix . The operation is defined as:
///
/// ```zig
///     B = alpha * op(A) * B,
/// ```
///
/// or
///
/// ```zig
///     B = alpha * B * op(A),
/// ```
///
/// where `op(X)` is `X`, `X^T`, `conj(X)`, or `X^H`, `alpha` is a scalar,
/// `A` is a unit, or non-unit, upper or lower triangular matrix, and `B` is an
/// `m`-by-`n` matrix.
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`[*]f64`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtrmm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    b: [*]f64,
    ldb: isize,
) void {
    return trmm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is triangular.
///
/// The `ctrmm` routines compute a scalar-matrix-matrix product with one
/// triangular matrix . The operation is defined as:
///
/// ```zig
///     B = alpha * op(A) * B,
/// ```
///
/// or
///
/// ```zig
///     B = alpha * B * op(A),
/// ```
///
/// where `op(X)` is `X`, `X^T`, `conj(X)`, or `X^H`, `alpha` is a scalar,
/// `A` is a unit, or non-unit, upper or lower triangular matrix, and `B` is an
/// `m`-by-`n` matrix.
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`[*]cf32`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctrmm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    b: [*]cf32,
    ldb: isize,
) void {
    return trmm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

/// Computes a matrix-matrix product where one input matrix is triangular.
///
/// The `ztrmm` routines compute a scalar-matrix-matrix product with one
/// triangular matrix . The operation is defined as:
///
/// ```zig
///     B = alpha * op(A) * B,
/// ```
///
/// or
///
/// ```zig
///     B = alpha * B * op(A),
/// ```
///
/// where `op(X)` is `X`, `X^T`, `conj(X)`, or `X^H`, `alpha` is a scalar,
/// `A` is a unit, or non-unit, upper or lower triangular matrix, and `B` is an
/// `m`-by-`n` matrix.
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`[*]cf64`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztrmm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    b: [*]cf64,
    ldb: isize,
) void {
    return trmm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `a` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `lda * k`, where
/// `k` is `m` if `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (many-item pointer to `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least `ldb * n` if
/// `order = col_major` or `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
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
pub inline fn trsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: anytype,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
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
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and types.canCoerce(Al, A) and opts.link_cblas != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return ci.cblas_strsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb));
                } else if (comptime A == f64) {
                    return ci.cblas_dtrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), scast(A, alpha), a, scast(c_int, lda), b, scast(c_int, ldb));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_ctrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb));
                } else if (comptime Scalar(A) == f64) {
                    const alpha_casted: A = scast(A, alpha);
                    return ci.cblas_ztrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transa), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha_casted, a, scast(c_int, lda), b, scast(c_int, ldb));
                }
            },
            else => {},
        }
    }

    return @import("blas/trsm.zig").trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, ctx);
}

/// Solves a triangular matrix equation.
///
/// The `strsm` routine solves one of the following matrix equations:
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f32`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`f32`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn strsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: f32,
    a: [*]const f32,
    lda: isize,
    b: [*]f32,
    ldb: isize,
) void {
    return trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

/// Solves a triangular matrix equation.
///
/// The `dtrsm` routine solves one of the following matrix equations:
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const f64`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`f64`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn dtrsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: f64,
    a: [*]const f64,
    lda: isize,
    b: [*]f64,
    ldb: isize,
) void {
    return trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

/// Solves a triangular matrix equation.
///
/// The `ctrsm` routine solves one of the following matrix equations:
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`cf32`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ctrsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: cf32,
    a: [*]const cf32,
    lda: isize,
    b: [*]cf32,
    ldb: isize,
) void {
    return trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

/// Solves a triangular matrix equation.
///
/// The `ztrsm` routine solves one of the following matrix equations:
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
/// `m` (`isize`): Specifies the number of rows of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `n` (`isize`): Specifies the number of columns of the matrix `B`. Must be
/// greater than or equal to 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * k`, where `k` is `m` if
/// `side = left` and `n` if `side = right`.
///
/// `lda` (`isize`): Specifies the leading dimension of `a` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `side = left` or `max(1, n)` if `side = right`.
///
/// `b` (`cf64`): Array, size at least `ldb * n` if `order = col_major` or
/// `ldb * m` if `order = row_major`.
///
/// `ldb` (`isize`): Specifies the leading dimension of `b` as declared in the
/// calling (sub)program. Must be greater than or equal to `max(1, m)` if
/// `order = col_major` or `max(1, n)` if `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub inline fn ztrsm(
    order: Order,
    side: Side,
    uplo: Uplo,
    transa: Transpose,
    diag: Diag,
    m: isize,
    n: isize,
    alpha: cf64,
    a: [*]const cf64,
    lda: isize,
    b: [*]cf64,
    ldb: isize,
) void {
    return trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch {};
}

pub const Error = error{
    InvalidArgument,
};
