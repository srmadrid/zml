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

const ci = @import("../c.zig");

pub const Order = enum(c_uint) {
    RowMajor = 101,
    ColumnMajor = 102,
};

pub const Transpose = enum(c_uint) {
    NoTrans = 111,
    Trans = 112,
    ConjTrans = 113,
    ConjNoTrans = 114,
};

pub const Uplo = enum(c_uint) {
    Upper = 121,
    Lower = 122,
};

pub const Diag = enum(c_uint) {
    NonUnit = 131,
    Unit = 132,
};

pub const Side = enum(c_uint) {
    Left = 141,
    Right = 142,
};

// Level 1 BLAS

/// Computes the sum of magnitudes of the vector elements.
///
/// The `asum_sub` routine computes the sum of the magnitudes of elements of a
/// real vector, or the sum of magnitudes of the real and imaginary parts of
/// elements of a complex vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n-1].re) + abs(x[n-1].im)
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer of `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// ret (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to the
/// sum of magnitudes of real and imaginary parts of all elements of the vector.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is a `bool` or not a numeric type, if the type of `ret` is not
/// a mutable one-item pointer, or if the child type of `ret` is not a numeric
/// type.
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
                if (X == f32) {
                    try ops.set(ret, ci.cblas_sasum(scast(c_int, n), x, scast(c_int, incx)), ctx);
                } else if (X == f64) {
                    try ops.set(ret, ci.cblas_dasum(scast(c_int, n), x, scast(c_int, incx)), ctx);
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    try ops.set(ret, ci.cblas_scasum(scast(c_int, n), x, scast(c_int, incx)), ctx);
                } else if (Scalar(X) == f64) {
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
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n-1])
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// ret (`*f32`): Pointer to the sum of magnitudes of all elements of the
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn sasum_sub(n: isize, x: [*]const f32, incx: isize, ret: *f32) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `dasum_sub` routine computes the sum of the magnitudes of elements of a
/// vector:
///
/// ```zig
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n-1])
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// ret (`*f64`): Pointer to the sum of magnitudes of all elements of the
/// vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn dasum_sub(n: isize, x: [*]const f64, incx: isize, ret: *f64) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `scasum_sub` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n-1].re) + abs(x[n-1].im)
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// ret (`*f32`): Pointer to the sum of magnitudes of the real and imaginary
/// parts of all elements of the vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn scasum_sub(n: isize, x: [*]const cf32, incx: isize, ret: *f32) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `dzasum_sub` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n-1].re) + abs(x[n-1].im)
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for indexing vector `x`. Must be
/// greater than 0.
///
/// ret (`*f64`): Pointer to the sum of magnitudes of the real and imaginary
/// parts of all elements of the vector.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn dzasum_sub(n: isize, x: [*]const cf64, incx: isize, ret: *f64) void {
    return asum_sub(n, x, incx, ret, .{}) catch {};
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `asum` routine computes the sum of the magnitudes of elements of a real
/// vector, or the sum of magnitudes of the real and imaginary parts of elements
/// of a complex vector:
///
/// ```zig
///     abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n-1].re) + abs(x[n-1].im)
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer of `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer or if the
/// child type of `x` is a `bool` or not a numeric type.
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

    comptime if (!types.isManyPointer(X) or types.isConstPointer(X))
        @compileError("zml.linalg.blas.asum requires x to be a mutable many-item pointer, got " ++ @typeName(X));

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
                if (X == f32) {
                    return ci.cblas_sasum(scast(c_int, n), x, scast(c_int, incx));
                } else if (X == f64) {
                    return ci.cblas_dasum(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_scasum(scast(c_int, n), x, scast(c_int, incx));
                } else if (Scalar(X) == f64) {
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
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n-1])
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn sasum(n: isize, x: [*]const f32, incx: isize) f32 {
    return asum(f32, n, x, incx, .{}) catch 0;
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `sasum` routine computes the sum of the magnitudes of elements of a
/// vector:
///
/// ```zig
///     ret = abs(x[0]) + abs(x[1]) + ... + abs(x[n-1])
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn dasum(n: isize, x: [*]const f64, incx: isize) f64 {
    return asum(f64, n, x, incx, .{}) catch 0;
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `scasum` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n-1].re) + abs(x[n-1].im)
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn scasum(n: isize, x: [*]const cf32, incx: isize) f32 {
    return asum(cf32, n, x, incx, .{}) catch 0;
}

/// Computes the sum of magnitudes of the vector elements.
///
/// The `scasum` routine computes the sum of magnitudes of the real and
/// imaginary parts of a vector:
///
/// ```zig
///     ret = abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n-1].re) + abs(x[n-1].im)
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn dzasum(n: isize, x: [*]const cf64, incx: isize) f64 {
    return asum(cf64, n, x, incx, .{}) catch 0;
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `axpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with a number of
/// elements that equals `n`.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `alpha` is not a numeric type, if the type
/// of `x` is not a many-item pointer, if the child type of `x` is not a numeric
/// type, if the type of `y` is not a mutable many-item pointer, if the child
/// type of `y` is not a numeric type, or if `alpha`, `x` and `y` are all
/// `bool`.
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

    if (comptime Al == X and Al == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(Al)) {
            .float => {
                if (Al == f32) {
                    return ci.cblas_saxpy(scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (Al == f64) {
                    return ci.cblas_daxpy(scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (Scalar(Al) == f32) {
                    return ci.cblas_caxpy(scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (Scalar(Al) == f64) {
                    return ci.cblas_zaxpy(scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy));
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
///     y = alpha * x + y
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with a number of
/// elements that equals `n`.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn saxpy(n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]f32, incy: isize) void {
    return axpy(n, alpha, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `daxpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with a number of
/// elements that equals `n`.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn daxpy(n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]f64, incy: isize) void {
    return axpy(n, alpha, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `caxpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with a number of
/// elements that equals `n`.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn caxpy(n: isize, alpha: cf32, x: [*]const cf32, incx: isize, y: [*]cf32, incy: isize) void {
    return axpy(n, alpha, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-scalar product and adds the result to a vector.
///
/// The `caxpy` routine performs a vector-vector operation defined as:
///
/// ```zig
///     y = alpha * x + y
/// ```
///
/// where `alpha` is a scalar, and `x` and `y` are vectors each with a number of
/// elements that equals `n`.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn zaxpy(n: isize, alpha: cf64, x: [*]const cf64, incx: isize, y: [*]cf64, incy: isize) void {
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
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a mutable
/// many-item pointer, or if the child type of `y` is not a numeric type.
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
                if (X == f32) {
                    return ci.cblas_scopy(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (X == f64) {
                    return ci.cblas_dcopy(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_ccopy(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (Scalar(X) == f64) {
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
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
pub fn scopy(n: isize, x: [*]const f32, incx: isize, y: [*]f32, incy: isize) void {
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
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
pub fn dcopy(n: isize, x: [*]const f64, incx: isize, y: [*]f64, incy: isize) void {
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
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
pub fn ccopy(n: isize, x: [*]const cf32, incx: isize, y: [*]cf32, incy: isize) void {
    return copy(n, x, incx, y, incy, .{}) catch {};
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
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
pub fn zcopy(n: isize, x: [*]const cf64, incx: isize, y: [*]cf64, incy: isize) void {
    return copy(n, x, incx, y, incy, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `dot_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incy))`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to the
/// result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, if `x`  and `y` are
/// both `bool`, if the type of `ret` is not a mutable one-item pointer, or if
/// the child type of `ret` is not a numeric type.
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
                if (X == f32) {
                    try ops.set(ret, ci.cblas_sdot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy)), ctx);
                } else if (X == f64) {
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
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*f32`): Pointer to the result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn sdot_sub(n: isize, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize, ret: *f32) void {
    return dot_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `ddot_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*f64`): Pointer to the result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn ddot_sub(n: isize, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize, ret: *f64) void {
    return dot_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `dot` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `integer`, `rational`,
/// `real` or `expression`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer of `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incy))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, or if `x` and `y`
/// are both `bool`.
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

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.dot requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

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
                if (X == f32) {
                    return ci.cblas_sdot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (X == f64) {
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
///     x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn sdot(n: isize, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize) f32 {
    return dot(n, x, incx, y, incy, .{}) catch 0;
}

/// Computes a vector-vector dot product.
///
/// The `ddot` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn ddot(n: isize, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize) f64 {
    return dot(n, x, incx, y, incy, .{}) catch 0;
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `dotc_sub` routine performs a vector-vector operation defined as:
///
/// ```zig
///     ret = conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n-1]) * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer of `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incy))`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to the
/// result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, if `x` and `y` are
/// both `bool`, if the type of `ret` is not a mutable one-item pointer or if
/// the child type of `ret` is not a numeric type.
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
                if (Scalar(X) == f32) {
                    return ci.cblas_cdotc_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), ret);
                } else if (Scalar(X) == f64) {
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
///     ret = conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n-1]) * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*cf32`): Pointer to the result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn cdotc_sub(n: isize, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize, ret: *cf32) void {
    return dotc_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `zdotc_sub` routine performs a vector-vector operation defined as:
///
/// ```zig
///     ret = conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n-1]) * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn zdotc_sub(n: isize, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize, ret: *cf64) void {
    return dotc_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `dotc` routine performs a vector-vector operation defined as:
///
/// ```zig
///     conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n-1]) * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer of `bool`, `int`, `float`, `cfloat` `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incy))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, if `x` and `y` are
/// both `bool`, if the type of `ret` is not a mutable one-item pointer or if
/// the child type of `ret` is not a numeric type.
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

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.dotc not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .cfloat => {
                if (Scalar(X) == f32) {
                    var temp: cf32 = undefined;
                    ci.cblas_cdotc_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), &temp);
                    return temp;
                } else if (Scalar(X) == f64) {
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
///     conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n-1]) * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn cdotc(n: isize, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize) cf32 {
    return dotc(n, x, incx, y, incy, .{}) catch .{ .re = 0, .im = 0 };
}

/// Computes a dot product of a conjugated vector with another vector.
///
/// The `zdotc` routine performs a vector-vector operation defined as:
///
/// ```zig
///     conj(x[0]) * y[0] + conj(x[1]) * y[1] + ... + conj(x[n-1]) * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn zdotc(n: isize, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize) cf64 {
    return dotc(n, x, incx, y, incy, .{}) catch .{ .re = 0, .im = 0 };
}

/// Computes a vector-vector dot product.
///
/// The `dotu_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incy))`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Pointer to the
/// result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Errors
/// ------
/// `linalg.blas.Error.InvalidArgument`: If `n` is less than or equal to 0.
///
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, if `x`  and `y` are
/// both `bool`, if the type of `ret` is not a mutable one-item pointer, or if
/// the child type of `ret` is not a numeric type.
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
                if (Scalar(X) == f32) {
                    return ci.cblas_cdotu_sub(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), ret);
                } else if (Scalar(X) == f64) {
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
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `ret` (`*cf32`): Pointer to the result of the dot product.
///
/// Returns
/// -------
/// `void`: The result is stored in `ret`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding CBLAS function.
pub fn cdotu_sub(n: isize, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize, ret: *cf32) void {
    return dotu_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `zdotu_sub` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     ret = x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn zdotu_sub(n: isize, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize, ret: *cf64) void {
    return dotu_sub(n, x, incx, y, incy, ret, .{}) catch {};
}

/// Computes a vector-vector dot product.
///
/// The `dotu` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incy))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, or if `x` and `y`
/// are both `bool`.
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

    comptime if (!types.isManyPointer(Y) or types.isConstPointer(Y))
        @compileError("zml.linalg.blas.dotu requires y to be a mutable many-item pointer, got " ++ @typeName(Y));

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
            .float => {
                if (X == f32) {
                    return ci.cblas_sdot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (X == f64) {
                    return ci.cblas_ddot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
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
///     x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn cdotu(n: isize, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize) cf32 {
    return dotu(n, x, incx, y, incy, .{}) catch .{ .re = 0, .im = 0 };
}

/// Computes a vector-vector dot product.
///
/// The `zdotu` routine performs a vector-vector reduction operation defined
/// as:
///
/// ```zig
///     x[0] * y[0] + x[1] * y[1] + ... + x[n-1] * y[n-1]
/// ```
///
/// where `x` and `y` are vectors.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn zdotu(n: isize, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize) cf64 {
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
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (many-item pointer of, `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer or if the
/// child type of `x` is a bool or not a numeric type.
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
                if (X == f32) {
                    return ci.cblas_snrm2(scast(c_int, n), x, scast(c_int, incx));
                } else if (X == f64) {
                    return ci.cblas_dnrm2(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_scnrm2(scast(c_int, n), x, scast(c_int, incx));
                } else if (Scalar(X) == f64) {
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
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn snrm2(n: isize, x: [*]const f32, incx: isize) f32 {
    return nrm2(f32, n, x, incx);
}

/// Computes the Euclidean norm of a vector.
///
/// The `dnrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn dnrm2(n: isize, x: [*]const f64, incx: isize) f64 {
    return nrm2(f64, n, x, incx);
}

/// Computes the Euclidean norm of a vector.
///
/// The `scnrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn scnrm2(n: isize, x: [*]const cf32, incx: isize) f32 {
    return nrm2(cf32, n, x, incx);
}

/// Computes the Euclidean norm of a vector.
///
/// The `dznrm2` routine performs a vector reduction operation defined as:
///
/// ```zig
///     ||x||,
/// ```
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vector `x`. Must be
/// greater than 0.
///
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn dznrm2(n: isize, x: [*]const cf64, incx: isize) f64 {
    return nrm2(cf64, n, x, incx);
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
/// `x` (mutable many-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, if `c` is not a
/// numeric type, if `c` is complex, if `s` is not a numeric type, if `s` is
/// complex, or if `x`, `y`, `c` and `s` are all `bool`.
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
                if (X == f32) {
                    return ci.cblas_srot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), c, s);
                } else if (X == f64) {
                    return ci.cblas_drot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), c, s);
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_csrot(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), c, s);
                } else if (Scalar(X) == f64) {
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
/// `x` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn srot(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, c: f32, s: f32) void {
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
/// `x` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn drot(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, c: f64, s: f64) void {
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
/// `x` (`[*]cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn csrot(n: isize, x: [*]cf32, incx: isize, y: [*]cf32, incy: isize, c: f32, s: f32) void {
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
/// `x` (`[*]cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn zdrot(n: isize, x: [*]cf64, incx: isize, y: [*]cf64, incy: isize, c: f64, s: f64) void {
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
/// `a` (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Provides the
/// `x`-coordinate of the point `p`. On return, it contains the parameter `r`
/// associated with the Givens rotation.
///
/// `b` (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Provides the
/// `y`-coordinate of the point `p`. On return, it contains the parameter `z`
/// associated with the Givens rotation.
///
/// `c` (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): On return, it
/// contains the parameter `c` associated with the Givens rotation.
///
/// `s` (mutable one-item pointer of `bool`, `int`, `float`, `cfloat`,
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
/// Raises
/// ------
/// `@compileError`: If the type of `a`, `b`, `c` or `s` is not a one-item
/// pointer, if the child type of `a`, `b`, `c` or `s` is not a numeric type,
/// if `c` is complex, or if `a`, `b`, `c` and `s` are all `bool`.
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
                if (A == f32) {
                    return ci.cblas_srotg(a, b, c, s);
                } else if (A == f64) {
                    return ci.cblas_drotg(a, b, c, s);
                }
            },
            .cfloat => {
                if (Scalar(A) == f32) {
                    return ci.cblas_crotg(a, b, c, s);
                } else if (Scalar(A) == f64) {
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
pub fn srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
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
pub fn drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
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
pub fn crotg(a: *cf32, b: *cf32, c: *f32, s: *cf32) void {
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
pub fn zrotg(a: *cf64, b: *cf64, c: *f64, s: *cf64) void {
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
/// `x` (mutable many-item pointer of `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`. On return, every element `x[i]` is replaced
/// by `h11 * x[i] + h12 * y[i]`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer of `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incy))`. On return, every element `y[i]` is replaced
/// by `h21 * x[i] + h22 * y[i]`.
///
/// `incy` (`isize`): Specifies the increment for the elements of `y`.
///
/// `param` (many-item pointer of `bool`, `int`, `float`, `integer`, `rational`,
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer, if the child
/// type of `x` is not a numeric type, if the type of `y` is not a many-item
/// pointer, if the child type of `y` is not a numeric type, if the type of
/// `param` is not a many-item pointer, if the child type of `param` is not a
/// numeric type, if `param` is complex, if `x`, `y`, and `param` are all
/// `bool`, or if `x`, `y`, or `param` are complex types.
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
                if (X == f32) {
                    return ci.cblas_srotm(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy), param);
                } else if (X == f64) {
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
/// `x` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return,
/// every element `x[i]` is replaced by `h11 * x[i] + h12 * y[i]`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incy))`. On return,
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
pub fn srotm(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, param: [*]const f32) void {
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
/// `x` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return,
/// every element `x[i]` is replaced by `h11 * x[i] + h12 * y[i]`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incy))`. On return,
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
pub fn drotm(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, param: [*]const f64) void {
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
/// `d1` (mutable one-item pointer of `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Provides the scaling factor for the
/// `x`-coordinate of the input vector. On return it provides the first diagonal
/// element of the updated matrix.
///
/// `d2` (mutable one-item pointer of `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Provides the scaling factor for the
/// `y`-coordinate of the input vector. On return it provides the second diagonal
/// element of the updated matrix.
///
/// `x1` (mutable one-item pointer of `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Provides the `x`-coordinate of the
/// input vector. On return it provides the `x`-coordinate of the rotated vector
/// before scaling.
///
/// `y1` (`bool`, `int`, `float`, `integer`, `rational`, `real` or
/// `expression`): Provides the `y`-coordinate of the input vector.
///
/// `param` (mutable many-item pointer of `bool`, `int`, `float`, `integer`,
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
/// Raises
/// ------
/// `@compileError`: If the type of `d1`, `d2` or `x1` is not a mutable one-item
/// pointer, if the type of `y1` is not numeric, if the type of `param` is not a
/// mutable many-item pointer, if the child type of `d1`, `d2`, `x1` or `param`
/// is not a numeric type, if `d1`, `d2`, `x1`, `y1` or `param` are complex,
/// or if `d1`, `d2`, `x1`, `y1` and `param` are all `bool`.
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

    if (comptime D1 == D2 and D1 == X1 and D1 == P and opts.link_cblas != null) {
        switch (comptime types.numericType(D1)) {
            .float => {
                if (D1 == f32) {
                    return ci.cblas_srotmg(d1, d2, x1, y1, param);
                } else if (D1 == f64) {
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
pub fn srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*]f32) void {
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
pub fn drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*]f64) void {
    return rotmg(f64, d1, d2, x1, y1, param, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): Specifies the scalar `alpha`.
///
/// `x` (many-item pointer of `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`. On return contains the updated vector `x`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `alpha` is not a numeric type, if the type
/// of `x` is not a many-item pointer, if the child type of `x` is not a numeric
/// type, or if `alpha` and `x` are both `bool`.
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

    if (comptime X == Al and opts.link_cblas != null) {
        switch (comptime types.numericType(X)) {
            .float => {
                if (X == f32) {
                    return ci.cblas_sscal(scast(c_int, n), alpha, x, scast(c_int, incx));
                } else if (X == f64) {
                    return ci.cblas_dscal(scast(c_int, n), alpha, x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_cscal(scast(c_int, n), &alpha, x, scast(c_int, incx));
                } else if (Scalar(X) == f64) {
                    return ci.cblas_zscal(scast(c_int, n), &alpha, x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/scal.zig").scal(n, alpha, x, incx, ctx);
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return
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
pub fn sscal(n: isize, alpha: f32, x: [*]f32, incx: isize) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return
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
pub fn dscal(n: isize, alpha: f64, x: [*]f64, incx: isize) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`cf32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return
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
pub fn cscal(n: isize, alpha: cf32, x: [*]cf32, incx: isize) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`cf64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return
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
pub fn zscal(n: isize, alpha: cf64, x: [*]cf64, incx: isize) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f32`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return
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
pub fn csscal(n: isize, alpha: f32, x: [*]cf32, incx: isize) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Computes the product of a vector by a scalar.
///
/// The `scal` routine performs a vector operation defined as:
///
/// ```zig
///     x = alpha * x
/// ```
///
/// where `alpha` is a scalar, and `x` is an `n`-element vector.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `alpha` (`f64`): Specifies the scalar `alpha`.
///
/// `x` (`[*]cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`. On return
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
pub fn zdscal(n: isize, alpha: f64, x: [*]cf64, incx: isize) void {
    return scal(n, alpha, x, incx, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `swap` routines return vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (mutable many-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer of `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `(1 + (n - 1) * abs(incy))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a mutable many-item pointer, if
/// the child type of `x` is not a numeric type, if the type of `y` is not a
/// mutable many-item pointer, or if the child type of `y` is not a numeric
/// type.
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
                if (X == f32) {
                    return ci.cblas_sswap(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (X == f64) {
                    return ci.cblas_dswap(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_cswap(scast(c_int, n), x, scast(c_int, incx), y, scast(c_int, incy));
                } else if (Scalar(X) == f64) {
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
/// Given two vectors `x` and `y`, the `swap` routines return vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn sswap(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `swap` routines return vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]f64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn dswap(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `swap` routines return vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf32`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn cswap(n: isize, x: [*]cf32, incx: isize, y: [*]cf32, incy: isize) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Swaps a vector with another vector.
///
/// Given two vectors `x` and `y`, the `swap` routines return vectors `y` and
/// `x` swapped, each replacing the other.
///
/// Parameters
/// ----------
/// `n` (`isize`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (`[*]cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
///
/// `incx` (`isize`): Specifies the increment for the elements of `x`.
///
/// `y` (`[*]cf64`): Array, size at least `(1 + (n - 1) * abs(incy))`.
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
pub fn zswap(n: isize, x: [*]cf64, incx: isize, y: [*]cf64, incy: isize) void {
    return swap(n, x, incx, y, incy, .{}) catch {};
}

/// Finds the index of the element with maximum absolute value.
///
/// Given a vector `x`, the `iamax` routine returns the position of the vector
/// element `x[i]` that has the largest absolute value for real flavors, or the
/// largest sum `|x[i].re| + |x[i].im|` for complex flavors.
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
/// `x` (many-item pointer of `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer or if the
/// child type of `x` is a `bool` or not a numeric type.
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
                if (X == f32) {
                    return ci.cblas_isamax(scast(c_int, n), x, scast(c_int, incx));
                } else if (X == f64) {
                    return ci.cblas_idamax(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_icamax(scast(c_int, n), x, scast(c_int, incx));
                } else if (Scalar(X) == f64) {
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
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn isamax(n: isize, x: [*]const f32, incx: isize) usize {
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
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn idamax(n: isize, x: [*]const f64, incx: isize) usize {
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
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn icamax(n: isize, x: [*]const cf32, incx: isize) usize {
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
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn izamax(n: isize, x: [*]const cf64, incx: isize) usize {
    return iamax(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `iamin` routine returns the position of the vector
/// element `x[i]` that has the smallest absolute value for real flavors, or the
/// smallest sum `|x[i].re| + |x[i].im|` for complex flavors.
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
/// `x` (many-item pointer of `int`, `float`, `cfloat`, `integer`, `rational`,
/// `real`, `complex` or `expression`): Array, size at least
/// `(1 + (n - 1) * abs(incx))`.
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
/// Raises
/// ------
/// `@compileError`: If the type of `x` is not a many-item pointer or if the
/// child type of `x` is a `bool` or not a numeric type.
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
                if (X == f32) {
                    return ci.cblas_isamax(scast(c_int, n), x, scast(c_int, incx));
                } else if (X == f64) {
                    return ci.cblas_idamax(scast(c_int, n), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(X) == f32) {
                    return ci.cblas_icamax(scast(c_int, n), x, scast(c_int, incx));
                } else if (Scalar(X) == f64) {
                    return ci.cblas_izamax(scast(c_int, n), x, scast(c_int, incx));
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
/// `x` (`[*]const f32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn isamin(n: isize, x: [*]const f32, incx: isize) usize {
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
/// `x` (`[*]const f64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn idamin(n: isize, x: [*]const f64, incx: isize) usize {
    return iamin(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `iamin` routine returns the position of the vector
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
/// `x` (`[*]const cf32`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn icamin(n: isize, x: [*]const cf32, incx: isize) usize {
    return iamin(n, x, incx, .{}) catch 0;
}

/// Finds the index of the element with the smallest absolute value.
///
/// Given a vector `x`, the `iamin` routine returns the position of the vector
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
/// `x` (`[*]const cf64`): Array, size at least `(1 + (n - 1) * abs(incx))`.
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
pub fn izamin(n: isize, x: [*]const cf64, incx: isize) usize {
    return iamin(n, x, incx, .{}) catch 0;
}

// Level 2 BLAS
pub inline fn gbmv(comptime T: type, order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_sgbmv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                } else if (T == f64) {
                    return ci.cblas_dgbmv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgbmv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgbmv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), scast(c_int, kl), scast(c_int, ku), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/gbmv.zig").gbmv(T, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn sgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return gbmv(f32, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return gbmv(f64, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn cgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: cf32, A: [*]const cf32, lda: isize, x: [*]const cf32, incx: isize, beta: cf32, y: [*]cf32, incy: isize) void {
    return gbmv(cf32, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: cf64, A: [*]const cf64, lda: isize, x: [*]const cf64, incx: isize, beta: cf64, y: [*]cf64, incy: isize) void {
    return gbmv(cf64, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn gemv(comptime T: type, order: Order, transA: Transpose, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_sgemv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                } else if (T == f64) {
                    return ci.cblas_dgemv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgemv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgemv(@intFromEnum(order), @intFromEnum(transA), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/gemv.zig").gemv(T, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn sgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return gemv(f32, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return gemv(f64, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn cgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: cf32, A: [*]const cf32, lda: isize, x: [*]const cf32, incx: isize, beta: cf32, y: [*]cf32, incy: isize) void {
    return gemv(cf32, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: cf64, A: [*]const cf64, lda: isize, x: [*]const cf64, incx: isize, beta: cf64, y: [*]cf64, incy: isize) void {
    return gemv(cf64, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn ger(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_sger(@intFromEnum(order), scast(c_int, m), scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                } else if (T == f64) {
                    return ci.cblas_dger(@intFromEnum(order), scast(c_int, m), scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/ger.zig").ger(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn sger(order: Order, m: isize, n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize, A: [*]f32, lda: isize) void {
    return ger(f32, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn dger(order: Order, m: isize, n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize, A: [*]f64, lda: isize) void {
    return ger(f64, order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn gerc(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgerc(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgerc(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/gerc.zig").gerc(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cgerc(order: Order, m: isize, n: isize, alpha: cf32, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize, A: [*]cf32, lda: isize) void {
    return gerc(cf32, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zgerc(order: Order, m: isize, n: isize, alpha: cf64, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize, A: [*]cf64, lda: isize) void {
    return gerc(cf64, order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn geru(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgeru(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgeru(@intFromEnum(order), scast(c_int, m), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/geru.zig").geru(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cgeru(order: Order, m: isize, n: isize, alpha: cf32, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize, A: [*]cf32, lda: isize) void {
    return geru(cf32, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zgeru(order: Order, m: isize, n: isize, alpha: cf64, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize, A: [*]cf64, lda: isize) void {
    return geru(cf64, order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn hbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_chbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                } else if (T == f64) {
                    return ci.cblas_zhbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hbmv.zig").hbmv(T, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn chbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: cf32, A: [*]const cf32, lda: isize, x: [*]const cf32, incx: isize, beta: cf32, y: [*]cf32, incy: isize) void {
    return hbmv(cf32, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zhbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: cf64, A: [*]const cf64, lda: isize, x: [*]const cf64, incx: isize, beta: cf64, y: [*]cf64, incy: isize) void {
    return hbmv(cf64, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn hemv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chemv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhemv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, A, scast(c_int, lda), x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hemv.zig").hemv(T, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn chemv(order: Order, uplo: Uplo, n: isize, alpha: cf32, A: [*]const cf32, lda: isize, x: [*]const cf32, incx: isize, beta: cf32, y: [*]cf32, incy: isize) void {
    return hemv(cf32, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zhemv(order: Order, uplo: Uplo, n: isize, alpha: cf64, A: [*]const cf64, lda: isize, x: [*]const cf64, incx: isize, beta: cf64, y: [*]cf64, incy: isize) void {
    return hemv(cf64, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn her(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cher(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), A, scast(c_int, lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zher(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), A, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/her.zig").her(T, order, uplo, n, alpha, x, incx, A, lda);
}
pub fn cher(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const cf32, incx: isize, A: [*]cf32, lda: isize) void {
    return her(cf32, order, uplo, n, alpha, x, incx, A, lda);
}
pub fn zher(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const cf64, incx: isize, A: [*]cf64, lda: isize) void {
    return her(cf64, order, uplo, n, alpha, x, incx, A, lda);
}

pub inline fn her2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cher2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zher2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/her2.zig").her2(T, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cher2(order: Order, uplo: Uplo, n: isize, alpha: cf32, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize, A: [*]cf32, lda: isize) void {
    return her2(cf32, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zher2(order: Order, uplo: Uplo, n: isize, alpha: cf64, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize, A: [*]cf64, lda: isize) void {
    return her2(cf64, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn hpmv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, Ap: [*]const T, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chpmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, Ap, x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhpmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, Ap, x, scast(c_int, incx), &beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hpmv.zig").hpmv(T, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn chpmv(order: Order, uplo: Uplo, n: isize, alpha: cf32, Ap: [*]const cf32, x: [*]const cf32, incx: isize, beta: cf32, y: [*]cf32, incy: isize) void {
    return hpmv(cf32, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn zhpmv(order: Order, uplo: Uplo, n: isize, alpha: cf64, Ap: [*]const cf64, x: [*]const cf64, incx: isize, beta: cf64, y: [*]cf64, incy: isize) void {
    return hpmv(cf64, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}

pub inline fn hpr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, Ap: [*]T) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chpr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), Ap);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhpr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/hpr.zig").hpr(T, order, uplo, n, alpha, x, incx, Ap);
}
pub fn chpr(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const cf32, incx: isize, Ap: [*]cf32) void {
    return hpr(cf32, order, uplo, n, alpha, x, incx, Ap);
}
pub fn zhpr(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const cf64, incx: isize, Ap: [*]cf64) void {
    return hpr(cf64, order, uplo, n, alpha, x, incx, Ap);
}

pub inline fn hpr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, Ap: [*]T) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chpr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), Ap);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhpr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), &alpha, x, scast(c_int, incx), y, scast(c_int, incy), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/hpr2.zig").hpr2(T, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn chpr2(order: Order, uplo: Uplo, n: isize, alpha: cf32, x: [*]const cf32, incx: isize, y: [*]const cf32, incy: isize, Ap: [*]cf32) void {
    return hpr2(cf32, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn zhpr2(order: Order, uplo: Uplo, n: isize, alpha: cf64, x: [*]const cf64, incx: isize, y: [*]const cf64, incy: isize, Ap: [*]cf64) void {
    return hpr2(cf64, order, uplo, n, alpha, x, incx, y, incy, Ap);
}

pub inline fn sbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_ssbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                } else if (T == f64) {
                    return ci.cblas_dsbmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/sbmv.zig").sbmv(T, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn ssbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return sbmv(f32, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dsbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return sbmv(f64, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn spmv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, Ap: [*]const T, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_sspmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, Ap, x, scast(c_int, incx), beta, y, scast(c_int, incy));
                } else if (T == f64) {
                    return ci.cblas_dspmv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, Ap, x, scast(c_int, incx), beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/spmv.zig").spmv(T, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn sspmv(order: Order, uplo: Uplo, n: isize, alpha: f32, Ap: [*]const f32, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return spmv(f32, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn dspmv(order: Order, uplo: Uplo, n: isize, alpha: f64, Ap: [*]const f64, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return spmv(f64, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}

pub inline fn spr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, Ap: [*]T) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_sspr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), Ap);
                } else if (T == f64) {
                    return ci.cblas_dspr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/spr.zig").spr(T, order, uplo, n, alpha, x, incx, Ap);
}
pub fn sspr(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, Ap: [*]f32) void {
    return spr(f32, order, uplo, n, alpha, x, incx, Ap);
}
pub fn dspr(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, Ap: [*]f64) void {
    return spr(f64, order, uplo, n, alpha, x, incx, Ap);
}

pub inline fn spr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, Ap: [*]T) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_sspr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy), Ap);
                } else if (T == f64) {
                    return ci.cblas_dspr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/spr2.zig").spr2(T, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn sspr2(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize, Ap: [*]f32) void {
    return spr2(f32, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn dspr2(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize, Ap: [*]f64) void {
    return spr2(f64, order, uplo, n, alpha, x, incx, y, incy, Ap);
}

pub inline fn symv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_ssymv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                } else if (T == f64) {
                    return ci.cblas_dsymv(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, A, scast(c_int, lda), x, scast(c_int, incx), beta, y, scast(c_int, incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/symv.zig").symv(T, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn ssymv(order: Order, uplo: Uplo, n: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return symv(f32, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dsymv(order: Order, uplo: Uplo, n: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return symv(f64, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn syr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_ssyr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), A, scast(c_int, lda));
                } else if (T == f64) {
                    return ci.cblas_dsyr(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), A, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr.zig").syr(T, order, uplo, n, alpha, x, incx, A, lda);
}
pub fn ssyr(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, A: [*]f32, lda: isize) void {
    return syr(f32, order, uplo, n, alpha, x, incx, A, lda);
}
pub fn dsyr(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, A: [*]f64, lda: isize) void {
    return syr(f64, order, uplo, n, alpha, x, incx, A, lda);
}

pub inline fn syr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_ssyr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                } else if (T == f64) {
                    return ci.cblas_dsyr2(@intFromEnum(order), @intFromEnum(uplo), scast(c_int, n), alpha, x, scast(c_int, incx), y, scast(c_int, incy), A, scast(c_int, lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr2.zig").syr2(T, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn ssyr2(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize, A: [*]f32, lda: isize) void {
    return syr2(f32, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn dsyr2(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize, A: [*]f64, lda: isize) void {
    return syr2(f64, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn tbmv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_stbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (T == f64) {
                    return ci.cblas_dtbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tbmv.zig").tbmv(T, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn stbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return tbmv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn dtbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return tbmv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ctbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const cf32, lda: isize, x: [*]cf32, incx: isize) void {
    return tbmv(cf32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ztbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const cf64, lda: isize, x: [*]cf64, incx: isize) void {
    return tbmv(cf64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}

pub inline fn tbsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_stbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (T == f64) {
                    return ci.cblas_dtbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), scast(c_int, k), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tbsv.zig").tbsv(T, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn stbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return tbsv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn dtbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return tbsv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ctbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const cf32, lda: isize, x: [*]cf32, incx: isize) void {
    return tbsv(cf32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ztbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const cf64, lda: isize, x: [*]cf64, incx: isize) void {
    return tbsv(cf64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}

pub inline fn tpmv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const T, x: [*]T, incx: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_stpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                } else if (T == f64) {
                    return ci.cblas_dtpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tpmv.zig").tpmv(T, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn stpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f32, x: [*]f32, incx: isize) void {
    return tpmv(f32, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn dtpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f64, x: [*]f64, incx: isize) void {
    return tpmv(f64, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ctpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const cf32, x: [*]cf32, incx: isize) void {
    return tpmv(cf32, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ztpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const cf64, x: [*]cf64, incx: isize) void {
    return tpmv(cf64, order, uplo, transA, diag, n, Ap, x, incx);
}

pub inline fn tpsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const T, x: [*]T, incx: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_stpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                } else if (T == f64) {
                    return ci.cblas_dtpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), Ap, x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tpsv.zig").tpsv(T, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn stpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f32, x: [*]f32, incx: isize) void {
    return tpsv(f32, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn dtpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f64, x: [*]f64, incx: isize) void {
    return tpsv(f64, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ctpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const cf32, x: [*]cf32, incx: isize) void {
    return tpsv(cf32, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ztpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const cf64, x: [*]cf64, incx: isize) void {
    return tpsv(cf64, order, uplo, transA, diag, n, Ap, x, incx);
}

pub inline fn trmv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_strmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (T == f64) {
                    return ci.cblas_dtrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/trmv.zig").trmv(T, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn strmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return trmv(f32, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn dtrmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return trmv(f64, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ctrmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const cf32, lda: isize, x: [*]cf32, incx: isize) void {
    return trmv(cf32, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ztrmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const cf64, lda: isize, x: [*]cf64, incx: isize) void {
    return trmv(cf64, order, uplo, transA, diag, n, A, lda, x, incx);
}

pub inline fn trsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_strsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (T == f64) {
                    return ci.cblas_dtrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, n), A, scast(c_int, lda), x, scast(c_int, incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/trsv.zig").trsv(T, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn strsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return trsv(f32, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn dtrsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return trsv(f64, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ctrsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const cf32, lda: isize, x: [*]cf32, incx: isize) void {
    return trsv(cf32, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ztrsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const cf64, lda: isize, x: [*]cf64, incx: isize) void {
    return trsv(cf64, order, uplo, transA, diag, n, A, lda, x, incx);
}

// Level 3 BLAS
pub inline fn gemm(comptime T: type, order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_sgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), scast(c_int, m), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                } else if (T == f64) {
                    return ci.cblas_dgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), scast(c_int, m), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), scast(c_int, m), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), scast(c_int, m), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/gemm.zig").gemm(T, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn sgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]const f32, ldb: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return gemm(f32, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn dgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]const f64, ldb: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return gemm(f64, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn cgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: cf32, A: [*]const cf32, lda: isize, B: [*]const cf32, ldb: isize, beta: cf32, C: [*]cf32, ldc: isize) void {
    return gemm(cf32, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: cf64, A: [*]const cf64, lda: isize, B: [*]const cf64, ldb: isize, beta: cf64, C: [*]cf64, ldc: isize) void {
    return gemm(cf64, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn hemm(comptime T: type, order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chemm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhemm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/hemm.zig").hemm(T, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn chemm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: cf32, A: [*]const cf32, lda: isize, B: [*]const cf32, ldb: isize, beta: cf32, C: [*]cf32, ldc: isize) void {
    return hemm(cf32, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zhemm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: cf64, A: [*]const cf64, lda: isize, B: [*]const cf64, ldb: isize, beta: cf64, C: [*]cf64, ldc: isize) void {
    return hemm(cf64, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn herk(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Scalar(T), A: [*]const T, lda: isize, beta: Scalar(T), C: [*]T, ldc: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cherk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), beta, C, scast(c_int, ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zherk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), beta, C, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/herk.zig").herk(T, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn cherk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f32, A: [*]const cf32, lda: isize, beta: f32, C: [*]cf32, ldc: isize) void {
    return herk(cf32, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn zherk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f64, A: [*]const cf64, lda: isize, beta: f64, C: [*]cf64, ldc: isize) void {
    return herk(cf64, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

pub inline fn her2k(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: Scalar(T), C: [*]T, ldc: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cher2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zher2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/her2k.zig").her2k(T, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn cher2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: cf32, A: [*]const cf32, lda: isize, B: [*]const cf32, ldb: isize, beta: f32, C: [*]cf32, ldc: isize) void {
    return her2k(cf32, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zher2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: cf64, A: [*]const cf64, lda: isize, B: [*]const cf64, ldb: isize, beta: f64, C: [*]cf64, ldc: isize) void {
    return her2k(cf64, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn symm(comptime T: type, order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_ssymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                } else if (T == f64) {
                    return ci.cblas_dsymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_csymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zsymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/symm.zig").symm(T, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn ssymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]const f32, ldb: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return symm(f32, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn dsymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]const f64, ldb: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return symm(f64, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn csymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: cf32, A: [*]const cf32, lda: isize, B: [*]const cf32, ldb: isize, beta: cf32, C: [*]cf32, ldc: isize) void {
    return symm(cf32, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zsymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: cf64, A: [*]const cf64, lda: isize, B: [*]const cf64, ldb: isize, beta: cf64, C: [*]cf64, ldc: isize) void {
    return symm(cf64, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn syrk(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_ssyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), beta, C, scast(c_int, ldc));
                } else if (T == f64) {
                    return ci.cblas_dsyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), beta, C, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_csyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), &beta, C, scast(c_int, ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zsyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), &beta, C, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/syrk.zig").syrk(T, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn ssyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return syrk(f32, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn dsyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return syrk(f64, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn csyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: cf32, A: [*]const cf32, lda: isize, beta: cf32, C: [*]cf32, ldc: isize) void {
    return syrk(cf32, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn zsyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: cf64, A: [*]const cf64, lda: isize, beta: cf64, C: [*]cf64, ldc: isize) void {
    return syrk(cf64, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

pub inline fn syr2k(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_ssyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                } else if (T == f64) {
                    return ci.cblas_dsyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), alpha, A, scast(c_int, lda), B, scast(c_int, ldb), beta, C, scast(c_int, ldc));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_csyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zsyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), scast(c_int, n), scast(c_int, k), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb), &beta, C, scast(c_int, ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr2k.zig").syr2k(T, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn ssyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]const f32, ldb: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return syr2k(f32, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn dsyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]const f64, ldb: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return syr2k(f64, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn csyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: cf32, A: [*]const cf32, lda: isize, B: [*]const cf32, ldb: isize, beta: cf32, C: [*]cf32, ldc: isize) void {
    return syr2k(cf32, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zsyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: cf64, A: [*]const cf64, lda: isize, B: [*]const cf64, ldb: isize, beta: cf64, C: [*]cf64, ldc: isize) void {
    return syr2k(cf64, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn trmm(comptime T: type, order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]T, ldb: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_strmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                } else if (T == f64) {
                    return ci.cblas_dtrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                }
            },
            else => {},
        }
    }

    return @import("blas/trmm.zig").trmm(T, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn strmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]f32, ldb: isize) void {
    return trmm(f32, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn dtrmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]f64, ldb: isize) void {
    return trmm(f64, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ctrmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: cf32, A: [*]const cf32, lda: isize, B: [*]cf32, ldb: isize) void {
    return trmm(cf32, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ztrmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: cf64, A: [*]const cf64, lda: isize, B: [*]cf64, ldb: isize) void {
    return trmm(cf64, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

pub inline fn trsm(comptime T: type, order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]T, ldb: isize) void {
    const supported = types.numericType(T);

    if (opts.link_cblas != null) {
        switch (supported) {
            .float => {
                if (T == f32) {
                    return ci.cblas_strsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                } else if (T == f64) {
                    return ci.cblas_dtrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                }
            },
            .cfloat => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), scast(c_int, m), scast(c_int, n), &alpha, A, scast(c_int, lda), B, scast(c_int, ldb));
                }
            },
            else => {},
        }
    }

    return @import("blas/trsm.zig").trsm(T, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn strsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]f32, ldb: isize) void {
    return trsm(f32, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn dtrsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]f64, ldb: isize) void {
    return trsm(f64, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ctrsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: cf32, A: [*]const cf32, lda: isize, B: [*]cf32, ldb: isize) void {
    return trsm(cf32, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ztrsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: cf64, A: [*]const cf64, lda: isize, B: [*]cf64, ldb: isize) void {
    return trsm(cf64, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

pub const Error = error{
    InvalidArgument,
};
