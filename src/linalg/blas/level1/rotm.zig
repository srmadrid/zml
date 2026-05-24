const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const blas = @import("../blas.zig");

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
/// `n` (`i32`): Specifies the number of elements in vectors `x` and `y`. Must
/// be greater than 0.
///
/// `x` (mutable many-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incx)`. On return, every element `x[i]` is replaced
/// by `h11 * x[i] + h12 * y[i]`.
///
/// `incx` (`i32`): Specifies the increment for the elements of `x`.
///
/// `y` (mutable many-item pointer to `bool`, `int`, `float`, `integer`,
/// `rational`, `real` or `expression`): Array, size at least
/// `1 + (n - 1) * abs(incy)`. On return, every element `y[i]` is replaced
/// by `h21 * x[i] + h22 * y[i]`.
///
/// `incy` (`i32`): Specifies the increment for the elements of `y`.
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
/// Signature
/// ---------
/// ```zig
/// fn rotm(n: i32, x: [*]X, incx: i32, y: [*]Y, incy: i32, param: [*]const P, ctx: anytype) !void
/// ```
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
pub fn rotm(
    n: i32,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
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
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime X == Y and X == P and options.link_cblas != null) {
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

    return _rotm(n, x, incx, y, incy, param, ctx);
}

fn _rotm(
    n: i32,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    param: anytype,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const P: type = types.Child(@TypeOf(param));
    //const C: type = types.Coerce(X, types.Coerce(Y, P));

    if (n == 0 or param[0] == -2) return blas.Error.InvalidArgument;

    if (comptime types.isArbitraryPrecision(X) or
        types.isArbitraryPrecision(Y) or
        types.isArbitraryPrecision(P))
    {
        @compileError("zml.linalg.blas.rotm not implemented for arbitrary precision types yet");
    } else {
        const flag = param[0];

        if (flag == -1) {
            const h11: P = param[1];
            const h21: P = param[2];
            const h12: P = param[3];
            const h22: P = param[4];
            var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
            var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
            for (0..scast(u32, n)) |_| {
                const x0: X = x[scast(u32, ix)];
                ops.add_( // x[ix] = h11 * x[ix] + h12 * y[iy]
                    &x[scast(u32, ix)],
                    ops.mul(
                        h11,
                        x[scast(u32, ix)],
                        ctx,
                    ) catch unreachable,
                    ops.mul(
                        h12,
                        y[scast(u32, iy)],
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ops.add_( // y[iy] = h21 * x[ix] + h22 * y[iy]
                    &y[scast(u32, iy)],
                    ops.mul(
                        h21,
                        x0,
                        ctx,
                    ) catch unreachable,
                    ops.mul(
                        h22,
                        y[scast(u32, iy)],
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        } else if (flag == 0) {
            const h21: P = param[2];
            const h12: P = param[3];
            var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
            var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
            for (0..scast(u32, n)) |_| {
                const x0 = x[scast(u32, ix)];
                ops.add_( // x[ix] = h12 * y[iy] + x[ix]
                    &x[scast(u32, ix)],
                    x[scast(u32, ix)],
                    ops.mul(
                        h12,
                        y[scast(u32, iy)],
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ops.add_( // y[iy] = h21 * x[ix] + y[iy]
                    &y[scast(u32, iy)],
                    y[scast(u32, iy)],
                    ops.mul(
                        h21,
                        x0,
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        } else if (flag == 1) {
            const h11 = param[1];
            const h22 = param[4];
            var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
            var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
            for (0..scast(u32, n)) |_| {
                const x0 = x[scast(u32, ix)];
                ops.add_( // x[ix] = h11 * x[ix] - h22 * y[iy]
                    &x[scast(u32, ix)],
                    ops.mul(
                        h11,
                        x[scast(u32, ix)],
                        ctx,
                    ) catch unreachable,
                    y[scast(u32, iy)],
                    ctx,
                ) catch unreachable;

                ops.sub_( // y[iy] = h22 * y[iy] - h11 * x[ix]
                    &y[scast(u32, iy)],
                    ops.mul(
                        h22,
                        y[scast(u32, iy)],
                        ctx,
                    ) catch unreachable,
                    x0,
                    ctx,
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        }
    }
}
