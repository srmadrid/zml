const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const blas = @import("../blas.zig");

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
/// Signature
/// ---------
/// ```zig
/// fn rotmg(d1: *D1, d2: *D2, x1: *X1, y1: Y1, param: [*]P, ctx: anytype) !void
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
pub fn rotmg(
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
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime D1 == D2 and D1 == X1 and D1 == P and types.canCoerce(Y1, D1) and options.link_cblas != null) {
        switch (comptime types.numericType(D1)) {
            .float => {
                if (comptime D1 == f32) {
                    return ci.cblas_srotmg(d1, d2, x1, scast(D1, y1), param);
                } else if (comptime D1 == f64) {
                    return ci.cblas_drotmg(d1, d2, x1, scast(D1, y1), param);
                }
            },
            else => {},
        }
    }

    return _rotmg(d1, d2, x1, y1, param, ctx);
}

fn _rotmg(
    d1: anytype,
    d2: anytype,
    x1: anytype,
    y1: anytype,
    param: anytype,
    ctx: anytype,
) !void {
    const D1: type = types.Child(@TypeOf(d1));
    const D2: type = types.Child(@TypeOf(d2));
    const X1: type = types.Child(@TypeOf(x1));
    const Y1: type = @TypeOf(y1);
    const P: type = types.Child(@TypeOf(param));
    const C: type = types.Coerce(D1, types.Coerce(D2, types.Coerce(X1, types.Coerce(Y1, P))));

    if (comptime types.isArbitraryPrecision(D1) or
        types.isArbitraryPrecision(D2) or
        types.isArbitraryPrecision(X1) or
        types.isArbitraryPrecision(Y1) or
        types.isArbitraryPrecision(P))
    {
        @compileError("zml.linalg.blas.rotmg not implemented for arbitrary precision types yet");
    } else {
        const gam: types.EnsureFloat(C) = 4096;
        const gamsq: types.EnsureFloat(C) = 16777216;
        const rgamsq: types.EnsureFloat(C) = 5.9604645e-8;

        var flag: types.EnsureFloat(C) = 0;
        var h11: types.EnsureFloat(C) = 0;
        var h12: types.EnsureFloat(C) = 0;
        var h21: types.EnsureFloat(C) = 0;
        var h22: types.EnsureFloat(C) = 0;
        var p1: types.EnsureFloat(C) = 0;
        var p2: types.EnsureFloat(C) = 0;
        var q1: types.EnsureFloat(C) = 0;
        var q2: types.EnsureFloat(C) = 0;
        var temp: types.EnsureFloat(C) = 0;
        var u: types.EnsureFloat(C) = 0;

        if (d1.* < 0) {
            flag = -1;
            h11 = 0;
            h12 = 0;
            h21 = 0;
            h22 = 0;

            d1.* = 0;
            d2.* = 0;
            x1.* = 0;
        } else {
            ops.mul_( // p2 = d2 * y1
                &p2,
                d2.*,
                y1,
                ctx,
            ) catch unreachable;
            if (p2 == 0) {
                flag = -2;

                ops.set( // param[0] = flag
                    &param[0],
                    flag,
                    ctx,
                ) catch unreachable;

                return;
            }

            ops.mul_( // p1 = d1 * x1
                &p1,
                d1.*,
                x1.*,
                ctx,
            ) catch unreachable;

            ops.mul_( // q2 = p2 * y1
                &q2,
                p2,
                y1,
                ctx,
            ) catch unreachable;

            ops.mul_( // q1 = p1 * x1
                &q1,
                p1,
                x1.*,
                ctx,
            ) catch unreachable;

            if (ops.abs(q1, ctx) catch unreachable > ops.abs(q2, ctx) catch unreachable) {
                ops.div_( // h21 = -y1 / x1
                    &h21,
                    -y1,
                    x1.*,
                    ctx,
                ) catch unreachable;

                ops.div_( // h12 = p1 / p2
                    &h12,
                    p1,
                    p2,
                    ctx,
                ) catch unreachable;

                ops.sub_( // u = 1 - h12 * h21
                    &u,
                    1,
                    ops.mul(
                        h12,
                        h21,
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                if (u > 0) {
                    flag = 0;

                    ops.div_( // d1 /= u
                        d1,
                        d1.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // d2 /= u
                        d2,
                        d2.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.mul_( // x1 *= u
                        x1,
                        x1.*,
                        u,
                        ctx,
                    ) catch unreachable;
                } else {
                    flag = -1;
                    h11 = 0;
                    h12 = 0;
                    h21 = 0;
                    h22 = 0;

                    d1.* = 0;
                    d2.* = 0;
                    x1.* = 0;
                }
            } else {
                if (q2 < 0) {
                    flag = -1;
                    h11 = 0;
                    h12 = 0;
                    h21 = 0;
                    h22 = 0;

                    d1.* = 0;
                    d2.* = 0;
                    x1.* = 0;
                } else {
                    flag = 1;

                    ops.div_( // h11 = p1 / p2
                        &h11,
                        p1,
                        p2,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // h22 = x1 / y1
                        &h22,
                        x1.*,
                        y1,
                        ctx,
                    ) catch unreachable;

                    ops.add_( // u = 1 + h11 * h22
                        &u,
                        1,
                        ops.mul(
                            h11,
                            h22,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // temp = d2 / u
                        &temp,
                        d2.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // d2 = d1 / u
                        d2,
                        d1.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.set( // d1 = temp
                        d1,
                        temp,
                        ctx,
                    ) catch unreachable;

                    ops.mul_( // x1 = y1 * u
                        x1,
                        y1,
                        u,
                        ctx,
                    ) catch unreachable;
                }
            }

            if (d1.* != 0) {
                while ((d1.* <= rgamsq) or (d1.* >= gamsq)) {
                    if (flag == 0) {
                        h11 = 1;
                        h22 = 1;
                        flag = -1;
                    } else {
                        h21 = -1;
                        h12 = 1;
                        flag = -1;
                    }
                    if (d1.* <= rgamsq) {
                        ops.mul_( // d1 *= gam * gam
                            d1,
                            d1.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // x1 /= gam
                            x1,
                            x1.*,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h11 /= gam
                            &h11,
                            h11,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h12 /= gam
                            &h12,
                            h12,
                            gam,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.div_( // d1 /= gam * gam
                            d1,
                            d1.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // x1 *= gam
                            x1,
                            x1.*,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h11 *= gam
                            &h11,
                            h11,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h12 *= gam
                            &h12,
                            h12,
                            gam,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }

            if (d2.* != 0) {
                while ((ops.abs(d2.*, ctx) catch unreachable <= rgamsq) or
                    (ops.abs(d2.*, ctx) catch unreachable >= gamsq))
                {
                    if (flag == 0) {
                        h11 = 1;
                        h22 = 1;
                        flag = -1;
                    } else {
                        h21 = -1;
                        h12 = 1;
                        flag = -1;
                    }
                    if (ops.abs(d2.*, ctx) catch unreachable <= rgamsq) {
                        ops.mul_( // d2 *= gam * gam
                            d2,
                            d2.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h21 /= gam
                            &h21,
                            h21,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h22 /= gam
                            &h22,
                            h22,
                            gam,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.div_( // d2 /= gam * gam
                            d2,
                            d2.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h21 *= gam
                            &h21,
                            h21,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h22 *= gam
                            &h22,
                            h22,
                            gam,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }
        }

        if (flag < 0) {
            ops.set( // param[1] = h11
                &param[1],
                h11,
                ctx,
            ) catch unreachable;

            ops.set( // param[2] = h21
                &param[2],
                h21,
                ctx,
            ) catch unreachable;

            ops.set( // param[3] = h12
                &param[3],
                h12,
                ctx,
            ) catch unreachable;

            ops.set( // param[4] = h22
                &param[4],
                h22,
                ctx,
            ) catch unreachable;
        } else if (flag == 0) {
            ops.set( // param[2] = h21
                &param[2],
                h21,
                ctx,
            ) catch unreachable;

            ops.set( // param[3] = h12
                &param[3],
                h12,
                ctx,
            ) catch unreachable;
        } else {
            ops.set( // param[1] = h11
                &param[1],
                h11,
                ctx,
            ) catch unreachable;

            ops.set( // param[4] = h22
                &param[4],
                h22,
                ctx,
            ) catch unreachable;
        }

        ops.set( // param[0] = flag
            &param[0],
            flag,
            ctx,
        ) catch unreachable;
    }
}
