const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");
const blas = @import("../blas.zig");

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
/// Signature
/// ---------
/// ```zig
/// fn rotg(a: *A, b: *B, c: *C, s: *S, ctx: anytype) !void
/// ```
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
pub fn rotg(
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
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and A == C and A == S and options.link_cblas != null) {
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

    return _rotg(a, b, c, s, ctx);
}

fn _rotg(
    a: anytype,
    b: anytype,
    c: anytype,
    s: anytype,
    ctx: anytype,
) !void {
    const A: type = types.Child(@TypeOf(b));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Child(@TypeOf(c));
    const S: type = types.Child(@TypeOf(s));
    const Ca: type = types.Coerce(A, types.Coerce(B, types.Coerce(C, S)));

    if (comptime types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B) or
        types.isArbitraryPrecision(C) or
        types.isArbitraryPrecision(S))
    {
        @compileError("rotg does not support arbitrary precision types.");
    } else {
        const safmin = std.math.floatMin(Scalar(Ca));
        const safmax = std.math.floatMax(Scalar(Ca));

        switch (comptime types.numericType(Ca)) {
            .int, .float => {
                const anorm: A = ops.abs(a.*, ctx) catch unreachable;
                const bnorm: B = ops.abs(b.*, ctx) catch unreachable;

                if (bnorm == 0) {
                    c.* = 1;
                    s.* = 0;
                    b.* = 0;
                } else if (anorm == 0) {
                    c.* = 0;
                    s.* = 1;
                    ops.set(a, b.*, ctx) catch unreachable;
                    b.* = 1;
                } else {
                    const scl = float.min(safmax, float.max(safmin, float.max(anorm, bnorm)));

                    const sigma = if (anorm > bnorm) std.math.sign(a.*) else std.math.sign(b.*); // switch for zml's float.sign when implemented

                    const r = sigma * scl * float.sqrt(((a.* / scl) * (a.* / scl)) + ((b.* / scl) * (b.* / scl)));
                    c.* = a.* / r;
                    s.* = b.* / r;

                    if (anorm > bnorm) {
                        b.* = s.*;
                    } else if (c.* != 0) {
                        b.* = 1 / c.*;
                    } else {
                        b.* = 1;
                    }

                    a.* = r;
                }
            },
            .cfloat => {
                const rtmin = float.sqrt(safmin);
                const rtmax = float.sqrt(safmax / 4);

                const f = a.*;
                var g = b.*;

                if (g.re == 0 and g.im == 0) {
                    c.* = 1.0;
                    s.* = Ca.init(0.0, 0.0);
                    a.* = f;
                } else if (f.re == 0 and f.im == 0) {
                    c.* = 0.0;

                    const g_abs = @sqrt(g.re * g.re + g.im * g.im);
                    const g_conj = g.conj();
                    s.* = Ca.init(g_conj.re / g_abs, g_conj.im / g_abs);
                    a.* = Ca.init(g_abs, 0.0);
                } else {
                    const f1 = @max(@abs(f.re), @abs(f.im));
                    const g1 = @max(@abs(g.re), @abs(g.im));

                    if (f1 > rtmin and f1 < rtmax and g1 > rtmin and g1 < rtmax) {
                        const f2 = f.re * f.re + f.im * f.im;
                        const g2 = g.re * g.re + g.im * g.im;
                        const h2 = f2 + g2;

                        if (f2 >= h2 * safmin) {
                            c.* = @sqrt(f2 / h2);
                            const h = @sqrt(f2 * h2);
                            const g_conj = g.conj();
                            s.* = Ca.init((g_conj.re * f.re + g_conj.im * f.im) / h, (g_conj.im * f.re - g_conj.re * f.im) / h);
                            a.* = Ca.init(f.re / c.*, f.im / c.*);
                        } else {
                            const h = @sqrt(f2 * h2);
                            c.* = f2 / h;
                            const g_conj = g.conj();
                            s.* = Ca.init((g_conj.re * f.re + g_conj.im * f.im) / h, (g_conj.im * f.re - g_conj.re * f.im) / h);
                            a.* = Ca.init(h, 0.0);
                        }
                    } else {
                        const u = @min(safmax, @max(safmin, @max(f1, g1)));
                        const fs = Ca.init(f.re / u, f.im / u);
                        const gs = Ca.init(g.re / u, g.im / u);
                        const f2 = fs.re * fs.re + fs.im * fs.im;
                        const g2 = gs.re * gs.re + gs.im * gs.im;
                        const h2 = f2 + g2;

                        c.* = @sqrt(f2 / h2);
                        const g_conj = gs.conj();
                        s.* = Ca.init((g_conj.re * fs.re + g_conj.im * fs.im) / @sqrt(f2 * h2), (g_conj.im * fs.re - g_conj.re * fs.im) / @sqrt(f2 * h2));
                        a.* = Ca.init(fs.re / c.*, fs.im / c.*);
                        a.* = Ca.init(a.*.re * u, a.*.im * u);
                    }
                }
            },
            else => unreachable,
        }
    }
}
