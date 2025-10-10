const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");
const blas = @import("../blas.zig");

pub fn rotg(
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
