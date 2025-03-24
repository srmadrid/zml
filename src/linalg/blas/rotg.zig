const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");

const Numeric = core.types.Numeric;

pub inline fn rotg(comptime T: type, a: *T, b: *T, c: *Numeric(T), s: *T) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    const safmin = std.math.floatMin(Numeric(T));
    const safmax = std.math.floatMax(Numeric(T));

    switch (numericType) {
        .bool => @compileError("blas.rotg does not support bool."),
        .int, .float => {
            const anorm: T = @abs(a.*);
            const bnorm: T = @abs(b.*);

            if (bnorm == 0) {
                c.* = 1;
                s.* = 0;
                b.* = 0;
            } else if (anorm == 0) {
                c.* = 0;
                s.* = 1;
                a.* = b.*;
                b.* = 1;
            } else {
                const scl = @min(safmax, @max(safmin, @max(anorm, bnorm)));

                const sigma = if (anorm > bnorm) std.math.sign(a.*) else std.math.sign(b.*);

                const r = sigma * scl * @sqrt(((a.* / scl) * (a.* / scl)) + ((b.* / scl) * (b.* / scl)));
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
            const rtmin = @sqrt(safmin);
            const rtmax = @sqrt(safmax / 4);

            const f = a.*;
            var g = b.*;

            if (g.re == 0 and g.im == 0) {
                c.* = 1.0;
                s.* = T.init(0.0, 0.0);
                a.* = f;
            } else if (f.re == 0 and f.im == 0) {
                c.* = 0.0;

                const g_abs = @sqrt(g.re * g.re + g.im * g.im);
                const g_conj = g.conjugate();
                s.* = T.init(g_conj.re / g_abs, g_conj.im / g_abs);
                a.* = T.init(g_abs, 0.0);
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
                        const g_conj = g.conjugate();
                        s.* = T.init((g_conj.re * f.re + g_conj.im * f.im) / h, (g_conj.im * f.re - g_conj.re * f.im) / h);
                        a.* = T.init(f.re / c.*, f.im / c.*);
                    } else {
                        const h = @sqrt(f2 * h2);
                        c.* = f2 / h;
                        const g_conj = g.conjugate();
                        s.* = T.init((g_conj.re * f.re + g_conj.im * f.im) / h, (g_conj.im * f.re - g_conj.re * f.im) / h);
                        a.* = T.init(h, 0.0);
                    }
                } else {
                    const u = @min(safmax, @max(safmin, @max(f1, g1)));
                    const fs = T.init(f.re / u, f.im / u);
                    const gs = T.init(g.re / u, g.im / u);
                    const f2 = fs.re * fs.re + fs.im * fs.im;
                    const g2 = gs.re * gs.re + gs.im * gs.im;
                    const h2 = f2 + g2;

                    c.* = @sqrt(f2 / h2);
                    const g_conj = gs.conjugate();
                    s.* = T.init((g_conj.re * fs.re + g_conj.im * fs.im) / @sqrt(f2 * h2), (g_conj.im * fs.re - g_conj.re * fs.im) / @sqrt(f2 * h2));
                    a.* = T.init(fs.re / c.*, fs.im / c.*);
                    a.* = T.init(a.*.re * u, a.*.im * u);
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.rotg only supports simple types."),
        .unsupported => unreachable,
    }
}

test "rotg" {
    var a: f64 = 2;
    var b: f64 = 2;
    var c: f64 = undefined;
    var s: f64 = undefined;

    blas.rotg(f64, &a, &b, &c, &s);

    try std.testing.expectApproxEqAbs(2.8284271247461903, a, 0.0000001);
    try std.testing.expectApproxEqAbs(1.4142135623730951, b, 0.0000001);
    try std.testing.expectApproxEqAbs(0.7071067811865475, c, 0.0000001);
    try std.testing.expectApproxEqAbs(0.7071067811865475, s, 0.0000001);

    const Complex = std.math.Complex;
    var a_c: Complex(f64) = Complex(f64).init(1, 2);
    var b_c: Complex(f64) = Complex(f64).init(3, 4);
    var c_c: f64 = undefined;
    var s_c: Complex(f64) = undefined;

    blas.rotg(Complex(f64), &a_c, &b_c, &c_c, &s_c);

    //try std.testing.expectApproxEqAbs(2.449489742783178, a_c.re, 0.0000001);
    //try std.testing.expectApproxEqAbs(-5.375387381226731e102, a_c.im, 0.0000001);
    //try std.testing.expectApproxEqAbs(3, b_c.re, 0.0000001);
    //try std.testing.expectApproxEqAbs(2, b_c.im, 0.0000001);
    //try std.testing.expectApproxEqAbs(0.6172133998483676, c_c, 0.0000001);
    //try std.testing.expectApproxEqAbs(0.7715167498104596, s_c.re, 0.0000001);
    //try std.testing.expectApproxEqAbs(0.15430334996209188, s_c.im, 0.0000001);
}
