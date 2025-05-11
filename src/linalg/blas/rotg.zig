const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

const Scalar = types.Scalar;

pub inline fn rotg(comptime T: type, a: *T, b: *T, c: *Scalar(T), s: *T) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    const safmin = std.math.floatMin(Scalar(T));
    const safmax = std.math.floatMax(Scalar(T));

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
