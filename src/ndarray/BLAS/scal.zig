const std = @import("std");
const core = @import("../../core/core.zig");
const BLAS = @import("BLAS.zig");
const options = @import("options");

pub inline fn scal(comptime T: type, n: isize, a: T, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.scal does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (a == 1) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            const nu = (n >> 3) << 3;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incx3 = incx * 3;
                const incx4 = incx * 4;
                const incx5 = incx * 5;
                const incx6 = incx * 6;
                const incx7 = incx * 7;
                const incx8 = incx * 8;
                while (ix != StX) {
                    x[@intCast(ix)] *= a;
                    x[@intCast(ix + incx)] *= a;
                    x[@intCast(ix + incx2)] *= a;
                    x[@intCast(ix + incx3)] *= a;
                    x[@intCast(ix + incx4)] *= a;
                    x[@intCast(ix + incx5)] *= a;
                    x[@intCast(ix + incx6)] *= a;
                    x[@intCast(ix + incx7)] *= a;

                    ix += incx8;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                x[@intCast(ix)] *= a;

                ix += incx;
            }
        },
        .Complex => {
            if (a.re == 1 and a.im == 0) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            const nu = (n >> 2) << 2;
            if (a.re == 0 and a.im == 0) {
                if (nu != 0) {
                    const StX = ix + nu * incx;
                    const incx2 = incx * 2;
                    const incx3 = incx * 3;
                    const incx4 = incx * 4;
                    while (ix != StX) {
                        x[@intCast(ix)].re = 0;
                        x[@intCast(ix)].im = 0;
                        x[@intCast(ix + incx)].re = 0;
                        x[@intCast(ix + incx)].im = 0;
                        x[@intCast(ix + incx2)].re = 0;
                        x[@intCast(ix + incx2)].im = 0;
                        x[@intCast(ix + incx3)].re = 0;
                        x[@intCast(ix + incx3)].im = 0;

                        ix += incx4;
                    }

                    for (@intCast(nu)..@intCast(n)) |_| {
                        x[@intCast(ix)].re = 0;
                        x[@intCast(ix)].im = 0;

                        ix += incx;
                    }
                }
            } else {
                if (nu != 0) {
                    const StX = ix + nu * incx;
                    const incx2 = incx * 2;
                    const incx3 = incx * 3;
                    const incx4 = incx * 4;
                    while (ix != StX) {
                        const x0 = x[@intCast(ix)];
                        const x1 = x[@intCast(ix + incx)];
                        const x2 = x[@intCast(ix + incx2)];
                        const x3 = x[@intCast(ix + incx3)];

                        x[@intCast(ix)].re = x0.re * a.re - x0.im * a.im;
                        x[@intCast(ix)].im = x0.re * a.im + x0.im * a.re;
                        x[@intCast(ix + incx)].re = x1.re * a.re - x1.im * a.im;
                        x[@intCast(ix + incx)].im = x1.re * a.im + x1.im * a.re;
                        x[@intCast(ix + incx2)].re = x2.re * a.re - x2.im * a.im;
                        x[@intCast(ix + incx2)].im = x2.re * a.im + x2.im * a.re;
                        x[@intCast(ix + incx3)].re = x3.re * a.re - x3.im * a.im;
                        x[@intCast(ix + incx3)].im = x3.re * a.im + x3.im * a.re;

                        ix += incx4;
                    }

                    for (@intCast(nu)..@intCast(n)) |_| {
                        x[@intCast(ix)].re *= a.re;
                        x[@intCast(ix)].im *= a.im;

                        ix += incx;
                    }
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.scal only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "scal" {
    const a: std.mem.Allocator = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);
    var x2 = try a.alloc(f64, n);
    defer a.free(x2);
    var x3 = try a.alloc(f64, n);
    defer a.free(x3);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
        x2[i] = @floatFromInt(n - i);
        x3[i] = @floatFromInt(i + 1);
    }

    BLAS.scal(f64, n, 2, x1.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1))), x1[i]);
    }
    if (options.use_cblas == null) {
        BLAS.scal(f64, n, 2, x2.ptr, -1);
        for (0..n) |i| {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i))), x2[i]);
        }
    }
    BLAS.scal(f64, n / 2, 2, x3.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1))), x3[i]);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x3[i]);
        }
    }

    var x4 = try a.alloc(Complex(f64), n);
    defer a.free(x4);
    var x5 = try a.alloc(Complex(f64), n);
    defer a.free(x5);
    var x6 = try a.alloc(Complex(f64), n);
    defer a.free(x6);

    for (0..n) |i| {
        x4[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
        x5[i] = Complex(f64).init(@floatFromInt(n - i), @floatFromInt(-@as(isize, @intCast((n - i)))));
        x6[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
    }

    BLAS.scal(Complex(f64), n, Complex(f64).init(2, 3), x4.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1))), x4[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)))), x4[i].im);
    }
    if (options.use_cblas == null) {
        BLAS.scal(Complex(f64), n, Complex(f64).init(2, 3), x5.ptr, -1);
        for (0..n) |i| {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i) + 3 * (n - i))), x5[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(n - i)) + 3 * @as(isize, @intCast(n - i)))), x5[i].im);
        }
    }
    BLAS.scal(Complex(f64), n / 2, Complex(f64).init(2, 3), x6.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1))), x6[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)))), x6[i].im);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x6[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x6[i].im);
        }
    }
}