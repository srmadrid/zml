const std = @import("std");
const core = @import("../../core/core.zig");
const BLAS = @import("BLAS.zig");

pub inline fn axpy(comptime T: type, n: isize, a: T, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.axpy does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (a == 0) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 2) << 2;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incx3 = incx * 3;
                const incx4 = incx * 4;
                const incy2 = incy * 2;
                const incy3 = incy * 3;
                const incy4 = incy * 4;
                while (ix != StX) {
                    y[@intCast(iy)] += a * x[@intCast(ix)];
                    y[@intCast(iy + incy)] += a * x[@intCast(ix + incx)];
                    y[@intCast(iy + incy2)] += a * x[@intCast(ix + incx2)];
                    y[@intCast(iy + incy3)] += a * x[@intCast(ix + incx3)];

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                y[@intCast(iy)] += a * x[@intCast(ix)];

                ix += incx;
                iy += incy;
            }
        },
        .Complex => {
            if (a.re == 0 and a.im == 0) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 1) << 1;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    y[@intCast(iy)].re += a.re * x[@intCast(ix)].re - a.im * x[@intCast(ix)].im;
                    y[@intCast(iy)].im += a.re * x[@intCast(ix)].im + a.im * x[@intCast(ix)].re;
                    y[@intCast(iy + incy)].re += a.re * x[@intCast(ix + incx)].re - a.im * x[@intCast(ix + incx)].im;
                    y[@intCast(iy + incy)].im += a.re * x[@intCast(ix + incx)].im + a.im * x[@intCast(ix + incx)].re;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                y[@intCast(iy)].re += a.re * x[@intCast(ix)].re - a.im * x[@intCast(ix)].im;
                y[@intCast(iy)].im += a.re * x[@intCast(ix)].im + a.im * x[@intCast(ix)].re;

                ix += incx;
                iy += incy;
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.axpy only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "axpy" {
    const a: std.mem.Allocator = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);
    var x2 = try a.alloc(f64, n);
    defer a.free(x2);
    var x3 = try a.alloc(f64, n);
    defer a.free(x3);
    var x4 = try a.alloc(f64, n);
    defer a.free(x4);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
        x2[i] = @floatFromInt(i + 1);
        x3[i] = @floatFromInt(n - i);
        x4[i] = @floatFromInt(i + 1);
    }

    BLAS.axpy(f64, n, 2, x1.ptr, 1, x2.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + (i + 1))), x2[i]);
    }
    BLAS.axpy(f64, n, 2, x1.ptr, 1, x3.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i) + (n - i))), x3[i]);
    }
    BLAS.axpy(f64, n / 2, 2, x1.ptr, 2, x4.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + (i + 1))), x4[i]);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x4[i]);
        }
    }

    var x5 = try a.alloc(Complex(f64), n);
    defer a.free(x5);
    var x6 = try a.alloc(Complex(f64), n);
    defer a.free(x6);
    var x7 = try a.alloc(Complex(f64), n);
    defer a.free(x7);
    var x8 = try a.alloc(Complex(f64), n);
    defer a.free(x8);

    for (0..n) |i| {
        x5[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
        x6[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
        x7[i] = Complex(f64).init(@floatFromInt(n - i), @floatFromInt(-@as(isize, @intCast((n - i)))));
        x8[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
    }

    BLAS.axpy(Complex(f64), n, Complex(f64).init(2, 3), x5.ptr, 1, x6.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1) + (i + 1))), x6[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)) - @as(isize, @intCast(i + 1)))), x6[i].im);
    }
    BLAS.axpy(Complex(f64), n, Complex(f64).init(2, 3), x5.ptr, 1, x7.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i) + 3 * (n - i) + (n - i))), x7[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(n - i)) + 3 * @as(isize, @intCast(n - i)) - @as(isize, @intCast(n - i)))), x7[i].im);
    }
    BLAS.axpy(Complex(f64), n / 2, Complex(f64).init(2, 3), x5.ptr, 2, x8.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1) + (i + 1))), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)) - @as(isize, @intCast(i + 1)))), x8[i].im);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x8[i].im);
        }
    }
}
