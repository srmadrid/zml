const std = @import("std");
const core = @import("../core/core.zig");
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn rot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: scalar(T), s: scalar(T)) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0 or (c == 1 and s == 0)) return;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.rot does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
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
                    var temp = c * x[@intCast(ix)] + s * y[@intCast(iy)];
                    y[@intCast(iy)] = c * y[@intCast(iy)] - s * x[@intCast(ix)];
                    x[@intCast(ix)] = temp;

                    temp = c * x[@intCast(ix + incx)] + s * y[@intCast(iy + incy)];
                    y[@intCast(iy + incy)] = c * y[@intCast(iy + incy)] - s * x[@intCast(ix + incx)];
                    x[@intCast(ix + incx)] = temp;

                    temp = c * x[@intCast(ix + incx2)] + s * y[@intCast(iy + incy2)];
                    y[@intCast(iy + incy2)] = c * y[@intCast(iy + incy2)] - s * x[@intCast(ix + incx2)];
                    x[@intCast(ix + incx2)] = temp;

                    temp = c * x[@intCast(ix + incx3)] + s * y[@intCast(iy + incy3)];
                    y[@intCast(iy + incy3)] = c * y[@intCast(iy + incy3)] - s * x[@intCast(ix + incx3)];
                    x[@intCast(ix + incx3)] = temp;

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const temp = c * x[@intCast(ix)] + s * y[@intCast(iy)];
                y[@intCast(iy)] = c * y[@intCast(iy)] - s * x[@intCast(ix)];
                x[@intCast(ix)] = temp;

                ix += incx;
                iy += incy;
            }
        },
        .Complex => {
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 2) << 2;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    var temp = c * x[@intCast(ix)].re + s * y[@intCast(iy)].re;
                    y[@intCast(iy)].re = c * y[@intCast(iy)].re - s * x[@intCast(ix)].re;
                    x[@intCast(ix)].re = temp;

                    temp = c * x[@intCast(ix)].im + s * y[@intCast(iy)].im;
                    y[@intCast(iy)].im = c * y[@intCast(iy)].im - s * x[@intCast(ix)].im;
                    x[@intCast(ix)].im = temp;

                    temp = c * x[@intCast(ix + incx)].re + s * y[@intCast(iy + incy)].re;
                    y[@intCast(iy + incy)].re = c * y[@intCast(iy + incy)].re - s * x[@intCast(ix + incx)].re;
                    x[@intCast(ix + incx)].re = temp;

                    temp = c * x[@intCast(ix + incx)].im + s * y[@intCast(iy + incy)].im;
                    y[@intCast(iy + incy)].im = c * y[@intCast(iy + incy)].im - s * x[@intCast(ix + incx)].im;
                    x[@intCast(ix + incx)].im = temp;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const temp = c * x[@intCast(ix)].re + s * y[@intCast(iy)].re;
                y[@intCast(iy)].re = c * y[@intCast(iy)].re - s * x[@intCast(ix)].re;
                x[@intCast(ix)].re = temp;

                const temp2 = c * x[@intCast(ix)].im + s * y[@intCast(iy)].im;
                y[@intCast(iy)].im = c * y[@intCast(iy)].im - s * x[@intCast(ix)].im;
                x[@intCast(ix)].im = temp2;

                ix += incx;
                iy += incy;
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.rot only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "rot" {
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

    BLAS.rot(f64, n, x1.ptr, 1, x2.ptr, 1, 2, 3);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1))), x1[i]);
        try std.testing.expectEqual(@as(f64, @floatFromInt(@as(isize, @intCast(2 * (i + 1))) - @as(isize, @intCast(3 * (i + 1))))), x2[i]);
        x1[i] = @floatFromInt(i + 1);
    }
    BLAS.rot(f64, n, x1.ptr, 1, x3.ptr, -1, 2, 3);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1))), x1[i]);
        try std.testing.expectEqual(@as(f64, @floatFromInt(@as(isize, @intCast(2 * (n - i))) - @as(isize, @intCast(3 * (n - i))))), x3[i]);
        x1[i] = @floatFromInt(i + 1);
    }
    BLAS.rot(f64, n / 2, x1.ptr, 2, x4.ptr, 2, 2, 3);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1))), x1[i]);
            try std.testing.expectEqual(@as(f64, @floatFromInt(@as(isize, @intCast(2 * (i + 1))) - @as(isize, @intCast(3 * (i + 1))))), x4[i]);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x1[i]);
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
