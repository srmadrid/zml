const std = @import("std");
const core = @import("../core/core.zig");
const blas = @import("blas.zig");
const options = @import("options");

pub inline fn swap(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    switch (supported) {
        .BuiltinBool => @compileError("blas.swap does not support bool."),
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
                    const x0 = x[@intCast(ix)];
                    const x1 = x[@intCast(ix + incx)];
                    const x2 = x[@intCast(ix + incx2)];
                    const x3 = x[@intCast(ix + incx3)];

                    x[@intCast(ix)] = y[@intCast(iy)];
                    y[@intCast(iy)] = x0;
                    x[@intCast(ix + incx)] = y[@intCast(iy + incy)];
                    y[@intCast(iy + incy)] = x1;
                    x[@intCast(ix + incx2)] = y[@intCast(iy + incy2)];
                    y[@intCast(iy + incy2)] = x2;
                    x[@intCast(ix + incx3)] = y[@intCast(iy + incy3)];
                    y[@intCast(iy + incy3)] = x3;

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const x0 = x[@intCast(ix)];

                x[@intCast(ix)] = y[@intCast(iy)];
                y[@intCast(iy)] = x0;

                ix += incx;
                iy += incy;
            }
        },
        .Complex => {
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 1) << 1;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    const x0 = x[@intCast(ix)];
                    const x1 = x[@intCast(ix + incx)];

                    x[@intCast(ix)] = y[@intCast(iy)];
                    y[@intCast(iy)] = x0;
                    x[@intCast(ix + incx)] = y[@intCast(iy + incy)];
                    y[@intCast(iy + incy)] = x1;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const x0 = x[@intCast(ix)];

                x[@intCast(ix)] = y[@intCast(iy)];
                y[@intCast(iy)] = x0;

                ix += incx;
                iy += incy;
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.swap only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "swap" {
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
        x2[i] = @floatFromInt(i + 2);
        x3[i] = @floatFromInt(n - i + 1);
        x4[i] = @floatFromInt(i + 2);
    }

    blas.swap(f64, n, x1.ptr, 1, x2.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x1[i]);
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x2[i]);
        x1[i] = @floatFromInt(i + 1);
    }
    blas.swap(f64, n, x1.ptr, 1, x3.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x1[i]);
        try std.testing.expectEqual(@as(f64, @floatFromInt(n - i)), x3[i]);
        x1[i] = @floatFromInt(i + 1);
    }
    blas.swap(f64, n / 2, x1.ptr, 2, x4.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x1[i]);
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x4[i]);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x1[i]);
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x4[i]);
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
        x6[i] = Complex(f64).init(@floatFromInt(i + 2), @floatFromInt(-@as(isize, @intCast((i + 2)))));
        x7[i] = Complex(f64).init(@floatFromInt(n - i + 1), @floatFromInt(-@as(isize, @intCast((n - i + 1)))));
        x8[i] = Complex(f64).init(@floatFromInt(i + 2), @floatFromInt(-@as(isize, @intCast((i + 2)))));
    }

    blas.swap(Complex(f64), n, x5.ptr, 1, x6.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x5[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 2)))), x5[i].im);
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x6[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x6[i].im);
        x5[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
    }
    blas.swap(Complex(f64), n, x5.ptr, 1, x7.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x5[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 2)))), x5[i].im);
        try std.testing.expectEqual(@as(f64, @floatFromInt(n - i)), x7[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(n - i)))), x7[i].im);
        x5[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
    }
    blas.swap(Complex(f64), n / 2, x5.ptr, 2, x8.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x5[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 2)))), x5[i].im);
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x8[i].im);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x5[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x5[i].im);
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 2)), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 2)))), x8[i].im);
        }
    }
}
