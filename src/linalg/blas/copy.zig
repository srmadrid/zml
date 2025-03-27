const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");

pub inline fn copy(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0) return;

    switch (numericType) {
        .bool, .int, .float, .cfloat => {
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
                    y[@intCast(iy)] = x[@intCast(ix)];
                    y[@intCast(iy + incy)] = x[@intCast(ix + incx)];
                    y[@intCast(iy + incy2)] = x[@intCast(ix + incx2)];
                    y[@intCast(iy + incy3)] = x[@intCast(ix + incx3)];

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                y[@intCast(iy)] = x[@intCast(ix)];

                ix += incx;
                iy += incy;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.copy only supports simple types."),
        .unsupported => unreachable,
    }
}

test copy {
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
        x2[i] = 0;
        x3[i] = 0;
        x4[i] = 0;
    }

    blas.copy(f64, n, x1.ptr, 1, x2.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x2[i]);
    }
    blas.copy(f64, n, x1.ptr, 1, x3.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(n - i)), x3[i]);
    }
    blas.copy(f64, n / 2, x1.ptr, 2, x4.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x4[i]);
        } else {
            try std.testing.expectEqual(0, x4[i]);
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
        x5[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
        x6[i] = Complex(f64).init(0, 0);
        x7[i] = Complex(f64).init(0, 0);
        x8[i] = Complex(f64).init(0, 0);
    }

    blas.copy(Complex(f64), n, x5.ptr, 1, x6.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x6[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x6[i].im);
    }
    blas.copy(Complex(f64), n, x5.ptr, 1, x7.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(n - i)), x7[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(n - i)))), x7[i].im);
    }
    blas.copy(Complex(f64), n / 2, x5.ptr, 2, x8.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x8[i].im);
        } else {
            try std.testing.expectEqual(0, x8[i].re);
            try std.testing.expectEqual(0, x8[i].im);
        }
    }
}
