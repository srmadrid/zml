const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");

pub inline fn dotc(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    var sum: T = T.init(0, 0);

    if (n <= 0) return sum;

    switch (numericType) {
        .bool => @compileError("blas.dotc does not support bool."),
        .int => @compileError("blas.dotc does not support integers. Use blas.dot instead."),
        .float => @compileError("blas.dotc does not support floats. Use blas.dot instead."),
        .cfloat => {
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 1) << 1;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    sum.re += x[@intCast(ix)].re * y[@intCast(iy)].re + x[@intCast(ix)].im * y[@intCast(iy)].im + x[@intCast(ix + incx)].re * y[@intCast(iy + incy)].re + x[@intCast(ix + incx)].im * y[@intCast(iy + incy)].im;
                    sum.im += x[@intCast(ix)].re * y[@intCast(iy)].im - x[@intCast(ix)].im * y[@intCast(iy)].re + x[@intCast(ix + incx)].re * y[@intCast(iy + incy)].im - x[@intCast(ix + incx)].im * y[@intCast(iy + incy)].re;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                sum.re += x[@intCast(ix)].re * y[@intCast(iy)].re + x[@intCast(ix)].im * y[@intCast(iy)].im;
                sum.im += x[@intCast(ix)].re * y[@intCast(iy)].im - x[@intCast(ix)].im * y[@intCast(iy)].re;

                ix += incx;
                iy += incy;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.dotc only supports simple types."),
        .unsupported => unreachable,
    }

    return sum;
}

test "dotc" {
    const a: std.mem.Allocator = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 1000;

    var x1 = try a.alloc(Complex(f64), n);
    defer a.free(x1);
    var x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);
    var x3 = try a.alloc(Complex(f64), n);
    defer a.free(x3);
    var x4 = try a.alloc(Complex(f64), n);
    defer a.free(x4);

    for (0..n) |i| {
        x1[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
        x2[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
        x3[i] = Complex(f64).init(@floatFromInt(n - i), @floatFromInt(-@as(isize, @intCast(n - i))));
        x4[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
    }

    const result1 = blas.dotc(Complex(f64), n, x1.ptr, 1, x2.ptr, 1);
    try std.testing.expectEqual(667667000, result1.re);
    try std.testing.expectEqual(0, result1.im);
    const result2 = blas.dotc(Complex(f64), n, x1.ptr, 1, x3.ptr, -1);
    try std.testing.expectEqual(667667000, result2.re);
    try std.testing.expectEqual(0, result2.im);
    const result3 = blas.dotc(Complex(f64), n / 2, x1.ptr, 2, x4.ptr, 2);
    try std.testing.expectEqual(333333000, result3.re);
    try std.testing.expectEqual(0, result3.im);
}
