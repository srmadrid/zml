const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn axpy(comptime T: type, n: usize, a: T, x: [*]T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (n == 0) {
        return;
    }

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.axpy does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (a == 0) {
                return;
            }

            if (incx == 1 and incy == 1) {

                // Code for both increments equal to 1.
                const m: usize = n % 4;

                // Clean-up loop for the remainder.
                for (0..m) |i| {
                    y[i] += a * x[i];
                }

                if (n < 4) {
                    return;
                }

                // Process in blocks of 4 for efficiency.
                var i: usize = m;
                while (i < n) : (i += 4) {
                    y[i] += a * x[i];
                    y[i + 1] += a * x[i + 1];
                    y[i + 2] += a * x[i + 2];
                    y[i + 3] += a * x[i + 3];
                }
            } else {
                // Code for unequal increments or increments not equal to 1.
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..n) |_| {
                    y[@intCast(iy)] += a * x[@intCast(ix)];
                    ix += incx;
                    iy += incy;
                }
            }
        },
        .Complex => {
            if (a.re == 0 and a.im == 0) {
                return;
            }

            if (incx == 1 and incy == 1) {
                // Code for both increments equal to 1.
                const m: usize = n % 4;

                // Clean-up loop for the remainder.
                for (0..m) |i| {
                    y[i].re += a.re * x[i].re - a.im * x[i].im;
                    y[i].im += a.re * x[i].im + a.im * x[i].re;
                }

                if (n < 4) {
                    return;
                }

                // Process in blocks of 4 for efficiency.
                var i: usize = m;
                while (i < n) : (i += 4) {
                    y[i].re += a.re * x[i].re - a.im * x[i].im;
                    y[i].im += a.re * x[i].im + a.im * x[i].re;
                    y[i + 1].re += a.re * x[i + 1].re - a.im * x[i + 1].im;
                    y[i + 1].im += a.re * x[i + 1].im + a.im * x[i + 1].re;
                    y[i + 2].re += a.re * x[i + 2].re - a.im * x[i + 2].im;
                    y[i + 2].im += a.re * x[i + 2].im + a.im * x[i + 2].re;
                    y[i + 3].re += a.re * x[i + 3].re - a.im * x[i + 3].im;
                    y[i + 3].im += a.re * x[i + 3].im + a.im * x[i + 3].re;
                }
            } else {
                // Code for unequal increments or increments not equal to 1.
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..n) |_| {
                    y[@intCast(iy)].re += a.re * x[@intCast(ix)].re - a.im * x[@intCast(ix)].im;
                    y[@intCast(iy)].im += a.re * x[@intCast(ix)].im + a.im * x[@intCast(ix)].re;
                    ix += incx;
                    iy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.axpy only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "axpy" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer A.deinit();

    A.setAll(1);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer B.deinit();

    try NDArray(f64).BLAS.axpy(2, A.flatten(), @constCast(&B.flatten()));

    for (0..B.size) |i| {
        try std.testing.expect(B.data[i] == 2);
    }

    const Complex = std.math.Complex;
    var C: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer C.deinit();

    C.setAll(Complex(f64).init(1, -1));

    var D: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer D.deinit();

    try NDArray(Complex(f64)).BLAS.axpy(Complex(f64).init(2, 2), C.flatten(), @constCast(&D.flatten()));

    for (0..D.size) |i| {
        try std.testing.expect(std.math.approxEqAbs(f64, D.data[i].re, 4, 0.0001) and std.math.approxEqAbs(f64, D.data[i].im, 0, 0.0001));
    }
}
