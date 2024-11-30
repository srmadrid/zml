const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

const scalar = core.supported.scalar;

pub inline fn rot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: scalar(T), s: scalar(T)) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.nrm2 does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (incx == 1 and incy == 1) {
                for (0..@intCast(n)) |i| {
                    const temp = c * x[i] + s * y[i];
                    y[i] = c * y[i] - s * x[i];
                    x[i] = temp;
                }
            } else {
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..@intCast(n)) |_| {
                    const temp = c * x[@intCast(ix)] + s * y[@intCast(iy)];
                    y[@intCast(iy)] = c * y[@intCast(iy)] - s * x[@intCast(ix)];
                    x[@intCast(ix)] = temp;
                    ix += incx;
                    iy += incy;
                }
            }
        },
        .Complex => {
            if (incx == 1 and incy == 1) {
                for (0..@intCast(n)) |i| {
                    var temp = c * x[i].re + s * y[i].re;
                    y[i].re = c * y[i].re - s * x[i].re;
                    x[i].re = temp;

                    temp = c * x[i].im + s * y[i].im;
                    y[i].im = c * y[i].im - s * x[i].im;
                    x[i].im = temp;
                }
            } else {
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..@intCast(n)) |_| {
                    var temp = c * x[@intCast(ix)].re + s * y[@intCast(iy)].re;
                    y[@intCast(iy)].re = c * y[@intCast(iy)].re - s * x[@intCast(ix)].re;
                    x[@intCast(ix)].re = temp;

                    temp = c * x[@intCast(ix)].im + s * y[@intCast(iy)].im;
                    y[@intCast(iy)].im = c * y[@intCast(iy)].im - s * x[@intCast(ix)].im;
                    x[@intCast(ix)].im = temp;

                    ix += incx;
                    iy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.nrm2 only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "rot" {
    const a = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(1);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(1);

    try NDArray(f64).BLAS.rot(@constCast(&A.flatten()), @constCast(&B.flatten()), 0.7071067811865475, 0.7071067811865475);

    for (0..A.size) |i| {
        try std.testing.expectApproxEqAbs(0.7071067811865475 * 2, A.data[i], 0.0001);
        try std.testing.expectApproxEqAbs(0, B.data[i], 0.0001);
    }

    const Complex = std.math.Complex;
    var C: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer C.deinit();

    C.setAll(Complex(f64).init(1, -1));

    var D: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer D.deinit();

    D.setAll(Complex(f64).init(1, -1));

    try NDArray(Complex(f64)).BLAS.rot(@constCast(&C.flatten()), @constCast(&D.flatten()), 0.7071067811865475, 0.7071067811865475);

    for (0..C.size) |i| {
        try std.testing.expectApproxEqAbs(0.7071067811865475 * 2, C.data[i].re, 0.0001);
        try std.testing.expectApproxEqAbs(-0.7071067811865475 * 2, C.data[i].im, 0.0001);
        try std.testing.expectApproxEqAbs(0, D.data[i].re, 0.0001);
        try std.testing.expectApproxEqAbs(0, D.data[i].im, 0.0001);
    }
}
