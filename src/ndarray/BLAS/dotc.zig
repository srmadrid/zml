const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn dotc(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    var temp: T = T.init(0, 0);

    if (n <= 0) return temp;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.dotc does not support bool."),
        .BuiltinInt => @compileError("BLAS.dotc does not support integers. Use BLAS.dot instead."),
        .BuiltinFloat => @compileError("BLAS.dotc does not support floats. Use BLAS.dot instead."),
        .Complex => {
            if (incx == 1 and incy == 1) {
                for (0..@intCast(n)) |i| {
                    temp.re += x[i].re * y[i].re + x[i].im * y[i].im;
                    temp.im += x[i].re * y[i].im - x[i].im * y[i].re;
                }
            } else {
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..@intCast(n)) |_| {
                    temp.re += x[@intCast(ix)].re * y[@intCast(iy)].re + x[@intCast(ix)].im * y[@intCast(iy)].im;
                    temp.im += x[@intCast(ix)].re * y[@intCast(iy)].im - x[@intCast(ix)].im * y[@intCast(iy)].re;
                    ix += incx;
                    iy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.dotc only supports simple types."),
        .Unsupported => unreachable,
    }

    return temp;
}

test "dotc" {
    const a = std.testing.allocator;

    const Complex = std.math.Complex;
    var A: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(Complex(f64).init(1, -1));

    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(2, 2));

    const result1 = try NDArray(Complex(f64)).BLAS.dotc(A.flatten(), B.flatten());
    try std.testing.expect(result1.re == 0 and result1.im == 4032);
}
