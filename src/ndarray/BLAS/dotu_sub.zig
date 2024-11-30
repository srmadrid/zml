const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn dotu_sub(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize, ret: *T) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    ret.* = T.init(0, 0);

    if (n <= 0) return;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.dotu does not support bool."),
        .BuiltinInt => @compileError("BLAS.dotu does not support integers. Use BLAS.dot instead."),
        .BuiltinFloat => @compileError("BLAS.dotu does not support floats. Use BLAS.dot instead."),
        .Complex => {
            if (incx == 1 and incy == 1) {
                for (0..@intCast(n)) |i| {
                    ret.*.re += x[i].re * y[i].re - x[i].im * y[i].im;
                    ret.*.im += x[i].re * y[i].im + x[i].im * y[i].re;
                }
            } else {
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..@intCast(n)) |_| {
                    ret.*.re += x[@intCast(ix)].re * y[@intCast(iy)].re - x[@intCast(ix)].im * y[@intCast(iy)].im;
                    ret.*.im += x[@intCast(ix)].re * y[@intCast(iy)].im + x[@intCast(ix)].im * y[@intCast(iy)].re;
                    ix += incx;
                    iy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.dotu only supports simple types."),
        .Unsupported => unreachable,
    }

    return;
}

test "dotu_sub" {
    const a = std.testing.allocator;

    const Complex = std.math.Complex;
    var A: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(Complex(f64).init(1, -1));

    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(2, 2));

    var result1: Complex(f64) = undefined;
    try NDArray(Complex(f64)).BLAS.dotu_sub(A.flatten(), B.flatten(), &result1);
    try std.testing.expect(result1.re == 4032 and result1.im == 0);
}
