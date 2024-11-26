const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

pub inline fn dotu(comptime T: type, x: NDArray(T), y: NDArray(T)) !T {
    const supported = core.supported.whatSupportedNumericType(T);

    if (x.size != y.size) {
        return Error.IncompatibleSize;
    }

    var res: T = undefined;
    switch (supported) {
        .BuiltinBool => @compileError("BLAS.dotu does not support bool."),
        .BuiltinInt => @compileError("BLAS.dotu does not support integers. Use BLAS.dot instead."),
        .BuiltinFloat => @compileError("BLAS.dotu does not support floats. Use BLAS.dot instead."),
        .Complex => {
            res = T.init(0, 0);
            for (0..x.size) |i| {
                res.re += x.data[i].re * y.data[i].re - x.data[i].im * y.data[i].im;
                res.im += x.data[i].re * y.data[i].im + x.data[i].im * y.data[i].re;
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.dotu only supports simple types."),
        .Unsupported => unreachable,
    }

    return res;
}

test "dotu" {
    const a = std.testing.allocator;

    const Complex = std.math.Complex;
    var A: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(Complex(f64).init(1, -1));

    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(2, 2));

    const result1 = try NDArray(Complex(f64)).BLAS.dotu(A, B);
    try std.testing.expect(result1.re == 4032 and result1.im == 0);
}
