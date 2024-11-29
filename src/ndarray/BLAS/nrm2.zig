const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

const scalar = core.supported.scalar;

pub inline fn nrm2(comptime T: type, x: NDArray(T)) !scalar(T) {
    if (x.ndim != 1) {
        return Error.IncompatibleDimensions;
    }

    const supported = core.supported.whatSupportedNumericType(T);

    var res: scalar(T) = 0;
    switch (supported) {
        .BuiltinBool => @compileError("BLAS.nrm2 does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            for (0..x.size) |i| {
                res += x.data[i] * x.data[i];
            }
        },
        .Complex => {
            for (0..x.size) |i| {
                res += x.data[i].re * x.data[i].re + x.data[i].im * x.data[i].im;
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.nrm2 only supports simple types."),
        .Unsupported => unreachable,
    }

    res = @sqrt(res);

    return res;
}

test "nrm2" {
    const a = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(1);

    const result1 = try NDArray(f64).BLAS.nrm2(A.flatten());
    try std.testing.expect(result1 == @sqrt(1008.0));

    const Complex = std.math.Complex;
    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(1, -1));

    const result2: f64 = try NDArray(Complex(f64)).BLAS.nrm2(B.flatten());
    try std.testing.expect(result2 == @sqrt(2016.0));
}
