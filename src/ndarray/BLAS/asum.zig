const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

const scalar = core.supported.scalar;

pub inline fn asum(comptime T: type, x: NDArray(T)) scalar(T) {
    const supported = core.supported.whatSupportedNumericType(T);

    var res: scalar(T) = 0;
    switch (supported) {
        .BuiltinBool => @compileError("BLAS.asum does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            for (0..x.size) |i| {
                res += @abs(x.data[i]);
            }
        },
        .Complex => {
            for (0..x.size) |i| {
                res += @abs(x.data[i].re) + @abs(x.data[i].im);
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.asum only supports simple types."),
        .Unsupported => unreachable,
    }

    return res;
}

test "asum" {
    const a = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(1);

    const result1 = NDArray(f64).BLAS.asum(A);
    try std.testing.expect(result1 == 1008);

    const Complex = std.math.Complex;
    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(1, -1));

    const result2: f64 = NDArray(Complex(f64)).BLAS.asum(B);
    try std.testing.expect(result2 == 2016);
}
