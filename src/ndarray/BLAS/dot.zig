const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

pub inline fn dot(comptime T: type, x: NDArray(T), y: NDArray(T)) !T {
    if (x.ndim != 1 or y.ndim != 1) {
        return Error.IncompatibleDimensions;
    }

    if (x.shape[0] != y.shape[0]) {
        return Error.IncompatibleDimensions;
    }

    const supported = core.supported.whatSupportedNumericType(T);

    var res: T = undefined;
    switch (supported) {
        .BuiltinBool => @compileError("BLAS.dot does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            res = 0;
            for (0..x.size) |i| {
                res += x.data[i] * y.data[i];
            }
        },
        .Complex => @compileError("BLAS.dot does not support complex numbers. Use BLAS.dotc or BLAS.dotu instead."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.dot only supports simple types."),
        .Unsupported => unreachable,
    }

    return res;
}

test "dot" {
    const a = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(1);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(2);

    const result1 = try NDArray(f64).BLAS.dot(A.flatten(), B.flatten());
    try std.testing.expect(result1 == 2016);
}
