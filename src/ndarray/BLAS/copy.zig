const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

pub inline fn copy(comptime T: type, x: NDArray(T), y: *NDArray(T)) !void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (x.size != y.size) {
        return Error.IncompatibleSize;
    }

    switch (supported) {
        .BuiltinBool, .BuiltinInt, .BuiltinFloat, .Complex => {
            for (0..x.size) |i| {
                y.data[i] = x.data[i];
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.copy only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "copy" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer A.deinit();

    A.setAll(1);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer B.deinit();

    try NDArray(f64).BLAS.copy(A, &B);

    for (0..B.size) |i| {
        try std.testing.expect(B.data[i] == 1);
    }

    const Complex = std.math.Complex;
    var C: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer C.deinit();

    C.setAll(Complex(f64).init(1, -1));

    var D: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer D.deinit();

    try NDArray(Complex(f64)).BLAS.copy(C, &D);

    for (0..D.size) |i| {
        try std.testing.expect(std.math.approxEqAbs(f64, D.data[i].re, 1, 0.0001) and std.math.approxEqAbs(f64, D.data[i].im, -1, 0.0001));
    }
}
