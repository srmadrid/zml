const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

pub inline fn rot(comptime T: type, x: *NDArray(T), y: *NDArray(T), c: T, s: T) !void {
    if (x.ndim != 1 or y.ndim != 1) {
        return Error.IncompatibleDimensions;
    }

    if (x.shape[0] != y.shape[0]) {
        return Error.IncompatibleDimensions;
    }

    const supported = core.supported.whatSupportedNumericType(T);

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.nrm2 does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            for (0..x.size) |i| {
                const temp = c * x.data[i] + s * y.data[i];
                y.data[i] = -s * x.data[i] + c * y.data[i];
                x.data[i] = temp;
            }
        },
        .Complex => {
            for (0..x.size) |i| {
                const temp = T.add(T.mul(c, x.data[i]), T.mul(s, y.data[i]));
                y.data[i] = T.sub(T.mul(c, y.data[i]), T.mul(T.conjugate(s), x.data[i]));
                x.data[i] = temp;
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

    try NDArray(Complex(f64)).BLAS.rot(@constCast(&C.flatten()), @constCast(&D.flatten()), Complex(f64).init(0.7071067811865475, 0.7071067811865475), Complex(f64).init(0.7071067811865475, -0.7071067811865475));

    for (0..C.size) |i| {
        try std.testing.expectApproxEqAbs(0.7071067811865475 * 2, C.data[i].re, 0.0001);
        try std.testing.expectApproxEqAbs(-0.7071067811865475 * 2, C.data[i].im, 0.0001);
        try std.testing.expectApproxEqAbs(0, D.data[i].re, 0.0001);
        try std.testing.expectApproxEqAbs(0, D.data[i].im, 0.0001);
    }
}
