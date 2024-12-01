const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn swap(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.swap does not support bool."),
        .BuiltinInt, .BuiltinFloat, .Complex => {
            if (incx == 1 and incy == 1) {
                const m: usize = @as(usize, @intCast(n)) % 3;
                for (0..m) |i| {
                    const temp = x[i];
                    x[i] = y[i];
                    y[i] = temp;
                }

                if (n < 3) return;

                var i: usize = m;
                while (i < n) : (i += 3) {
                    var temp = x[i];
                    x[i] = y[i];
                    y[i] = temp;
                    temp = x[i + 1];
                    x[i + 1] = y[i + 1];
                    y[i + 1] = temp;
                    temp = x[i + 2];
                    x[i + 2] = y[i + 2];
                    y[i + 2] = temp;
                }
            } else {
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..@intCast(n)) |_| {
                    const temp = x[@intCast(ix)];
                    x[@intCast(ix)] = y[@intCast(iy)];
                    y[@intCast(iy)] = temp;
                    ix += incx;
                    iy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.swap only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "swap" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer A.deinit();

    A.setAll(1);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 2, 3, 4 }, .{});
    defer B.deinit();

    B.setAll(2);

    try NDArray(f64).BLAS.swap(@constCast(&A.flatten()), @constCast(&B.flatten()));

    for (0..A.size) |i| {
        try std.testing.expect(A.data[i] == 2);
        try std.testing.expect(B.data[i] == 1);
    }

    const Complex = std.math.Complex;
    var C: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer C.deinit();

    C.setAll(Complex(f64).init(1, -1));

    var D: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer D.deinit();

    D.setAll(Complex(f64).init(2, 1));

    try NDArray(Complex(f64)).BLAS.swap(@constCast(&C.flatten()), @constCast(&D.flatten()));

    for (0..B.size) |i| {
        try std.testing.expect(std.math.approxEqAbs(f64, C.data[i].re, 2, 0.0001) and std.math.approxEqAbs(f64, C.data[i].im, 1, 0.0001));
        try std.testing.expect(std.math.approxEqAbs(f64, D.data[i].re, 1, 0.0001) and std.math.approxEqAbs(f64, D.data[i].im, -1, 0.0001));
    }
}
