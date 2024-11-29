const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const Error = @import("../ndarray.zig").Error;
const core = @import("../../core/core.zig");

pub inline fn copy(comptime T: type, n: usize, x: [*]T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (n == 0) return;

    switch (supported) {
        .BuiltinBool, .BuiltinInt, .BuiltinFloat, .Complex => {
            if (incx == 1 and incy == 1) {
                // Code for both increments equal to 1
                const m: usize = n % 7;

                // Clean-up loop for the remainder
                for (0..m) |i| {
                    y[i] = x[i];
                }

                if (n < 7) return;

                // Process in blocks of 7 for efficiency
                var i: usize = m;
                while (i < n) : (i += 7) {
                    y[i] = x[i];
                    y[i + 1] = x[i + 1];
                    y[i + 2] = x[i + 2];
                    y[i + 3] = x[i + 3];
                    y[i + 4] = x[i + 4];
                    y[i + 5] = x[i + 5];
                    y[i + 6] = x[i + 6];
                }
            } else {
                // Code for unequal increments or increments not equal to 1
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..n) |_| {
                    y[@intCast(iy)] = x[@intCast(ix)];
                    ix += incx;
                    iy += incy;
                }
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

    try NDArray(f64).BLAS.copy(A.flatten(), @constCast(&B.flatten()));

    for (0..B.size) |i| {
        try std.testing.expect(B.data[i] == 1);
    }

    const Complex = std.math.Complex;
    var C: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer C.deinit();

    C.setAll(Complex(f64).init(1, -1));

    var D: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 2, 3, 4 }, .{});
    defer D.deinit();

    try NDArray(Complex(f64)).BLAS.copy(C.flatten(), @constCast(&D.flatten()));

    for (0..D.size) |i| {
        try std.testing.expect(std.math.approxEqAbs(f64, D.data[i].re, 1, 0.0001) and std.math.approxEqAbs(f64, D.data[i].im, -1, 0.0001));
    }
}
