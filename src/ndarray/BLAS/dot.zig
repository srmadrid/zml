const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

pub inline fn dot(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return 0.0;

    var temp: T = 0;
    switch (supported) {
        .BuiltinBool => @compileError("BLAS.dot does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (incx == 1 and incy == 1) {
                const m: usize = @as(usize, @intCast(n)) % 5;
                for (0..m) |i| {
                    temp += x[i] * y[i];
                }

                if (n < 5) return temp;

                var mp1: usize = m;
                while (mp1 < n) : (mp1 += 5) {
                    temp += x[mp1] * y[mp1] + x[mp1 + 1] * y[mp1 + 1] + x[mp1 + 2] * y[mp1 + 2] + x[mp1 + 3] * y[mp1 + 3] + x[mp1 + 4] * y[mp1 + 4];
                }
            } else {
                var ix: isize = 0;
                var iy: isize = 0;

                if (incx < 0) ix = (-@as(isize, @intCast(n)) + 1) * incx;
                if (incy < 0) iy = (-@as(isize, @intCast(n)) + 1) * incy;

                for (0..@intCast(n)) |_| {
                    temp += x[@intCast(ix)] * y[@intCast(iy)];
                    ix += incx;
                    iy += incy;
                }
            }
        },
        .Complex => @compileError("BLAS.dot does not support complex numbers. Use BLAS.dotc or BLAS.dotu instead."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.dot only supports simple types."),
        .Unsupported => unreachable,
    }

    return temp;
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
