const std = @import("std");
const NDArray = @import("../ndarray.zig").NDArray;
const core = @import("../../core/core.zig");

const scalar = core.supported.scalar;

pub inline fn asum(comptime T: type, n: usize, x: [*]T, incx: isize) scalar(T) {
    const supported = core.supported.whatSupportedNumericType(T);

    if (n == 0 or incx == 0) {
        return 0;
    }

    var temp: scalar(T) = 0;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.asum does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (incx == 1) {
                // Code for increment equal to 1.
                const m: usize = n % 6;

                // Clean-up loop for elements not divisible by 6.
                for (0..m) |i| {
                    temp += @abs(x[i]);
                }

                if (n < 6) {
                    return temp;
                }

                // Main loop, unrolling in chunks of 6.
                var i: usize = m;
                while (i < n) : (i += 6) {
                    temp += @abs(x[i]) + @abs(x[i + 1]) +
                        @abs(x[i + 2]) + @abs(x[i + 3]) +
                        @abs(x[i + 4]) + @abs(x[i + 5]);
                }
            } else {
                // Code for increment not equal to 1.
                var idx: isize = 0;
                for (0..n) |_| {
                    temp += @abs(x[@intCast(idx)]);
                    idx += incx;
                }
            }
        },
        .Complex => {
            if (incx == 1) {
                // Code for increment equal to 1.
                const m: usize = n % 6;

                // Clean-up loop for elements not divisible by 6.
                for (0..m) |i| {
                    temp += @abs(x[i].re) + @abs(x[i].im);
                }

                if (n < 6) {
                    return temp;
                }

                // Main loop, unrolling in chunks of 6.
                var i: usize = m;
                while (i < n) : (i += 6) {
                    temp += @abs(x[i].re) + @abs(x[i].im) + @abs(x[i + 1].re) + @abs(x[i + 1].im) +
                        @abs(x[i + 2].re) + @abs(x[i + 2].im) + @abs(x[i + 3].re) + @abs(x[i + 3].im) +
                        @abs(x[i + 4].re) + @abs(x[i + 4].im) + @abs(x[i + 5].re) + @abs(x[i + 5].im);
                }
            } else {
                // Code for increment not equal to 1.
                var idx: isize = 0;
                for (0..n) |_| {
                    temp += @abs(x[@intCast(idx)].re) + @abs(x[@intCast(idx)].im);
                    idx += incx;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.asum only supports simple types."),
        .Unsupported => unreachable,
    }

    return temp;
}

test "asum" {
    const a = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 8, 18, 7 }, .{});
    defer A.deinit();

    A.setAll(1);

    const result1 = try NDArray(f64).BLAS.asum(A.flatten());
    try std.testing.expect(result1 == 1008);

    const Complex = std.math.Complex;
    var B: NDArray(Complex(f64)) = try NDArray(Complex(f64)).init(a, &.{ 8, 18, 7 }, .{});
    defer B.deinit();

    B.setAll(Complex(f64).init(1, -1));

    const result2: f64 = try NDArray(Complex(f64)).BLAS.asum(B.flatten());
    try std.testing.expect(result2 == 2016);
}
