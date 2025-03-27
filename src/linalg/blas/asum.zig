const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");

const Numeric = core.types.Numeric;

pub inline fn asum(comptime T: type, n: isize, x: [*]const T, incx: isize) Numeric(T) {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0 or incx < 0) return 0;

    var sum: Numeric(T) = 0;

    switch (numericType) {
        .bool => @compileError("blas.asum does not support bool."),
        .int, .float => {
            var ix: isize = 0;
            const nu = (n >> 3) << 3;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incx3 = incx * 3;
                const incx4 = incx * 4;
                const incx5 = incx * 5;
                const incx6 = incx * 6;
                const incx7 = incx * 7;
                const incx8 = incx * 8;
                while (ix != StX) {
                    sum += @abs(x[@intCast(ix)]) + @abs(x[@intCast(ix + incx)]) + @abs(x[@intCast(ix + incx2)]) + @abs(x[@intCast(ix + incx3)]) + @abs(x[@intCast(ix + incx4)]) + @abs(x[@intCast(ix + incx5)]) + @abs(x[@intCast(ix + incx6)]) + @abs(x[@intCast(ix + incx7)]);

                    ix += incx8;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                sum += @abs(x[@intCast(ix)]);

                ix += incx;
            }
        },
        .cfloat => {
            var ix: isize = 0;
            const nu = (n >> 2) << 2;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incx3 = incx * 3;
                const incx4 = incx * 4;
                while (ix != StX) {
                    sum += @abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im) + @abs(x[@intCast(ix + incx)].re) + @abs(x[@intCast(ix + incx)].im) + @abs(x[@intCast(ix + incx2)].re) + @abs(x[@intCast(ix + incx2)].im) + @abs(x[@intCast(ix + incx3)].re) + @abs(x[@intCast(ix + incx3)].im);

                    ix += incx4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                sum += @abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im);

                ix += incx;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.asum only supports simple types."),
        .unsupported => unreachable,
    }

    return sum;
}

test asum {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
    }

    const result1 = blas.asum(f64, n, x1.ptr, 1);
    try std.testing.expectEqual(500500, result1);
    const result2 = blas.asum(f64, n / 2, x1.ptr, 2);
    try std.testing.expectEqual(250000, result2);

    var x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);

    for (0..n) |i| {
        x2[i] = Complex(f64).init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
    }

    const result3 = blas.asum(Complex(f64), n, x2.ptr, 1);
    try std.testing.expectEqual(1001000, result3);
    const result4 = blas.asum(Complex(f64), n / 2, x2.ptr, 2);
    try std.testing.expectEqual(500000, result4);
}
