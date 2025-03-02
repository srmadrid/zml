const std = @import("std");
const core = @import("../../core/core.zig");
const blas = @import("../blas.zig");

const scalar = core.supported.scalar;

pub inline fn iamin(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0 or incx <= 0) return 0;

    if (n == 1) return 0;

    var imin: usize = std.math.maxInt(usize);
    var i: usize = 0;

    switch (supported) {
        .BuiltinBool => @compileError("blas.iamin does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            var min: T = std.math.floatMax(T);

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
                    var absx = @abs(x[@intCast(ix)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx2)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx3)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx4)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx5)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx6)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx7)]);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    ix += incx8;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const absx = @abs(x[@intCast(ix)]);
                if (absx < min) {
                    min = absx;
                    imin = i;
                }
                i += 1;

                ix += incx;
            }
        },
        .Complex => {
            var min: scalar(T) = std.math.floatMax(scalar(T));

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
                    var absx = @abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx)].re) + @abs(x[@intCast(ix + incx)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx2)].re) + @abs(x[@intCast(ix + incx2)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx3)].re) + @abs(x[@intCast(ix + incx3)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx4)].re) + @abs(x[@intCast(ix + incx4)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx5)].re) + @abs(x[@intCast(ix + incx5)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx6)].re) + @abs(x[@intCast(ix + incx6)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx7)].re) + @abs(x[@intCast(ix + incx7)].im);
                    if (absx < min) {
                        min = absx;
                        imin = i;
                    }
                    i += 1;

                    ix += incx8;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const absx = @abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im);
                if (absx < min) {
                    min = absx;
                    imin = i;
                }
                i += 1;

                ix += incx;
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.iamin only supports simple types."),
        .Unsupported => unreachable,
    }

    return imin;
}

test "iamin" {
    const a: std.mem.Allocator = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);

    for (0..n) |i| {
        x1[i] = std.math.floatMax(f64);
    }

    x1[127] = 0;
    x1[456] = 0;

    const result1 = blas.iamin(f64, n, x1.ptr, 1);
    try std.testing.expectEqual(127, result1);
    const result2 = blas.iamin(f64, n, x1.ptr, -1);
    try std.testing.expectEqual(0, result2);
    const result3 = blas.iamin(f64, n / 2, x1.ptr, 2);
    try std.testing.expectEqual(228, result3);

    var x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);

    for (0..n) |i| {
        x2[i] = Complex(f64).init(std.math.floatMax(f64), std.math.floatMax(f64));
    }

    x2[127] = Complex(f64).init(0, 0);
    x2[456] = Complex(f64).init(0, 0);

    const result4 = blas.iamin(Complex(f64), n, x2.ptr, 1);
    try std.testing.expectEqual(127, result4);
    const result5 = blas.iamin(Complex(f64), n, x2.ptr, -1);
    try std.testing.expectEqual(0, result5);
    const result6 = blas.iamin(Complex(f64), n / 2, x2.ptr, 2);
    try std.testing.expectEqual(228, result6);
}
