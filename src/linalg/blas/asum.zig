const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

const Scalar = types.Scalar;

pub inline fn asum(comptime T: type, n: isize, x: [*]const T, incx: isize) Scalar(T) {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0 or incx < 0) return 0;

    var sum: Scalar(T) = 0;

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
    }

    return sum;
}
