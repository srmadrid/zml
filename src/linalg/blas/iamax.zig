const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

const Scalar = types.Scalar;

pub inline fn iamax(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0 or incx <= 0) return 0;

    if (n == 1) return 0;

    var imax: usize = 0;
    var i: usize = 0;

    switch (numericType) {
        .bool => @compileError("blas.iamax does not support bool."),
        .int, .float => {
            var max: T = 0;

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
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx)]);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx2)]);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx3)]);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx4)]);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx5)]);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx6)]);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx7)]);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    ix += incx8;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const absx = @abs(x[@intCast(ix)]);
                if (absx > max) {
                    max = absx;
                    imax = i;
                }
                i += 1;

                ix += incx;
            }
        },
        .cfloat => {
            var max: Scalar(T) = 0;

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
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx)].re) + @abs(x[@intCast(ix + incx)].im);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx2)].re) + @abs(x[@intCast(ix + incx2)].im);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx3)].re) + @abs(x[@intCast(ix + incx3)].im);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx4)].re) + @abs(x[@intCast(ix + incx4)].im);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx5)].re) + @abs(x[@intCast(ix + incx5)].im);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx6)].re) + @abs(x[@intCast(ix + incx6)].im);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    absx = @abs(x[@intCast(ix + incx7)].re) + @abs(x[@intCast(ix + incx7)].im);
                    if (absx > max) {
                        max = absx;
                        imax = i;
                    }
                    i += 1;

                    ix += incx8;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                const absx = @abs(x[@intCast(ix)].re) + @abs(x[@intCast(ix)].im);
                if (absx > max) {
                    max = absx;
                    imax = i;
                }
                i += 1;

                ix += incx;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.iamax only supports simple types."),
        .unsupported => unreachable,
    }

    return imax;
}
