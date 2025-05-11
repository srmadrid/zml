const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");
const options = @import("options");

pub inline fn scal(comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    switch (numericType) {
        .bool => @compileError("blas.scal does not support bool."),
        .int, .float => {
            if (alpha == 1) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
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
                    x[@intCast(ix)] *= alpha;
                    x[@intCast(ix + incx)] *= alpha;
                    x[@intCast(ix + incx2)] *= alpha;
                    x[@intCast(ix + incx3)] *= alpha;
                    x[@intCast(ix + incx4)] *= alpha;
                    x[@intCast(ix + incx5)] *= alpha;
                    x[@intCast(ix + incx6)] *= alpha;
                    x[@intCast(ix + incx7)] *= alpha;

                    ix += incx8;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                x[@intCast(ix)] *= alpha;

                ix += incx;
            }
        },
        .cfloat => {
            if (alpha.re == 1 and alpha.im == 0) return;

            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            const nu = (n >> 2) << 2;
            if (alpha.re == 0 and alpha.im == 0) {
                if (nu != 0) {
                    const StX = ix + nu * incx;
                    const incx2 = incx * 2;
                    const incx3 = incx * 3;
                    const incx4 = incx * 4;
                    while (ix != StX) {
                        x[@intCast(ix)].re = 0;
                        x[@intCast(ix)].im = 0;
                        x[@intCast(ix + incx)].re = 0;
                        x[@intCast(ix + incx)].im = 0;
                        x[@intCast(ix + incx2)].re = 0;
                        x[@intCast(ix + incx2)].im = 0;
                        x[@intCast(ix + incx3)].re = 0;
                        x[@intCast(ix + incx3)].im = 0;

                        ix += incx4;
                    }

                    for (@intCast(nu)..@intCast(n)) |_| {
                        x[@intCast(ix)].re = 0;
                        x[@intCast(ix)].im = 0;

                        ix += incx;
                    }
                }
            } else if (alpha.re != 1 or alpha.im != 0) {
                if (nu != 0) {
                    const StX = ix + nu * incx;
                    const incx2 = incx * 2;
                    const incx3 = incx * 3;
                    const incx4 = incx * 4;
                    while (ix != StX) {
                        const x0 = x[@intCast(ix)];
                        const x1 = x[@intCast(ix + incx)];
                        const x2 = x[@intCast(ix + incx2)];
                        const x3 = x[@intCast(ix + incx3)];

                        x[@intCast(ix)].re = x0.re * alpha.re - x0.im * alpha.im;
                        x[@intCast(ix)].im = x0.re * alpha.im + x0.im * alpha.re;
                        x[@intCast(ix + incx)].re = x1.re * alpha.re - x1.im * alpha.im;
                        x[@intCast(ix + incx)].im = x1.re * alpha.im + x1.im * alpha.re;
                        x[@intCast(ix + incx2)].re = x2.re * alpha.re - x2.im * alpha.im;
                        x[@intCast(ix + incx2)].im = x2.re * alpha.im + x2.im * alpha.re;
                        x[@intCast(ix + incx3)].re = x3.re * alpha.re - x3.im * alpha.im;
                        x[@intCast(ix + incx3)].im = x3.re * alpha.im + x3.im * alpha.re;

                        ix += incx4;
                    }

                    for (@intCast(nu)..@intCast(n)) |_| {
                        const x0 = x[@intCast(ix)];
                        x[@intCast(ix)].re = x0.re * alpha.re - x0.im * alpha.im;
                        x[@intCast(ix)].im = x0.re * alpha.im + x0.im * alpha.re;

                        ix += incx;
                    }
                }
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.scal only supports simple types."),
        .unsupported => unreachable,
    }
}
