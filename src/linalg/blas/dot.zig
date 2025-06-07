const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

pub inline fn dot(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return 0;

    var sum: T = 0;
    switch (numericType) {
        .bool => @compileError("blas.dot does not support bool."),
        .int, .float => {
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 2) << 2;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incx3 = incx * 3;
                const incx4 = incx * 4;
                const incy2 = incy * 2;
                const incy3 = incy * 3;
                const incy4 = incy * 4;
                while (ix != StX) {
                    sum += x[@intCast(ix)] * y[@intCast(iy)] + x[@intCast(ix + incx)] * y[@intCast(iy + incy)] + x[@intCast(ix + incx2)] * y[@intCast(iy + incy2)] + x[@intCast(ix + incx3)] * y[@intCast(iy + incy3)];

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                sum += x[@intCast(ix)] * y[@intCast(iy)];

                ix += incx;
                iy += incy;
            }
        },
        .cfloat => @compileError("blas.dot does not support complex numbers. Use blas.dotc or blas.dotu instead."),
        .integer, .rational, .real, .complex, .expression => @compileError("blas.dot only supports simple types."),
    }

    return sum;
}
