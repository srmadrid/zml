const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

pub inline fn copy(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    if (n <= 0) return;

    switch (numericType) {
        .bool, .int, .float, .cfloat => {
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
                    y[@intCast(iy)] = x[@intCast(ix)];
                    y[@intCast(iy + incy)] = x[@intCast(ix + incx)];
                    y[@intCast(iy + incy2)] = x[@intCast(ix + incx2)];
                    y[@intCast(iy + incy3)] = x[@intCast(ix + incx3)];

                    ix += incx4;
                    iy += incy4;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                y[@intCast(iy)] = x[@intCast(ix)];

                ix += incx;
                iy += incy;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.copy only supports simple types."),
    }
}
