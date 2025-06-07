const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../blas.zig");

pub inline fn dotu_sub(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize, ret: *T) void {
    @setRuntimeSafety(false);
    const numericType = types.numericType(T);

    ret.* = T.init(0, 0);

    if (n <= 0) return;

    switch (numericType) {
        .bool => @compileError("blas.dotu_sub does not support bool."),
        .int => @compileError("blas.dotu_sub does not support integers. Use blas.dot instead."),
        .float => @compileError("blas.dotu_sub does not support floats. Use blas.dot instead."),
        .cfloat => {
            var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
            var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
            const nu = (n >> 1) << 1;
            if (nu != 0) {
                const StX = ix + nu * incx;
                const incx2 = incx * 2;
                const incy2 = incy * 2;
                while (ix != StX) {
                    ret.*.re += x[@intCast(ix)].re * y[@intCast(iy)].re - x[@intCast(ix)].im * y[@intCast(iy)].im + x[@intCast(ix + incx)].re * y[@intCast(iy + incy)].re - x[@intCast(ix + incx)].im * y[@intCast(iy + incy)].im;
                    ret.*.im += x[@intCast(ix)].re * y[@intCast(iy)].im + x[@intCast(ix)].im * y[@intCast(iy)].re + x[@intCast(ix + incx)].re * y[@intCast(iy + incy)].im + x[@intCast(ix + incx)].im * y[@intCast(iy + incy)].re;

                    ix += incx2;
                    iy += incy2;
                }
            }

            for (@intCast(nu)..@intCast(n)) |_| {
                ret.*.re += x[@intCast(ix)].re * y[@intCast(iy)].re - x[@intCast(ix)].im * y[@intCast(iy)].im;
                ret.*.im += x[@intCast(ix)].re * y[@intCast(iy)].im + x[@intCast(ix)].im * y[@intCast(iy)].re;

                ix += incx;
                iy += incy;
            }
        },
        .integer, .rational, .real, .complex, .expression => @compileError("blas.dotu only supports simple types."),
    }

    return;
}
