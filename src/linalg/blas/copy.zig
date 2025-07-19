const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn copy(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !void {
    if (n <= 0) return blas.Error.InvalidArgument;

    var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
    for (0..scast(usize, n)) |_| {
        try ops.set( // y[iy] = x[ix]
            &y[scast(usize, iy)],
            x[scast(usize, ix)],
            ctx,
        );

        ix += incx;
        iy += incy;
    }
}
