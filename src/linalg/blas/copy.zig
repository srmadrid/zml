const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn copy(
    n: i32,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    ctx: anytype,
) !void {
    if (n < 0) return blas.Error.InvalidArgument;

    if (n == 0) return;

    var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
    var iy: i32 = if (incy < 0) (-n + 1) * incy else 0;
    var i: u32 = 0;
    while (i < n) : (i += 1) {
        try ops.set( // y[iy] = x[ix]
            &y[scast(u32, iy)],
            x[scast(u32, ix)],
            ctx,
        );

        ix += incx;
        iy += incy;
    }
}
