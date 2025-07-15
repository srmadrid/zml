const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn axpy(
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C: type = types.Coerce(Al, X);

    if (n <= 0) return blas.Error.InvalidArgument;

    if (ops.eq(alpha, 0, .{}) catch unreachable) return;

    var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
    switch (comptime types.isArbitraryPrecision(C)) {
        true => {
            var temp: C = try ops.init(C, .{ .allocator = options.allocator });
            defer ops.deinit(temp, .{ .allocator = options.allocator });
            for (0..scast(usize, n)) |_| {
                try ops.mul_(
                    &temp,
                    alpha,
                    x[scast(usize, ix)],
                    .{ .allocator = options.allocator },
                );
                try ops.add_(
                    &y[scast(usize, iy)],
                    y[scast(usize, iy)],
                    temp,
                    .{ .allocator = options.allocator },
                );

                ix += incx;
                iy += incy;
            }
        },
        false => {
            for (0..scast(usize, n)) |_| {
                try ops.add_(
                    &y[scast(usize, iy)],
                    y[scast(usize, iy)],
                    ops.mul(alpha, x[scast(usize, ix)], .{}) catch unreachable,
                    .{ .allocator = options.allocator },
                );

                ix += incx;
                iy += incy;
            }
        },
    }
}
