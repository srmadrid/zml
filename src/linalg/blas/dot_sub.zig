const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");

const blas = @import("../blas.zig");

pub fn dot_sub(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ret: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = types.Coerce(X, Y);

    try ops.set(ret, 0, .{ .allocator = options.allocator });

    if (n <= 0) return blas.Error.InvalidArgument;

    var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
    if (comptime types.isArbitraryPrecision(C)) {
        var temp: C = try ops.init(C, .{ .allocator = options.allocator });
        defer ops.deinit(temp, .{ .allocator = options.allocator });
        for (0..scast(usize, n)) |_| {
            try ops.mul_(
                &temp,
                x[scast(usize, ix)],
                y[scast(usize, iy)],
                .{ .allocator = options.allocator },
            );
            try ops.add_(
                ret,
                ret.*,
                temp,
                .{ .allocator = options.allocator },
            );

            ix += incx;
            iy += incy;
        }
    } else {
        for (0..scast(usize, n)) |_| {
            try ops.add_(
                ret,
                ret.*,
                ops.mul(x[scast(usize, ix)], y[scast(usize, iy)], .{}) catch unreachable,
                .{ .allocator = options.allocator },
            );

            ix += incx;
            iy += incy;
        }
    }
}
