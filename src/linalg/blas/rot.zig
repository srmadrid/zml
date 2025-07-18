const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const blas = @import("../blas.zig");

pub fn rot(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    c: anytype,
    s: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = @TypeOf(c);
    const S: type = @TypeOf(s);
    const Ca: type = types.Coerce(X, types.Coerce(Y, types.Coerce(C, S)));

    if (n <= 0) return blas.Error.InvalidArgument;

    if (ops.eq(c, 1, .{}) catch unreachable and ops.eq(s, 0, .{}) catch unreachable) return;

    var ix: isize = if (incx < 0) (-n + 1) * incx else 0;
    var iy: isize = if (incy < 0) (-n + 1) * incy else 0;
    switch (comptime types.isArbitraryPrecision(Ca)) {
        true => {
            var temp1: Ca = try ops.init(Ca, .{ .allocator = options.allocator });
            defer ops.deinit(&temp1, .{ .allocator = options.allocator });
            var temp2: Ca = try ops.init(Ca, .{ .allocator = options.allocator });
            defer ops.deinit(&temp2, .{ .allocator = options.allocator });
            for (0..scast(usize, n)) |_| {
                try ops.mul_(
                    &temp1,
                    c,
                    x[scast(usize, ix)],
                    .{ .allocator = options.allocator },
                );
                try ops.mul_(
                    &temp2,
                    s,
                    y[scast(usize, iy)],
                    .{ .allocator = options.allocator },
                );
                try ops.add_(
                    &temp1,
                    temp1,
                    temp2,
                    .{ .allocator = options.allocator },
                );

                try ops.mul_(
                    &temp2,
                    s,
                    x[scast(usize, ix)],
                    .{ .allocator = options.allocator },
                );

                try ops.set(
                    &x[scast(usize, ix)],
                    temp1,
                    .{ .allocator = options.allocator },
                );

                try ops.mul_(
                    &temp1,
                    c,
                    y[scast(usize, iy)],
                    .{ .allocator = options.allocator },
                );
                try ops.sub_(
                    &temp1,
                    temp1,
                    temp2,
                    .{ .allocator = options.allocator },
                );
                try ops.set(
                    &y[scast(usize, iy)],
                    temp1,
                    .{ .allocator = options.allocator },
                );

                ix += incx;
                iy += incy;
            }
        },
        false => {
            var temp: Ca = try ops.init(Ca, .{});
            for (0..scast(usize, n)) |_| {
                ops.add_(
                    &temp,
                    ops.mul(c, x[scast(usize, ix)], .{}) catch unreachable,
                    ops.mul(s, y[scast(usize, iy)], .{}) catch unreachable,
                    .{},
                ) catch unreachable;

                ops.sub_(
                    &y[scast(usize, iy)],
                    ops.mul(c, y[scast(usize, iy)], .{}) catch unreachable,
                    ops.mul(s, x[scast(usize, ix)], .{}) catch unreachable,
                    .{},
                ) catch unreachable;

                ops.set(
                    &x[scast(usize, ix)],
                    temp,
                    .{},
                ) catch unreachable;

                ix += incx;
                iy += incy;
            }
        },
    }
}
