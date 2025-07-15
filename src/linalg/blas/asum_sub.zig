const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const float = @import("../../float.zig");

const blas = @import("../blas.zig");

pub fn asum_sub(
    n: isize,
    x: anytype,
    incx: isize,
    ret: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    const X: type = types.Child(@TypeOf(x));

    try ops.set(ret, 0, .{ .allocator = options.allocator });

    if (n <= 0 or incx <= 0) return blas.Error.InvalidArgument;

    var ix: isize = 0;
    switch (comptime types.numericType(X)) {
        .cfloat => {
            for (0..scast(usize, n)) |_| {
                try ops.add_(ret, ret.*, float.abs(x[scast(usize, ix)].re) + float.abs(x[scast(usize, ix)].im), .{ .allocator = options.allocator });

                ix += incx;
            }
        },
        .complex => {
            var temp: Scalar(X) = try ops.init(Scalar(X), .{ .allocator = options.allocator });
            defer ops.deinit(temp, .{ .allocator = options.allocator });
            for (0..scast(usize, n)) |_| {
                try ops.add_(
                    &temp,
                    ops.abs(x[scast(usize, ix)].re, .{ .copy = false }) catch unreachable,
                    ops.abs(x[scast(usize, ix)].im, .{ .copy = false }) catch unreachable,
                    .{ .allocator = options.allocator },
                );
                try ops.add_(ret, ret.*, temp, .{ .allocator = options.allocator });

                ix += incx;
            }
        },
        else => {
            for (0..scast(usize, n)) |_| {
                try ops.add_(ret, ret.*, ops.abs(x[scast(usize, ix)], .{ .copy = false }) catch unreachable, .{ .allocator = options.allocator });

                ix += incx;
            }
        },
    }
}
