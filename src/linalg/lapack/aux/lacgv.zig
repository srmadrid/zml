const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn lacgv(
    n: i32,
    x: anytype,
    incx: i32,
) !void {
    if (incx == 1) {
        var i: u32 = 0;
        while (i < types.scast(u32, n)) : (i += 1) {
            try ops.conj_( // x[i] = conj(x[i])
                &x[i],
                x[i],
                if (comptime types.isArbitraryPrecision(types.Child(@TypeOf(x)))) .{ .copy = false } else .{},
            );
        }
    } else {
        var ix: i32 = if (incx < 0) (-n + 1) * incx else 0;
        for (0..types.scast(u32, n)) |_| {
            try ops.conj_( // x[ix] = conj(x[ix])
                &x[types.scast(u32, ix)],
                x[types.scast(u32, ix)],
                if (comptime types.isArbitraryPrecision(types.Child(@TypeOf(x)))) .{ .copy = false } else .{},
            );

            ix += incx;
        }
    }
}
