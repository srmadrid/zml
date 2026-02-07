const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");

pub fn check_aliasing(o: *const integer.Integer, x: anytype) bool {
    const X: type = if (comptime types.isPointer(@TypeOf(x))) types.Child(@TypeOf(x)) else @TypeOf(x);

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .dyadic, .cfloat => return false,
        .integer => {
            const ostart: usize = @intFromPtr(o.limbs);
            const oend: usize = ostart + o.size * 4;
            const xstart: usize = @intFromPtr(x.limbs);
            const xend: usize = xstart + x.size * 4;

            return ostart < xend and xstart < oend;
        },
        else => unreachable,
    }
}
