const std = @import("std");

const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");

pub fn acos(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.acos: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    const w: @TypeOf(z) = cfloat.asin(z);
    return .{
        .re = 1.570796326794896619231321691639751442098585 - w.re,
        .im = -w.im,
    };
}
