const types = @import("../types.zig");
const float = @import("../float.zig");

pub fn neg(z: anytype) @TypeOf(z) {
    comptime if (types.numericType(@TypeOf(z)) != .cfloat)
        @compileError("cfloat.abs: z must be a cfloat, got " ++ @typeName(@TypeOf(z)));

    return .{ .re = -z.re, .im = -z.im };
}
