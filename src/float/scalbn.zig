const std = @import("std");

const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const float = @import("../float.zig");

const dbl64 = @import("dbl64.zig");
const ldbl128 = @import("ldbl128.zig");

pub inline fn scalbn(x: anytype, n: i32) EnsureFloat(@TypeOf(x)) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.scalbn: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    return std.math.scalbn(x, n);
}
