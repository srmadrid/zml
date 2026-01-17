const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn Abs(comptime Z: type) type {
    comptime if (!types.isNumeric(Z) or !types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.abs: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return types.Scalar(Z);
}

pub fn abs(z: anytype) Abs(@TypeOf(z)) {
    return ops.hypot(z.re, z.im, .{}) catch unreachable;
}
