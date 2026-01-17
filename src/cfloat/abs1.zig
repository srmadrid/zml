const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn Abs1(comptime Z: type) type {
    comptime if (!types.isNumeric(Z) or types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.abs1: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return types.Scalar(Z);
}

pub fn abs1(z: anytype) Abs1(@TypeOf(z)) {
    return ops.add(
        ops.abs(z.re, .{}) catch unreachable,
        ops.abs(z.im, .{}) catch unreachable,
        .{},
    ) catch unreachable;
}
