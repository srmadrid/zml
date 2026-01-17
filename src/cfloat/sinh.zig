const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn sinh(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or !types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.sinh: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return .{
        .re = ops.mul(
            ops.sinh(z.re, .{}) catch unreachable,
            ops.cos(z.im, .{}) catch unreachable,
            .{},
        ) catch unreachable,
        .im = ops.mul(
            ops.cosh(z.re, .{}) catch unreachable,
            ops.sin(z.im, .{}) catch unreachable,
            .{},
        ) catch unreachable,
    };
}
