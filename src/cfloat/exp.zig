const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn exp(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or !types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.exp: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    const r: @TypeOf(z.re) = ops.exp(z.re, .{}) catch unreachable;
    return .{
        .re = ops.mul(
            r,
            ops.cos(z.im, .{}) catch unreachable,
            .{},
        ) catch unreachable,
        .im = ops.mul(
            r,
            ops.sin(z.im, .{}) catch unreachable,
            .{},
        ) catch unreachable,
    };
}
