const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn Abs2(comptime Z: type) type {
    comptime if (!types.isNumeric(Z) or !types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.abs2: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return types.Scalar(Z);
}

pub fn abs2(z: anytype) Abs2(@TypeOf(z)) {
    return ops.add(
        ops.mul(
            z.re,
            z.re,
            .{},
        ) catch unreachable,
        ops.mul(
            z.im,
            z.im,
            .{},
        ) catch unreachable,
        .{},
    ) catch unreachable;
}
