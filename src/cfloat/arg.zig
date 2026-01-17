const types = @import("../types.zig");
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");

pub fn Arg(comptime Z: type) type {
    comptime if (!types.isNumeric(Z) or !types.numericType(Z).le(.cfloat))
        @compileError("zml.cfloat.arg: z must be a bool, an int, a float or a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return types.EnsureFloat(types.Scalar(Z));
}

pub fn arg(z: anytype) Arg(@TypeOf(z)) {
    switch (types.numericType(@TypeOf(z))) {
        .bool, .int, .float, .dyadic => {
            return if (ops.ge(z, 0, .{}) catch unreachable)
                constants.zero(Arg(@TypeOf(z)), .{}) catch unreachable
            else
                constants.pi(Arg(@TypeOf(z)), .{}) catch unreachable;
        },
        .cfloat => return ops.atan2(z.im, z.re, .{}) catch unreachable,
        else => unreachable,
    }
}
