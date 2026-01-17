const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");
const ops = @import("../ops.zig");

pub fn Log(comptime Z: type) type {
    comptime if (!types.isNumeric(Z) or !types.numericType(Z).le(.cfloat))
        @compileError("zml.cfloat.log: z must be a bool, an int, a float or a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return cfloat.Cfloat(types.EnsureFloat(types.Scalar(Z)));
}

pub fn log(z: anytype) Log(@TypeOf(z)) {
    const zz = types.scast(Log(@TypeOf(z)), z);

    var p: types.Scalar(@TypeOf(zz)) = cfloat.abs(zz);
    p = ops.log(p, .{}) catch unreachable;
    return .{
        .re = p,
        .im = ops.atan2(zz.im, zz.re, .{}) catch unreachable,
    };
}
