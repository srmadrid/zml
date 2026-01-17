const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn neg(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.neg: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return .{
        .re = ops.neg(z.re, .{}) catch unreachable,
        .im = ops.neg(z.im, .{}) catch unreachable,
    };
}
