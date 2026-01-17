const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");
const constants = @import("../constants.zig");

pub fn acosh(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.acosh: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return cfloat.acos(z).mulImag(constants.one(types.Scalar(Z), .{}) catch unreachable);
}
