const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");

pub fn asinh(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.asinh: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return cfloat.asin(
        z.mulImag(
            constants.one(types.Scalar(Z), .{}) catch unreachable,
        ),
    ).mulImag(
        ops.neg(
            constants.one(types.Scalar(Z), .{}) catch unreachable,
            .{},
        ) catch unreachable,
    );
}
