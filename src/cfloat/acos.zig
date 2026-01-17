const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");

pub fn acos(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or !types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.acos: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    const w: @TypeOf(z) = cfloat.asin(z);
    return .{
        .re = 1.570796326794896619231321691639751442098585 - w.re,
        .im = -w.im,
    };
}
