const types = @import("../types.zig");
const cast = types.cast;
const Scalar = types.Scalar;
const Coerce = types.Coerce;

pub fn mul(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)) or (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right))))
        @compileError("At least one of left or right must be cfloat, the other must be an int, float, or cfloat");

    const R: type = Coerce(@TypeOf(left), @TypeOf(right));

    switch (types.numericType(@TypeOf(left))) {
        .int, .float => {
            switch (types.numericType(@TypeOf(right))) {
                .cfloat => {
                    return cast(R, right, .{}).mulReal(cast(Scalar(R), left, .{}));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(@TypeOf(right))) {
                .int, .float => {
                    return cast(R, left, .{}).mulReal(cast(Scalar(R), right, .{}));
                },
                .cfloat => {
                    return cast(R, left, .{}).mul(cast(R, right, .{}));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
