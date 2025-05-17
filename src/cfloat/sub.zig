const types = @import("../types.zig");
const cast = types.cast;
const Scalar = types.Scalar;
const Coerce = types.Coerce;

pub fn sub(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)) or (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right))))
        @compileError("At least one of left or right must be cfloat, the other must be an int, float, or cfloat");

    const R: type = Coerce(@TypeOf(left), @TypeOf(right));

    switch (types.numericType(@TypeOf(left))) {
        .int, .float => {
            switch (types.numericType(@TypeOf(right))) {
                .cfloat => {
                    return .{
                        .re = cast(Scalar(R), left, .{}) - cast(Scalar(R), right.re, .{}),
                        .im = -cast(Scalar(R), right.im, .{}),
                    };
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(@TypeOf(right))) {
                .int, .float => {
                    return cast(R, left, .{}).subReal(cast(Scalar(R), right, .{}));
                },
                .cfloat => {
                    return cast(R, left, .{}).sub(cast(R, right, .{}));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
