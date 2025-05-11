const types = @import("../types.zig");
const Coerce = types.Coerce;

pub fn div(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)) or (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right))))
        @compileError("At least one of left or right must be cfloat, the other must be an int, float, or cfloat");

    switch (types.numericType(@TypeOf(left))) {
        .int => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}).inverse().mulReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), left, .{}));
                },
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}).inverse().mulReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), left, .{}));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(@TypeOf(right))) {
                .int => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).divReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), right, .{}));
                },
                .float => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).divReal(types.cast(types.Scalar(Coerce(@TypeOf(left), @TypeOf(right))), right, .{}));
                },
                .cfloat => {
                    return types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}).div(types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
