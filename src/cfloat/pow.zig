const cfloat = @import("../cfloat.zig");
const types = @import("../types.zig");
const Coerce = types.Coerce;
const cast = types.cast;

pub fn pow(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)) or (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right))))
        @compileError("At least one of left or right must be cfloat, the other must be an int, float, or cfloat");

    switch (types.numericType(@TypeOf(left))) {
        .int => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return pow(cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}), cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}));
                },
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(right))) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return pow(cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}), cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(@TypeOf(right))) {
                .int => {
                    return pow(cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}), cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}));
                },
                .float => {
                    return pow(cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}), cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{}));
                },
                .cfloat => {
                    const l: Coerce(@TypeOf(left), @TypeOf(right)) = cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{});
                    const r: Coerce(@TypeOf(left), @TypeOf(right)) = cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{});

                    return cfloat.exp(r.mul(cfloat.log(l)));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
