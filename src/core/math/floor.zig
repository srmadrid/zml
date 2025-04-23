const types = @import("../types.zig");
const EnsureFloat = types.EnsureFloat;
const cast = @import("../types.zig").cast;

pub fn floor(x: anytype) EnsureFloat(@TypeOf(x)) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            return @floor(cast(EnsureFloat(@TypeOf(x)), x, .{}));
        },
        .float => {
            switch (@TypeOf(x)) {
                f16 => return @floor(x),
                f32 => return @floor(x),
                f64 => return @floor(x),
                f80 => return @floor(x),
                f128 => return @floor(x),
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
