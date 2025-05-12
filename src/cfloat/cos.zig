const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;
const Cfloat = @import("../cfloat.zig").Cfloat;
const cast = types.cast;

pub fn cos(z: anytype) Cfloat(EnsureFloat(Scalar(@TypeOf(z)))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("z must be an int, float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return cos(cast(Cfloat(EnsureFloat(Scalar(@TypeOf(z)))), z, .{}));
        },
        .float => {
            return cos(cast(Cfloat(@TypeOf(z)), z, .{}));
        },
        .cfloat => {
            const y: @TypeOf(z) = .{
                .re = -z.im,
                .im = z.re,
            };

            return cfloat.cosh(y);
        },
        else => unreachable,
    }
}
