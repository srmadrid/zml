const types = @import("../types.zig");
const float = @import("../float.zig");
const Scalar = types.Scalar;
const EnsureFloat = types.EnsureFloat;

pub fn arg(z: anytype) EnsureFloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)))
        @compileError("cfloat.arg: z must be a bool, int, float or cfloat, got " ++ @typeName(@TypeOf(z)));

    switch (types.numericType(@TypeOf(z))) {
        .int => {
            return if (z >= 0) 0 else float.pi(EnsureFloat(@TypeOf(z)));
        },
        .float => {
            return if (z >= 0) 0 else float.pi(@TypeOf(z));
        },
        .cfloat => {
            return float.atan2(z.im, z.re);
        },
        else => unreachable,
    }
}
