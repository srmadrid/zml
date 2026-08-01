const meta = @import("../../meta.zig");

const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

pub fn isNan(x: anytype) bool {
    const X: type = @TypeOf(x);

    switch (comptime meta.numericType(X)) {
        .bool => return false,
        .int => return false,
        .float => return float.isNan(x),
        .dyadic => return dyadic.isNan(x),
        .complex => return complex.isNan(),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{X},
                "isNan",
                fn (X) bool,
                &.{X},
            ) orelse
                @compileError("zsl.numeric.isNan: " ++ @typeName(X) ++ " must implement `fn isNan(" ++ @typeName(X) ++ ") " ++ "bool`");

            return Impl.isNan(x);
        },
    }
}
