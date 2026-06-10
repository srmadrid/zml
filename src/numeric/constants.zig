const meta = @import("../meta.zig");

const int = @import("../int.zig");
const float = @import("../float.zig");
const dyadic = @import("../dyadic.zig");
const complex = @import("../complex.zig");

/// Returns the additive identity (zero) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the zero value for.
///
/// ## Returns
/// `N`: The zero value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the required `zero: N` declaration.
pub fn zero(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.zero: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return false,
        .int => return 0,
        .float => return 0.0,
        .dyadic => return .zero,
        .complex => return .zero,
        .custom => {
            comptime if (!@hasDecl(N, "zero") or @TypeOf(N.zero) != N)
                @compileError("zsl.numeric.zero: " ++ @typeName(N) ++ " must expose a `zero: " ++ @typeName(N) ++ "` declaration");

            return N.zero;
        },
    }
}

/// Returns the multiplicative identity (one) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the one value for.
///
/// ## Returns
/// `N`: The one value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the required `one: N` declaration.
pub fn one(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.one: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return 1,
        .float => return 1.0,
        .dyadic => return .one,
        .complex => return .one,
        .custom => {
            comptime if (!@hasDecl(N, "one") or @TypeOf(N.one) != N)
                @compileError("zsl.numeric.one: " ++ @typeName(N) ++ " must expose a `one: " ++ @typeName(N) ++ "` declaration");

            return N.one;
        },
    }
}

/// Returns the numeric constant two for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the two value for.
///
/// ## Returns
/// `N`: The two value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the required `two: N` declaration.
pub fn two(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.two: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return 2,
        .float => return 2.0,
        .dyadic => return .two,
        .complex => return .two,
        .custom => {
            comptime if (!@hasDecl(N, "two") or @TypeOf(N.two) != N)
                @compileError("zsl.numeric.two: " ++ @typeName(N) ++ " must expose a `two: " ++ @typeName(N) ++ "` declaration");

            return N.two;
        },
    }
}
