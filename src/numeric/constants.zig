const meta = @import("../meta.zig");

const int = @import("../int.zig");
const float = @import("../float.zig");
const dyadic = @import("../dyadic.zig");
const complex = @import("../complex.zig");

/// Returns the maximum representable finite value (highest point on the number
/// line) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the highest value for.
///
/// ## Returns
/// `N`: The maximum representable value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `highest: N` declaration, or implement the `highest`
/// method. The expected signature and behavior of `highest` are as follows:
/// * `fn highest(anytype) N`: Returns the highest representable value.
pub fn highest(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.highest: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return int.highest(N),
        .float => return float.highest(N),
        .dyadic => return .highest,
        .complex => @compileError("zsl.numeric.highest: not defined for " ++ @typeName(N) ++ "."),
        .custom => {
            if (comptime @hasDecl(N, "highest") and @TypeOf(N.highest) == N)
                return N.highest
            else if (comptime meta.hasMethod(N, "highest", fn () N, &.{}))
                return N.highest()
            else
                @compileError("zsl.numeric.highest: " ++ @typeName(N) ++ " must expose a `highest: " ++ @typeName(N) ++ "` declaration, or implement `fn highest() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the minimum representable finite value (lowest point on the number
/// line) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the lowest value for.
///
/// ## Returns
/// `N`: The minimum representable value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `lowest: N` declaration, or implement the `lowest`
/// method. The expected signature and behavior of `lowest` are as follows:
/// * `fn lowest(anytype) N`: Returns the lowest representable value.
pub fn lowest(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.lowest: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return int.lowest(N),
        .float => return float.lowest(N),
        .dyadic => return .lowest,
        .complex => @compileError("zsl.numeric.lowest: not defined for " ++ @typeName(N) ++ "."),
        .custom => {
            if (comptime @hasDecl(N, "lowest") and @TypeOf(N.lowest) == N)
                return N.lowest
            else if (comptime meta.hasMethod(N, "lowest", fn () N, &.{}))
                return N.lowest()
            else
                @compileError("zsl.numeric.lowest: " ++ @typeName(N) ++ " must expose a `lowest: " ++ @typeName(N) ++ "` declaration, or implement `fn lowest() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the smallest positive magnitude strictly greater than zero (closest
/// to zero) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the smallest positive value
///   for.
///
/// ## Returns
/// `N`: The smallest positive non-zero magnitude.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `smallest: N` declaration, or implement the `smallest`
/// method. The expected signature and behavior of `smallest` are as follows:
/// * `fn smallest(anytype) N`: Returns the smallest positive non-zero magnitude.
pub fn smallest(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.smallest: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return int.smallest(N),
        .float => return float.smallest(N),
        .dyadic => return .smallest,
        .complex => @compileError("zsl.numeric.smallest: not defined for " ++ @typeName(N) ++ "."),
        .custom => {
            if (comptime @hasDecl(N, "smallest") and @TypeOf(N.smallest) == N)
                return N.smallest
            else if (comptime meta.hasMethod(N, "smallest", fn () N, &.{}))
                return N.smallest()
            else
                @compileError("zsl.numeric.smallest: " ++ @typeName(N) ++ " must expose a `smallest: " ++ @typeName(N) ++ "` declaration, or implement `fn smallest() " ++ @typeName(N) ++ "`");
        },
    }
}
