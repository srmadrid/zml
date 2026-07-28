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
/// `N` must expose the `zero: N` declaration, or implement the `zero` method.
/// The expected signature and behavior of `zero` are as follows:
/// * `fn zero(anytype) N`: Returns the zero value.
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
            if (comptime @hasDecl(N, "zero") and @TypeOf(N.zero) == N)
                return N.zero
            else if (comptime meta.hasMethod(N, "zero", fn () N, &.{}))
                return N.zero()
            else
                @compileError("zsl.numeric.zero: " ++ @typeName(N) ++ " must expose a `zero: " ++ @typeName(N) ++ "` declaration, or implement `fn zero() " ++ @typeName(N) ++ "`");
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
/// `N` must expose the `one: N` declaration, or implement the `one` method.
/// The expected signature and behavior of `one` are as follows:
/// * `fn one(anytype) N`: Returns the one value.
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
            if (comptime @hasDecl(N, "one") and @TypeOf(N.one) == N)
                return N.one
            else if (comptime meta.hasMethod(N, "one", fn () N, &.{}))
                return N.one()
            else
                @compileError("zsl.numeric.one: " ++ @typeName(N) ++ " must expose a `one: " ++ @typeName(N) ++ "` declaration, or implement `fn one() " ++ @typeName(N) ++ "`");
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
/// `N` must expose the `two: N` declaration, or implement the `two` method.
/// The expected signature and behavior of `two` are as follows:
/// * `fn two(anytype) N`: Returns the two value.
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
            if (comptime @hasDecl(N, "two") and @TypeOf(N.two) == N)
                return N.two
            else if (comptime meta.hasMethod(N, "two", fn () N, &.{}))
                return N.two()
            else
                @compileError("zsl.numeric.two: " ++ @typeName(N) ++ " must expose a `two: " ++ @typeName(N) ++ "` declaration, or implement `fn two() " ++ @typeName(N) ++ "`");
        },
    }
}
