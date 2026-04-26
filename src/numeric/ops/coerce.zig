const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

/// Coerces any two numeric types to the smallest type that can represent all
/// values representable by either type.
///
/// For two ints, if they have different signedness, the result is a signed int.
/// The bit-width of the result is either the larger of the two bit-widths (if
/// the signed type is larger) or the larger of the two bit-widths plus one (if
/// the unsigned type is larger). If both ints are "standard" (see
/// `meta.standard_integer_types`), the result is the next larger standard type
/// that can hold both values.
///
/// ## Arguments
/// * `X` (`comptime type`): The first type to coerce. Must be a numeric type.
/// * `Y` (`comptime type`): The second type to coerce. Must be a numeric type.
///
/// ## Returns
/// `type`: The coerced type that can represent all values of both `X` and `Y`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Coerce` method. The expected
/// signature and behavior of `Coerce` are as follows:
/// * `fn Coerce(type, type) type`: Returns the smallest type that can represent
///   all values of types `X` and `Y`.
pub fn Coerce(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.Coerce: X and Y must be numeric types, got \n\tX = " ++ @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "Coerce",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.Coerce: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Coerce(type, type) type`");

            return Impl.Coerce(X, Y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "Coerce", fn (type, type) type, &.{ X, Y }))
                @compileError("zsl.numeric.Coerce: " ++ @typeName(X) ++ " must implement `fn Coerce(type, type) type`");

            return X.Coerce(X, Y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Coerce", fn (type, type) type, &.{ X, Y }))
            @compileError("zsl.numeric.Coerce: " ++ @typeName(Y) ++ " must implement `fn Coerce(type, type) type`");

        return Y.Coerce(X, Y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => return bool,
            .int => return int.Coerce(X, Y),
            .float => return float.Coerce(X, Y),
            .dyadic => return dyadic.Coerce(X, Y),
            .complex => return complex.Coerce(X, Y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.Coerce(X, Y),
            .float => return float.Coerce(X, Y),
            .dyadic => return dyadic.Coerce(X, Y),
            .complex => return complex.Coerce(X, Y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.Coerce(X, Y),
            .dyadic => return dyadic.Coerce(X, Y),
            .complex => return complex.Coerce(X, Y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Coerce(X, Y),
            .complex => return complex.Coerce(X, Y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.Coerce(X, Y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}
