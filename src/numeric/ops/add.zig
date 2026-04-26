const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Add(X: type, Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y))
        @compileError("zsl.numeric.add: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ X, Y },
                "Add",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.add: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Add(type, type) type`");

            return Impl.Add(X, Y);
        } else { // only X custom
            comptime if (!meta.hasMethod(X, "Add", fn (type, type) type, &.{ X, Y }))
                @compileError("zsl.numeric.add: " ++ @typeName(X) ++ " must implement `fn Add(type, type) type`");

            return X.Add(X, Y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Add", fn (type, type) type, &.{ X, Y }))
            @compileError("zsl.numeric.add: " ++ @typeName(Y) ++ " must implement `fn Add(type, type) type`");

        return Y.Add(X, Y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => @compileError("zsl.numeric.add: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int => return int.Add(X, Y),
            .float => return float.Add(X, Y),
            .dyadic => return dyadic.Add(X, Y),
            .complex => return complex.Add(X, Y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.Add(X, Y),
            .float => return float.Add(X, Y),
            .dyadic => return dyadic.Add(X, Y),
            .complex => return complex.Add(X, Y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.Add(X, Y),
            .dyadic => return dyadic.Add(X, Y),
            .complex => return complex.Add(X, Y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Add(X, Y),
            .complex => return complex.Add(X, Y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.Add(X, Y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs addition between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.add(x: X, y: Y) numeric.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `numeric.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Add` method. The expected
/// signature and behavior of `Add` are as follows:
/// * `fn Add(type, type) type`: Returns the type of `x + y`.
///
/// `numeric.Add(X, Y)`, `X` or `Y` must implement the required `add` method.
/// The expected signatures and behavior of `add` are as follows:
/// * `fn add(X, Y) numeric.Add(X, Y)`: Returns the addition of `x` and `y`.
pub fn add(x: anytype, y: anytype) numeric.Add(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Add(X, Y);

    if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X, Y },
                "add",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.add: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn add(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.add(x, y);
        } else { // only X custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "add",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zsl.numeric.add: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn add(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.add(x, y);
        }
    } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, Y },
            "add",
            fn (X, Y) R,
            &.{ X, Y },
        ) orelse
            @compileError("zsl.numeric.add: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn add(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.add(x, y);
    }

    switch (comptime meta.numericType(X)) {
        .bool => switch (comptime meta.numericType(Y)) {
            .bool => unreachable,
            .int => return int.add(x, y),
            .float => return float.add(x, y),
            .dyadic => return dyadic.add(x, y),
            .complex => return complex.add(x, y),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(Y)) {
            .bool, .int => return int.add(x, y),
            .float => return float.add(x, y),
            .dyadic => return dyadic.add(x, y),
            .complex => return complex.add(x, y),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float => return float.add(x, y),
            .dyadic => return dyadic.add(x, y),
            .complex => return complex.add(x, y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.add(x, y),
            .complex => return complex.add(x, y),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic, .complex => return complex.add(x, y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}
