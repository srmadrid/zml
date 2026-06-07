const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Atan2(Y: type, X: type) type {
    comptime if (!meta.isNumeric(Y) or !meta.isNumeric(X))
        @compileError("zsl.numeric.atan2: x and y must be numerics, got \n\ty: " ++ @typeName(Y) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    if (comptime meta.isCustomNumeric(Y)) {
        if (comptime meta.isCustomNumeric(X)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ Y, X },
                "Atan2",
                fn (type, type) type,
                &.{ Y, X },
            ) orelse
                @compileError("zsl.numeric.atan2: " ++ @typeName(Y) ++ " or " ++ @typeName(X) ++ " must implement `fn Atan2(type, type) type`");

            return Impl.Atan2(Y, X);
        } else { // only Y custom
            comptime if (!meta.hasMethod(Y, "Atan2", fn (type, type) type, &.{ Y, X }))
                @compileError("zsl.numeric.atan2: " ++ @typeName(Y) ++ " must implement `fn Atan2(type, type) type`");

            return Y.Atan2(Y, X);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only Y custom
        comptime if (!meta.hasMethod(Y, "Atan2", fn (type, type) type, &.{ Y, X }))
            @compileError("zsl.numeric.atan2: " ++ @typeName(X) ++ " must implement `fn Atan2(type, type) type`");

        return X.Atan2(Y, X);
    }

    switch (comptime meta.numericType(Y)) {
        .bool => switch (comptime meta.numericType(X)) {
            .bool => @compileError("zsl.numeric.atan2: not defined for " ++ @typeName(Y) ++ " and " ++ @typeName(X) ++ "."),
            .int => @compileError("zsl.numeric.atan2: not defined for " ++ @typeName(Y) ++ " and " ++ @typeName(X) ++ "."),
            .float => return float.Atan2(Y, X),
            .dyadic => return dyadic.Atan2(Y, X),
            .complex => return complex.Atan2(Y, X),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(X)) {
            .bool, .int => @compileError("zsl.numeric.atan2: not defined for " ++ @typeName(Y) ++ " and " ++ @typeName(X) ++ "."),
            .float => return float.Atan2(Y, X),
            .dyadic => return dyadic.Atan2(Y, X),
            .complex => return complex.Atan2(Y, X),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(X)) {
            .bool, .int, .float => return float.Atan2(Y, X),
            .dyadic => return dyadic.Atan2(Y, X),
            .complex => return complex.Atan2(Y, X),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(X)) {
            .bool, .int, .float, .dyadic => return dyadic.Atan2(Y, X),
            .complex => return complex.Atan2(Y, X),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(X)) {
            .bool, .int, .float, .dyadic, .complex => return complex.Atan2(Y, X),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Computes the arctangent `tan⁻¹(y/x)` of any two numeric operands, using the
/// signs of both arguments to determine the correct quadrant of the result.
///
/// ## Signature
/// ```zig
/// numeric.atan2(y: Y, x: X) numeric.Atan2(Y, X)
/// ```
///
/// ## Arguments
/// * `y` (`anytype`): The `y` coordinate.
/// * `x` (`anytype`): The `x` coordinate.
///
/// ## Returns
/// `numeric.Atan2(@TypeOf(y), @TypeOf(x))`: The arctangent of `y/x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `Y` or `X` must implement the required `Atan2` method. The expected
/// signature and behavior of `Atan2` are as follows:
/// * `fn Atan2(type, type) type`: Returns the type of the arctangent of `y/x`.
///
/// `numeric.Atan2(Y, X)`, `Y` or `X` must implement the required `atan2` method.
/// The expected signatures and behavior of `atan2` are as follows:
/// * `fn atan2(Y, X) numeric.Atan2(Y, X)`: Returns the arctangent of `y/x`.
pub fn atan2(y: anytype, x: anytype) numeric.Atan2(@TypeOf(y), @TypeOf(x)) {
    const Y: type = @TypeOf(y);
    const X: type = @TypeOf(x);
    const R: type = numeric.Atan2(Y, X);

    if (comptime meta.isCustomNumeric(Y)) {
        if (comptime meta.isCustomNumeric(X)) { // X and Y both custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, Y, X },
                "atan2",
                fn (Y, X) R,
                &.{ Y, X },
            ) orelse
                @compileError("zsl.numeric.atan2: " ++ @typeName(R) ++ ", " ++ @typeName(Y) ++ " or " ++ @typeName(X) ++ " must implement `fn atan2(" ++ @typeName(Y) ++ ", " ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.atan2(y, x);
        } else { // only Y custom
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, Y },
                "atan2",
                fn (Y, X) R,
                &.{ Y, X },
            ) orelse
                @compileError("zsl.numeric.atan2: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn atan2(" ++ @typeName(Y) ++ ", " ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.atan2(y, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        const Impl: type = comptime meta.anyHasMethod(
            &.{ R, X },
            "atan2",
            fn (Y, X) R,
            &.{ Y, X },
        ) orelse
            @compileError("zsl.numeric.atan2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn atan2(" ++ @typeName(Y) ++ ", " ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

        return Impl.atan2(y, x);
    }

    switch (comptime meta.numericType(Y)) {
        .bool => switch (comptime meta.numericType(X)) {
            .bool => unreachable,
            .int => unreachable,
            .float => return float.atan2(y, x),
            .dyadic => return dyadic.atan2(y, x),
            .complex => return complex.atan2(y, x),
            .custom => unreachable,
        },
        .int => switch (comptime meta.numericType(X)) {
            .bool, .int => unreachable,
            .float => return float.atan2(y, x),
            .dyadic => return dyadic.atan2(y, x),
            .complex => return complex.atan2(y, x),
            .custom => unreachable,
        },
        .float => switch (comptime meta.numericType(X)) {
            .bool, .int, .float => return float.atan2(y, x),
            .dyadic => return dyadic.atan2(y, x),
            .complex => return complex.atan2(y, x),
            .custom => unreachable,
        },
        .dyadic => switch (comptime meta.numericType(X)) {
            .bool, .int, .float, .dyadic => return dyadic.atan2(y, x),
            .complex => return complex.atan2(y, x),
            .custom => unreachable,
        },
        .complex => switch (comptime meta.numericType(X)) {
            .bool, .int, .float, .dyadic, .complex => return complex.atan2(y, x),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs computation of the arctangent `tan⁻¹(y/x)` of two numerics `y` and
/// `x` into a numeric `o`, using the signs of both arguments to determine the
/// correct quadrant of the result.
///
/// ## Signature
/// ```zig
/// numeric.atan2Into(o: *O, y: Y, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `y` (`anytype`): The `y` coordinate.
/// * `x` (`anytype`): The `x` coordinate.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O`, `Y` or `X` should implement the required `atan2Into` method. The
/// expected signature and behavior of `atan2Into` are as follows:
/// * `fn atan2Into(*O, Y, X) void`: Computes the arctangent of `y/x` and stores
///   it in `o`.
///
/// If none of `O`, `Y` and `X` implement the required `atan2Into` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.atan2`, potentially resulting in a less efficient implementation. In
/// this case, `O`, `Y` and `X` must adhere to the requirements of these
/// functions.
pub fn atan2Into(o: anytype, y: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const Y: type = @TypeOf(y);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(Y) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.atan2Into: o must be a mutable one-item pointer to a numeric, and y and x must be numerics, got \n\to: " ++ @typeName(O) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(Y)) {
            if (comptime meta.isCustomNumeric(X)) { // O, Y and X all custom
                if (comptime meta.anyHasMethod(&.{ O, Y, X }, "atan2Into", fn (*O, Y, X) void, &.{ *O, Y, X })) |Impl|
                    return Impl.atan2Into(o, y, x);
            } else { // only O and Y custom
                if (comptime meta.anyHasMethod(&.{ O, Y }, "atan2Into", fn (*O, Y, X) void, &.{ *O, Y, X })) |Impl|
                    return Impl.atan2Into(o, y, x);
            }
        } else {
            if (comptime meta.isCustomNumeric(X)) { // only O and X custom
                if (comptime meta.anyHasMethod(&.{ O, X }, "atan2Into", fn (*O, Y, X) void, &.{ *O, Y, X })) |Impl|
                    return Impl.atan2Into(o, y, x);
            } else { // only O custom
                if (comptime meta.hasMethod(O, "atan2Into", fn (*O, Y, X) void, &.{ *O, Y, X }))
                    return O.atan2Into(o, y, x);
            }
        }
    } else {
        if (comptime meta.isCustomNumeric(Y)) {
            if (comptime meta.isCustomNumeric(X)) { // only Y and X custom
                if (comptime meta.anyHasMethod(&.{ Y, X }, "atan2Into", fn (*O, Y, X) void, &.{ *O, Y, X })) |Impl|
                    return Impl.atan2Into(o, y, x);
            } else { // only Y custom
                if (comptime meta.hasMethod(Y, "atan2Into", fn (*O, Y, X) void, &.{ *O, Y, X }))
                    return Y.atan2Into(o, y, x);
            }
        } else if (comptime meta.isCustomNumeric(X)) { // only X custom
            if (comptime meta.hasMethod(Y, "atan2Into", fn (*O, Y, X) void, &.{ *O, Y, X }))
                return Y.atan2Into(o, y, x);
        }
    }

    return numeric.set(o, numeric.atan2(y, x));
}
