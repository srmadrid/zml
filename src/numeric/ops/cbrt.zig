const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Cbrt(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.cbrt: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.cbrt: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.cbrt: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Cbrt", fn (type) type, &.{X}))
                @compileError("zsl.numeric.cbrt: " ++ @typeName(X) ++ " must implement `fn Cbrt(type) type`");

            return X.Cbrt(X);
        },
    }
}

/// Returns the cube root `∛x` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.cbrt(x: X) numeric.Cbrt(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the cube root of.
///
/// ## Returns
/// `numeric.Cbrt(@TypeOf(x))`: The cube root of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Cbrt` method. The expected signature and
/// behavior of `Cbrt` are as follows:
/// * `fn Cbrt(type) type`: Returns the type of the cube root of `x`.
///
/// `numeric.Cbrt(X)` or `X` must implement the required `cbrt` method. The
/// expected signature and behavior of `cbrt` are as follows:
/// * `fn cbrt(X) numeric.Cbrt(X)`: Returns the cube root of `x`.
pub fn cbrt(x: anytype) numeric.Cbrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Cbrt(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.cbrt(x),
        .dyadic => return dyadic.cbrt(x),
        .complex => return complex.cbrt(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "cbrt",
                fn (X) numeric.Cbrt(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.cbrt: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn cbrt(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.cbrt(x);
        },
    }
}

/// Performs computation of the cube root `∛x` of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.cbrtInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the cube root of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `cbrtInto` method. The expected
/// signature and behavior of `cbrtInto` are as follows:
/// * `fn cbrtInto(*O, X) void`: Computes the cube root of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `cbrtInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.cbrt`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn cbrtInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.cbrtInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "cbrtInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.cbrtInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "cbrtInto", fn (*O, X) void, &.{ *O, X }))
                return O.cbrtInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "cbrtInto", fn (*O, X) void, &.{ *O, X }))
            return X.cbrtInto(o, x);
    }

    return numeric.set(o, numeric.cbrt(x));
}
