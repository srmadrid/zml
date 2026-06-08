const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Sqrt(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Sqrt: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Sqrt: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Sqrt: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Sqrt", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Sqrt: " ++ @typeName(X) ++ " must implement `fn Sqrt(type) type`");

            return X.Sqrt(X);
        },
    }
}

/// Returns the square root `√x` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.sqrt(x: X) numeric.Sqrt(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the square root of.
///
/// ## Returns
/// `numeric.Sqrt(@TypeOf(x))`: The square root of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Sqrt` method. The expected signature and
/// behavior of `Sqrt` are as follows:
/// * `fn Sqrt(type) type`: Returns the type of the square root of `x`.
///
/// `numeric.Sqrt(X)` or `X` must implement the required `sqrt` method. The
/// expected signature and behavior of `sqrt` are as follows:
/// * `fn sqrt(X) numeric.Sqrt(X)`: Returns the square root of `x`.
pub fn sqrt(x: anytype) numeric.Sqrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Sqrt(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.sqrt(x),
        .dyadic => return dyadic.sqrt(x),
        .complex => return complex.sqrt(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "sqrt",
                fn (X) numeric.Sqrt(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.sqrt: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sqrt(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.sqrt(x);
        },
    }
}

/// Performs computation of the square root `√x` of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.sqrtInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the square root of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `sqrtInto` method. The expected
/// signature and behavior of `sqrtInto` are as follows:
/// * `fn sqrtInto(*O, X) void`: Computes the square root of `x` and stores it
///   in `o`.
///
/// If neither `O` nor `X` implement the required `sqrtInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.sqrt`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn sqrtInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.sqrtInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "sqrtInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.sqrtInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "sqrtInto", fn (*O, X) void, &.{ *O, X }))
                return O.sqrtInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "sqrtInto", fn (*O, X) void, &.{ *O, X }))
            return X.sqrtInto(o, x);
    }

    return numeric.set(o, numeric.sqrt(x));
}
