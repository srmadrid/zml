const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Atanh(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.atanh: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.atanh: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.atanh: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Atanh", fn (type) type, &.{X}))
                @compileError("zsl.numeric.atanh: " ++ @typeName(X) ++ " must implement `fn Atanh(type) type`");

            return X.Atanh(X);
        },
    }
}

/// Returns the hyperbolic arctangent `tanh⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.atanh(x: X) numeric.Atanh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic arctangent of.
///
/// ## Returns
/// `numeric.Atanh(@TypeOf(x))`: The hyperbolic arctangent of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Atanh` method. The expected signature
/// and behavior of `Atanh` are as follows:
/// * `fn Atanh(type) type`: Returns the type of the hyperbolic arctangent of
///   `x`.
///
/// `numeric.Atanh(X)` or `X` must implement the required `atanh` method. The
/// expected signature and behavior of `atanh` are as follows:
/// * `fn atanh(X) numeric.Atanh(X)`: Returns the hyperbolic arctangent of `x`.
pub fn atanh(x: anytype) numeric.Atanh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Atanh(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.atanh(x),
        .dyadic => return dyadic.atanh(x),
        .complex => return complex.atanh(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "atanh",
                fn (X) numeric.Atanh(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.atanh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn atanh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.atanh(x);
        },
    }
}

/// Performs computation of the hyperbolic arctangent `tanh⁻¹(x)` of a numeric
/// `x` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.atanhInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the hyperbolic arctangent of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `atanhInto` method. The expected
/// signature and behavior of `atanhInto` are as follows:
/// * `fn atanhInto(*O, X) void`: Computes the hyperbolic arctangent of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `atanhInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.atanh`, potentially resulting in a less efficient implementation.
/// In this case, `O` and `X` must adhere to the requirements of these functions.
pub fn atanhInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.atanhInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "atanhInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.atanhInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "atanhInto", fn (*O, X) void, &.{ *O, X }))
                return O.atanhInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "atanhInto", fn (*O, X) void, &.{ *O, X }))
            return X.atanhInto(o, x);
    }

    return numeric.set(o, numeric.atanh(x));
}
