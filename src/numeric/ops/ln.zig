const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Ln(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.ln: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.ln: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.ln: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Ln", fn (type) type, &.{X}))
                @compileError("zsl.numeric.ln: " ++ @typeName(X) ++ " must implement `fn Ln(type) type`");

            return X.Ln(X);
        },
    }
}

/// Returns the natural logarithm `logₑ(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.ln(x: X) numeric.Ln(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the natural logarithm of.
///
/// ## Returns
/// `numeric.Ln(@TypeOf(x))`: The natural logarithm of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Ln` method. The expected signature and
/// behavior of `Ln` are as follows:
/// * `fn Ln(type) type`: Returns the type of the natural logarithm of `x`.
///
/// `numeric.Ln(X)` or `X` must implement the required `ln` method. The
/// expected signature and behavior of `ln` are as follows:
/// * `fn ln(X) numeric.Ln(X)`: Returns the natural logarithm of `x`.
pub fn ln(x: anytype) numeric.Ln(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Ln(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.ln(x),
        .dyadic => return dyadic.ln(x),
        .complex => return complex.ln(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "ln",
                fn (X) numeric.Ln(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.ln: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn ln(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.ln(x);
        },
    }
}

/// Performs computation of the natural logarithm `logₑ(x)` of a numeric `x`
/// into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.lnInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the natural logarithm of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `lnInto` method. The expected
/// signature and behavior of `lnInto` are as follows:
/// * `fn lnInto(*O, X) void`: Computes the natural logarithm of `x` and stores
///   it in `o`.
///
/// If neither `O` nor `X` implement the required `lnInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.ln`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn lnInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.lnInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "lnInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.lnInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "lnInto", fn (*O, X) void, &.{ *O, X }))
                return O.lnInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "lnInto", fn (*O, X) void, &.{ *O, X }))
            return X.lnInto(o, x);
    }

    return numeric.set(o, numeric.ln(x));
}
