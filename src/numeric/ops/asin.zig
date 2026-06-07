const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Asin(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.asin: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.asin: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.asin: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Asin", fn (type) type, &.{X}))
                @compileError("zsl.numeric.asin: " ++ @typeName(X) ++ " must implement `fn Asin(type) type`");

            return X.Asin(X);
        },
    }
}

/// Returns the arcsine `sin⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.asin(x: X) numeric.Asin(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the arcsine of.
///
/// ## Returns
/// `numeric.Asin(@TypeOf(x))`: The arcsine of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Asin` method. The expected signature and
/// behavior of `Asin` are as follows:
/// * `fn Asin(type) type`: Returns the type of the arcsine of `x`.
///
/// `numeric.Asin(X)` or `X` must implement the required `asin` method. The
/// expected signature and behavior of `asin` are as follows:
/// * `fn asin(X) numeric.Asin(X)`: Returns the arcsine of `x`.
pub fn asin(x: anytype) numeric.Asin(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Asin(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.asin(x),
        .dyadic => return dyadic.asin(x),
        .complex => return complex.asin(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "asin",
                fn (X) numeric.Asin(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.asin: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn asin(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.asin(x);
        },
    }
}

/// Performs computation of the arcsine `sin⁻¹(x)` of a numeric `x` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.asinInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the arcsine of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `asinInto` method. The expected
/// signature and behavior of `asinInto` are as follows:
/// * `fn asinInto(*O, X) void`: Computes the arcsine of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `asinInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.asin`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn asinInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.asinInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "asinInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.asinInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "asinInto", fn (*O, X) void, &.{ *O, X }))
                return O.asinInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "asinInto", fn (*O, X) void, &.{ *O, X }))
            return X.asinInto(o, x);
    }

    return numeric.set(o, numeric.asin(x));
}
