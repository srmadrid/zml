const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Abs(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Abs: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return complex.Abs(X),
        .custom => {
            if (comptime !meta.hasMethod(X, "Abs", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Abs: " ++ @typeName(X) ++ " must implement `fn Abs(type) type`");

            return X.Abs(X);
        },
    }
}

/// Returns the absolute value of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.abs(x: X) numeric.Abs(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the absolute value of.
///
/// ## Returns
/// `numeric.Abs(@TypeOf(x))`: The absolute value of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `abs` method. The expected signature and
/// behavior of `abs` are as follows:
/// * `fn abs(type) type`: Returns the type of the absolute value of `x`.
///
/// `numeric.Abs(X)` or `X` must implement the required `abs` method. The
/// expected signature and behavior of `abs` are as follows:
/// * `fn abs(X) numeric.Abs(X)`: Returns the absolute value of `x`.
pub fn abs(x: anytype) numeric.Abs(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Abs(X);

    switch (comptime meta.numericType(X)) {
        .bool => return x,
        .int => return int.abs(x),
        .float => return float.abs(x),
        .dyadic => return dyadic.abs(x),
        .complex => return complex.abs(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "abs",
                fn (X) numeric.Abs(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.abs: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.abs(x);
        },
    }
}

/// Performs computation of the absolute value of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.absInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the absolute value of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `absInto` method. The expected
/// signature and behavior of `absInto` are as follows:
/// * `fn absInto(*O, X) void`: Computes the absolute value of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `absInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.abs`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn absInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.absInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "absInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.absInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "absInto", fn (*O, X) void, &.{ *O, X }))
                return O.absInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "absInto", fn (*O, X) void, &.{ *O, X }))
            return X.absInto(o, x);
    }

    return numeric.set(o, numeric.abs(x));
}
