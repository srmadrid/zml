const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Exp(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Exp: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Exp: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Exp: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Exp", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Exp: " ++ @typeName(X) ++ " must implement `fn Exp(type) type`");

            return X.Exp(X);
        },
    }
}

/// Returns the exponential `eˣ` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.exp(x: X) numeric.Exp(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the exponential of.
///
/// ## Returns
/// `numeric.Exp(@TypeOf(x))`: The exponential of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Exp` method. The expected signature and
/// behavior of `Exp` are as follows:
/// * `fn Exp(type) type`: Returns the type of the exponential of `x`.
///
/// `numeric.Exp(X)` or `X` must implement the required `exp` method. The
/// expected signature and behavior of `exp` are as follows:
/// * `fn exp(X) numeric.Exp(X)`: Returns the exponential of `x`.
pub fn exp(x: anytype) numeric.Exp(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Exp(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.exp(x),
        .dyadic => return dyadic.exp(x),
        .complex => return complex.exp(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "exp",
                fn (X) numeric.Exp(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.exp: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn exp(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.exp(x);
        },
    }
}

/// Performs computation of the exponential `eˣ` of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.expInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the exponential of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `expInto` method. The expected
/// signature and behavior of `expInto` are as follows:
/// * `fn expInto(*O, X) void`: Computes the exponential of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `expInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.exp`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn expInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.expInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "expInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.expInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "expInto", fn (*O, X) void, &.{ *O, X }))
                return O.expInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "expInto", fn (*O, X) void, &.{ *O, X }))
            return X.expInto(o, x);
    }

    return numeric.set(o, numeric.exp(x));
}
