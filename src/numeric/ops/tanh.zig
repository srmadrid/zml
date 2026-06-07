const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Tanh(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.tanh: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.tanh: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.tanh: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Tanh", fn (type) type, &.{X}))
                @compileError("zsl.numeric.tanh: " ++ @typeName(X) ++ " must implement `fn Tanh(type) type`");

            return X.Tanh(X);
        },
    }
}

/// Returns the hyperbolic tangent `tanh(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.tanh(x: X) numeric.Tanh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic tangent of.
///
/// ## Returns
/// `numeric.Tanh(@TypeOf(x))`: The hyperbolic tangent of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Tanh` method. The expected signature and
/// behavior of `Tanh` are as follows:
/// * `fn Tanh(type) type`: Returns the type of the hyperbolic tangent of `x`.
///
/// `numeric.Tanh(X)` or `X` must implement the required `tanh` method. The
/// expected signature and behavior of `tanh` are as follows:
/// * `fn tanh(X) numeric.Tanh(X)`: Returns the hyperbolic tangent of `x`.
pub fn tanh(x: anytype) numeric.Tanh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Tanh(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.tanh(x),
        .dyadic => return dyadic.tanh(x),
        .complex => return complex.tanh(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "tanh",
                fn (X) numeric.Tanh(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.tanh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn tanh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.tanh(x);
        },
    }
}

/// Performs computation of the hyperbolic tangent `tanh(x)` of a numeric `x`
/// into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.tanhInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the hyperbolic tangent of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `tanhInto` method. The expected
/// signature and behavior of `tanhInto` are as follows:
/// * `fn tanhInto(*O, X) void`: Computes the hyperbolic tangent of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `tanhInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.tanh`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn tanhInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.tanhInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "tanhInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.tanhInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "tanhInto", fn (*O, X) void, &.{ *O, X }))
                return O.tanhInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "tanhInto", fn (*O, X) void, &.{ *O, X }))
            return X.tanhInto(o, x);
    }

    return numeric.set(o, numeric.tanh(x));
}
