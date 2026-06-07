const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Cos(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.cos: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.cos: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.cos: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Cos", fn (type) type, &.{X}))
                @compileError("zsl.numeric.cos: " ++ @typeName(X) ++ " must implement `fn Cos(type) type`");

            return X.Cos(X);
        },
    }
}

/// Returns the cosine `cos(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.cos(x: X) numeric.Cos(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the cosine of.
///
/// ## Returns
/// `numeric.Cos(@TypeOf(x))`: The cosine of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Cos` method. The expected signature and
/// behavior of `Cos` are as follows:
/// * `fn Cos(type) type`: Returns the type of the cosine of `x`.
///
/// `numeric.Cos(X)` or `X` must implement the required `cos` method. The
/// expected signature and behavior of `cos` are as follows:
/// * `fn Cos(X) numeric.Cos(X)`: Returns the cosine of `x`.
pub fn cos(x: anytype) numeric.Cos(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Cos(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.cos(x),
        .dyadic => return dyadic.cos(x),
        .complex => return complex.cos(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "cos",
                fn (X) numeric.Cos(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.cos: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn cos(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.cos(x);
        },
    }
}

/// Performs computation of the cosine `cos(x)` of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.cosInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the cosine of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `cosInto` method. The expected
/// signature and behavior of `cosInto` are as follows:
/// * `fn cosInto(*O, X) void`: Computes the cosine of `x` and stores it in `o`.
///
/// If neither `O` nor `X` implement the required `cosInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.cos`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn cosInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.cosInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "cosInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.cosInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "cosInto", fn (*O, X) void, &.{ *O, X }))
                return O.cosInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "cosInto", fn (*O, X) void, &.{ *O, X }))
            return X.cosInto(o, x);
    }

    return numeric.set(o, numeric.cos(x));
}
