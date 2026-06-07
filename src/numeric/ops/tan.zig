const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Tan(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.tan: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.tan: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.tan: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Tan", fn (type) type, &.{X}))
                @compileError("zsl.numeric.tan: " ++ @typeName(X) ++ " must implement `fn Tan(type) type`");

            return X.Tan(X);
        },
    }
}

/// Returns the tangent `tan(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.tan(x: X) numeric.Tan(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the tangent of.
///
/// ## Returns
/// `numeric.Tan(@TypeOf(x))`: The tangent of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Tan` method. The expected signature and
/// behavior of `Tan` are as follows:
/// * `fn Tan(type) type`: Returns the type of the tangent of `x`.
///
/// `numeric.Tan(X)` or `X` must implement the required `tan` method. The
/// expected signature and behavior of `tan` are as follows:
/// * `fn tan(X) numeric.Tan(X)`: Returns the tangent of `x`.
pub fn tan(x: anytype) numeric.Tan(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Tan(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.tan(x),
        .dyadic => return dyadic.tan(x),
        .complex => return complex.tan(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "tan",
                fn (X) numeric.Tan(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.tan: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn tan(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.tan(x);
        },
    }
}

/// Performs computation of the tangent `tan(x)` of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.tanInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the tangent of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `tanInto` method. The expected
/// signature and behavior of `tanInto` are as follows:
/// * `fn tanInto(*O, X) void`: Computes the tangent of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `tanInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.tan`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn tanInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.tanInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "tanInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.tanInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "tanInto", fn (*O, X) void, &.{ *O, X }))
                return O.tanInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "tanInto", fn (*O, X) void, &.{ *O, X }))
            return X.tanInto(o, x);
    }

    return numeric.set(o, numeric.tan(x));
}
