const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Abs2(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.abs2: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return complex.Abs2(X),
        .custom => {
            if (comptime !meta.hasMethod(X, "Abs2", fn (type) type, &.{X}))
                @compileError("zsl.numeric.abs2: " ++ @typeName(X) ++ " must implement `fn Abs2(type) type`");

            return X.Abs2(X);
        },
    }
}

/// Returns the squared absolute value of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.abs2(x: X) numeric.Abs2(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the squared absolute value of.
///
/// ## Returns
/// `numeric.Abs2(@TypeOf(x))`: The squared absolute value of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Abs2` method. The expected signature and
/// behavior of `Abs2` are as follows:
/// * `fn Abs2(type) type`: Returns the type of the squared absolute value of
///   `x`.
///
/// `numeric.Abs2(X)` or `X` must implement the required `abs2` method. The
/// expected signature and behavior of `abs2` are as follows:
/// * `fn abs2(X) numeric.Abs2(X)`: Returns the squared absolute value of
///   `x`.
pub fn abs2(x: anytype) numeric.Abs2(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Abs2(X);

    switch (comptime meta.numericType(X)) {
        .bool => return x,
        .int => return int.mul(x, x),
        .float => return float.mul(x, x),
        .dyadic => return dyadic.mul(x, x),
        .complex => return complex.abs2(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "abs2",
                fn (X) numeric.Abs2(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.abs2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs2(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.abs2(x);
        },
    }
}

/// Performs computation of the squared absolute value of a numeric `x` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.abs2Into(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the squared absolute value of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `abs2Into` method. The expected
/// signature and behavior of `abs2Into` are as follows:
/// * `fn abs2Into(*O, X) void`: Computes the squared absolute value of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `abs2Into` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.abs2`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn abs2Into(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.abs2Into: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "abs2Into", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.abs2Into(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "abs2Into", fn (*O, X) void, &.{ *O, X }))
                return O.abs2Into(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "abs2Into", fn (*O, X) void, &.{ *O, X }))
            return X.abs2Into(o, x);
    }

    return numeric.set(o, numeric.abs2(x));
}
