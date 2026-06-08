const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Abs1(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Abs1: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return complex.Abs1(X),
        .custom => {
            if (comptime !meta.hasMethod(X, "Abs1", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Abs1: " ++ @typeName(X) ++ " must implement `fn Abs1(type) type`");

            return X.Abs1(X);
        },
    }
}

/// Returns the 1-norm of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.abs1(x: X) numeric.Abs1(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the 1-norm of.
///
/// ## Returns
/// `numeric.Abs1(@TypeOf(x))`: The 1-norm of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Abs1` method. The expected signature and
/// behavior of `Abs1` are as follows:
/// * `fn Abs1(type) type`: Returns the type of the 1-norm of `x`.
///
/// `numeric.Abs1(X)` or `X` must implement the required `abs1` method. The
/// expected signature and behavior of `abs1` are as follows:
/// * `fn abs1(X) numeric.Abs1(X)`: Returns the 1-norm of `x`.
pub fn abs1(x: anytype) numeric.Abs1(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Abs1(X);

    switch (comptime meta.numericType(X)) {
        .bool => return x,
        .int => return int.abs(x),
        .float => return float.abs(x),
        .dyadic => return dyadic.abs(x),
        .complex => return complex.abs1(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "abs1",
                fn (X) numeric.Abs1(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.abs1: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs1(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.abs1(x);
        },
    }
}

/// Performs computation of the 1-norm of a numeric `x` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.abs1Into(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the 1-norm of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `abs1Into` method. The expected
/// signature and behavior of `abs1Into` are as follows:
/// * `fn abs1Into(*O, X) void`: Computes the 1-norm of `x` and stores it in `o`.
///
/// If neither `O` nor `X` implement the required `abs1Into` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.abs1`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn abs1Into(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.abs1Into: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "abs1Into", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.abs1Into(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "abs1Into", fn (*O, X) void, &.{ *O, X }))
                return O.abs1Into(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "abs1Into", fn (*O, X) void, &.{ *O, X }))
            return X.abs1Into(o, x);
    }

    return numeric.set(o, numeric.abs1(x));
}
