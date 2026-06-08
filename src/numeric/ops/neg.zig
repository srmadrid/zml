const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Neg(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Neg: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Neg", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Neg: " ++ @typeName(X) ++ " must implement `fn Neg(type) type`");

            return X.Neg(X);
        },
    }
}

/// Returns the negation of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.neg(x: X) numeric.Neg(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the negation of.
///
/// ## Returns
/// `numeric.Neg(@TypeOf(x))`: The negation of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Neg` method. The expected signature and
/// behavior of `Neg` are as follows:
/// * `fn Neg(type) type`: Returns the type of the negation of `x`.
///
/// `numeric.Neg(X)` or `X` must implement the required `neg` method. The
/// expected signature and behavior of `neg` are as follows:
/// * `fn neg(X) numeric.Neg(X)`: Returns the negation of `x`.
pub fn neg(x: anytype) numeric.Neg(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Neg(X);

    switch (comptime meta.numericType(X)) {
        .bool => return !x,
        .int => return -x,
        .float => return -x,
        .dyadic => return x.neg(),
        .complex => return x.neg(),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "neg",
                fn (X) numeric.Neg(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.neg: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn neg(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.neg(x);
        },
    }
}

/// Performs computation of the negation of a numeric `x` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.negInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the negation of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `negInto` method. The expected
/// signature and behavior of `negInto` are as follows:
/// * `fn negInto(*O, X) void`: Computes the negation of `x` and stores it in `o`.
///
/// If neither `O` nor `X` implement the required `negInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.neg`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn negInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.negInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "negInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.negInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "negInto", fn (*O, X) void, &.{ *O, X }))
                return O.negInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "negInto", fn (*O, X) void, &.{ *O, X }))
            return X.negInto(o, x);
    }

    return numeric.set(o, numeric.neg(x));
}
