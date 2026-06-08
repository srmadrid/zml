const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Cosh(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Cosh: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Cosh: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Cosh: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Cosh", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Cosh: " ++ @typeName(X) ++ " must implement `fn Cosh(type) type`");

            return X.Cosh(X);
        },
    }
}

/// Returns the hyperbolic cosine `cosh(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.cosh(x: X) numeric.Cosh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic cosine of.
///
/// ## Returns
/// `numeric.Cosh(@TypeOf(x))`: The hyperbolic cosine of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Cosh` method. The expected signature and
/// behavior of `Cosh` are as follows:
/// * `fn Cosh(type) type`: Returns the type of the hyperbolic cosine of `x`.
///
/// `numeric.Cosh(X)` or `X` must implement the required `cosh` method. The
/// expected signature and behavior of `cosh` are as follows:
/// * `fn cosh(X) numeric.Cosh(X)`: Returns the hyperbolic cosine of `x`.
pub fn cosh(x: anytype) numeric.Cosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Cosh(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.cosh(x),
        .dyadic => return dyadic.cosh(x),
        .complex => return complex.cosh(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "cosh",
                fn (X) numeric.Cosh(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.cosh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn cosh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.cosh(x);
        },
    }
}

/// Performs computation of the hyperbolic cosine `cosh(x)` of a numeric `x`
/// into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.coshInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the hyperbolic cosine of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `coshInto` method. The expected
/// signature and behavior of `coshInto` are as follows:
/// * `fn coshInto(*O, X) void`: Computes the hyperbolic cosine of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `coshInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.cosh`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn coshInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.coshInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "coshInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.coshInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "coshInto", fn (*O, X) void, &.{ *O, X }))
                return O.coshInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "coshInto", fn (*O, X) void, &.{ *O, X }))
            return X.coshInto(o, x);
    }

    return numeric.set(o, numeric.cosh(x));
}
