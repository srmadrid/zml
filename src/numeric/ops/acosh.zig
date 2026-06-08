const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Acosh(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Acosh: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Acosh: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Acosh: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Acosh", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Acosh: " ++ @typeName(X) ++ " must implement `fn Acosh(type) type`");

            return X.Acosh(X);
        },
    }
}

/// Returns the hyperbolic arccosine `cos⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.acosh(x: X) numeric.Acosh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic arccosine of.
///
/// ## Returns
/// `numeric.Acosh(@TypeOf(x))`: The hyperbolic arccosine of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Acosh` method. The expected signature and
/// behavior of `Acosh` are as follows:
/// * `fn Acosh(type) type`: Returns the type of the hyperbolic arccosine of
///   `x`.
///
/// `numeric.Acosh(X)` or `X` must implement the required `acosh` method. The
/// expected signature and behavior of `acosh` are as follows:
/// * `fn acosh(X) numeric.Acosh(X)`: Returns the hyperbolic arccosine of `x`.
pub fn acosh(x: anytype) numeric.Acosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Acosh(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.acosh(x),
        .dyadic => return dyadic.acosh(x),
        .complex => return complex.acosh(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "acosh",
                fn (X) numeric.Acosh(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.acosh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn acosh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.acosh(x);
        },
    }
}

/// Performs computation of the hyperbolic arccosine `cosh⁻¹(x)` of a numeric
/// `x` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.acoshInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the hyperbolic arccosine of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `acoshInto` method. The expected
/// signature and behavior of `acoshInto` are as follows:
/// * `fn acoshInto(*O, X) void`: Computes the hyperbolic arccosine of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `acoshInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.acosh`, potentially resulting in a less efficient implementation.
/// In this case, `O` and `X` must adhere to the requirements of these functions.
pub fn acoshInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.acoshInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "acoshInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.acoshInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "acoshInto", fn (*O, X) void, &.{ *O, X }))
                return O.acoshInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "acoshInto", fn (*O, X) void, &.{ *O, X }))
            return X.acoshInto(o, x);
    }

    return numeric.set(o, numeric.acosh(x));
}
