const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Asinh(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Asinh: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Asinh: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Asinh: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Asinh", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Asinh: " ++ @typeName(X) ++ " must implement `fn Asinh(type) type`");

            return X.Asinh(X);
        },
    }
}

/// Returns the hyperbolic arcsine `sinh⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.asinh(x: X) numeric.Asinh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic arcsine of.
///
/// ## Returns
/// `numeric.Asinh(@TypeOf(x))`: The hyperbolic arcsine of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Asinh` method. The expected signature and
/// behavior of `Asinh` are as follows:
/// * `fn Asinh(type) type`: Returns the type of the hyperbolic arcsine of
///   `x`.
///
/// `numeric.Asinh(X)` or `X` must implement the required `asinh` method. The
/// expected signature and behavior of `asinh` are as follows:
/// * `fn asinh(X) numeric.Asinh(X)`: Returns the hyperbolic arcsine of `x`.
pub fn asinh(x: anytype) numeric.Asinh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Asinh(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.asinh(x),
        .dyadic => return dyadic.asinh(x),
        .complex => return complex.asinh(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "asinh",
                fn (X) numeric.Asinh(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.asinh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn asinh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.asinh(x);
        },
    }
}

/// Performs computation of the hyperbolic arcsine `sinh⁻¹(x)` of a numeric `x`
/// into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.asinhInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the hyperbolic arcsine of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `asinhInto` method. The expected
/// signature and behavior of `asinhInto` are as follows:
/// * `fn asinhInto(*O, X) void`: Computes the hyperbolic arcsine of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `asinhInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.asinh`, potentially resulting in a less efficient implementation.
/// In this case, `O` and `X` must adhere to the requirements of these functions.
pub fn asinhInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.asinhInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "asinhInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.asinhInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "asinhInto", fn (*O, X) void, &.{ *O, X }))
                return O.asinhInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "asinhInto", fn (*O, X) void, &.{ *O, X }))
            return X.asinhInto(o, x);
    }

    return numeric.set(o, numeric.asinh(x));
}
