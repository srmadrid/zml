const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Erf(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Erf: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Erf: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Erf: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Erf", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Erf: " ++ @typeName(X) ++ " must implement `fn Erf(type) type`");

            return X.Erf(X);
        },
    }
}

/// Returns the error function of a numeric `x`.
///
/// The error function is defined as:
/// $$
/// \mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.erf(x: X) numeric.Erf(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the error function of.
///
/// ## Returns
/// `numeric.Erf(@TypeOf(x))`: The error function of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Erf` method. The expected signature and
/// behavior of `Erf` are as follows:
/// * `fn Erf(type) type`: Returns the type of the error function of `x`.
///
/// `numeric.Erf(X)` or `X` must implement the required `erf` method. The
/// expected signature and behavior of `erf` are as follows:
/// * `fn erf(X) numeric.Erf(X)`: Returns the error function of `x`.
pub fn erf(x: anytype) numeric.Erf(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Erf(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.erf(x),
        .dyadic => return dyadic.erf(x),
        .complex => return complex.erf(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "erf",
                fn (X) numeric.Erf(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.erf: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn erf(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.erf(x);
        },
    }
}

/// Performs computation of the error function of a numeric `x` into a numeric
/// `o`.
///
/// The error function is defined as:
/// $$
/// \mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.erfInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the error function of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `erfInto` method. The expected
/// signature and behavior of `erfInto` are as follows:
/// * `fn erfInto(*O, X) void`: Computes the error function of `x` and stores it
///   in `o`.
///
/// If neither `O` nor `X` implement the required `erfInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.erf`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn erfInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.erfInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "erfInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.erfInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "erfInto", fn (*O, X) void, &.{ *O, X }))
                return O.erfInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "erfInto", fn (*O, X) void, &.{ *O, X }))
            return X.erfInto(o, x);
    }

    return numeric.set(o, numeric.erf(x));
}
