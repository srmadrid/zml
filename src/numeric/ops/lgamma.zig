const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Lgamma(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Lgamma: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Lgamma: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Lgamma: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Lgamma", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Lgamma: " ++ @typeName(X) ++ " must implement `fn Lgamma(type) type`");

            return X.Lgamma(X);
        },
    }
}

/// Returns the log-gamma function of a numeric `x`.
///
/// The log-gamma function is defined as:
/// $$
/// \log(\Gamma(x)) = \left(\int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t\right).
/// $$
///
/// ## Signature
/// ```zig
/// numeric.lgamma(x: X) numeric.Lgamma(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the log-gamma function of.
///
/// ## Returns
/// `numeric.Lgamma(@TypeOf(x))`: The log-gamma function  of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Lgamma` method. The expected signature
/// and behavior of `Lgamma` are as follows:
/// * `fn Lgamma(type) type`: Returns the type of the log-gamma function of
///   `x`.
///
/// `numeric.Lgamma(X)` or `X` must implement the required `lgamma` method.
/// The expected signature and behavior of `lgamma` are as follows:
/// * `fn lgamma(X) numeric.Lgamma(X)`: Returns the log-gamma function of `x`.
pub fn lgamma(x: anytype) numeric.Lgamma(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Lgamma(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.lgamma(x),
        .dyadic => return dyadic.lgamma(x),
        .complex => return complex.lgamma(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "lgamma",
                fn (X) numeric.Lgamma(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.lgamma: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn lgamma(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.lgamma(x);
        },
    }
}

/// Performs in-place computation of the log-gamma function of a numeric `x`
/// into a numeric `o`.
///
/// The log-gamma function is defined as:
/// $$
/// \log(\Gamma(x)) = \left(\int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t\right).
/// $$
///
/// ## Signature
/// ```zig
/// numeric.lgammaInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the log-gamma function of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `lgammaInto` method. The expected
/// signature and behavior of `lgammaInto` are as follows:
/// * `fn lgammaInto(*O, X) void`: Computes the log-gamma function of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `lgammaInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.lgamma`, potentially resulting in a less efficient implementation.
/// In this case, `O` and `X` must adhere to the requirements of these functions.
pub fn lgammaInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.lgammaInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "lgammaInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.lgammaInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "lgammaInto", fn (*O, X) void, &.{ *O, X }))
                return O.lgammaInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "lgammaInto", fn (*O, X) void, &.{ *O, X }))
            return X.lgammaInto(o, x);
    }

    return numeric.set(o, numeric.lgamma(x));
}
