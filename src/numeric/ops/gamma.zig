const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Gamma(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Gamma: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Gamma: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Gamma: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Gamma", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Gamma: " ++ @typeName(X) ++ " must implement `fn Gamma(type) type`");

            return X.Gamma(X);
        },
    }
}

/// Returns the gamma function of a numeric `x`.
///
/// The gamma function is defined as:
/// $$
/// \Gamma(x) = \int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.gamma(x: X) numeric.Gamma(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the gamma function of.
///
/// ## Returns
/// `numeric.Gamma(@TypeOf(x))`: The gamma function  of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Gamma` method. The expected signature
/// and behavior of `Gamma` are as follows:
/// * `fn Gamma(type) type`: Returns the type of the gamma function of `x`.
///
/// `numeric.Gamma(X)` or `X` must implement the required `gamma` method. The
/// expected signature and behavior of `gamma` are as follows:
/// * `fn gamma(X) numeric.Gamma(X)`: Returns the gamma function of `x`.
pub fn gamma(x: anytype) numeric.Gamma(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Gamma(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.gamma(x),
        .dyadic => return dyadic.gamma(x),
        .complex => return complex.gamma(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "gamma",
                fn (X) numeric.Gamma(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.gamma: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn gamma(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.gamma(x);
        },
    }
}

/// Performs computation of the gamma function of a numeric `x` into a numeric
/// `o`.
///
/// The gamma function is defined as:
/// $$
/// \Gamma(x) = \int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.gammaInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the gamma function of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `gammaInto` method. The expected
/// signature and behavior of `gammaInto` are as follows:
/// * `fn gammaInto(*O, X) void`: Computes the gamma function of `x` and stores
///   it in `o`.
///
/// If neither `O` nor `X` implement the required `gammaInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.gamma`, potentially resulting in a less efficient implementation.
/// In this case, `O` and `X` must adhere to the requirements of these functions.
pub fn gammaInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.gammaInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "gammaInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.gammaInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "gammaInto", fn (*O, X) void, &.{ *O, X }))
                return O.gammaInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "gammaInto", fn (*O, X) void, &.{ *O, X }))
            return X.gammaInto(o, x);
    }

    return numeric.set(o, numeric.gamma(x));
}
