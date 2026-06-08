const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Erfc(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Erfc: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Erfc: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Erfc: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Erfc", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Erfc: " ++ @typeName(X) ++ " must implement `fn Erfc(type) type`");

            return X.Erfc(X);
        },
    }
}

/// Returns the complementary error function of a numeric `x`.
///
/// The error function is defined as:
/// $$
/// \mathrm{erfc}(x) = 1 - \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.erfc(x: X) numeric.Erfc(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the complementary error function
///   of.
///
/// ## Returns
/// `numeric.Erfc(@TypeOf(x))`: The error function of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Erfc` method. The expected signature and
/// behavior of `Erfc` are as follows:
/// * `fn Erfc(type) type`: Returns the type of the error function of `x`.
///
/// `numeric.Erfc(X)` or `X` must implement the required `erfc` method. The
/// expected signature and behavior of `erfc` are as follows:
/// * `fn erfc(X) numeric.Erfc(X)`: Returns the error function of `x`.
pub fn erfc(x: anytype) numeric.Erfc(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Erfc(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.erfc(x),
        .dyadic => return dyadic.erfc(x),
        .complex => return complex.erfc(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "erfc",
                fn (X) numeric.Erfc(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.erfc: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn erfc(" ++ @typeName(X) ++ " ) " ++ @typeName(R) ++ "`");

            return Impl.erfc(x);
        },
    }
}

/// Performs computation of the complementary error function of a numeric `x`
/// into a numeric `o`.
///
/// The complementary error function is defined as:
/// $$
/// \mathrm{erfc}(x) = 1 - \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.erfcInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the complementary error function
///   of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `erfcInto` method. The expected
/// signature and behavior of `erfcInto` are as follows:
/// * `fn erfcInto(*O, X) void`: Computes the complementary error function of
///   `x` and stores it in `o`.
///
/// If neither `O` nor `X` implement the required `erfcInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.erfc`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn erfcInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.erfcInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "erfcInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.erfcInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "erfcInto", fn (*O, X) void, &.{ *O, X }))
                return O.erfcInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "erfcInto", fn (*O, X) void, &.{ *O, X }))
            return X.erfcInto(o, x);
    }

    return numeric.set(o, numeric.erfc(x));
}
