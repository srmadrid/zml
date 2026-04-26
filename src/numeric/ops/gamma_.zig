const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs in-place computation of the gamma function of a numeric `x` into a
/// numeric `o`.
///
/// The gamma function is defined as:
/// $$
/// \Gamma(x) = \int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.gamma_(o: *O, x: X) void
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
/// `O` or `X` should implement the required `gamma_` method. The expected
/// signature and behavior of `gamma_` are as follows:
/// * `fn gamma_(*O, X) void`: Computes the gamma function of `x` and stores
///   it in `o`.
///
/// If neither `O` nor `X` implement the required `gamma_` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.gamma`, potentially resulting in a less efficient implementation.
/// In this case, `O` and `X` must adhere to the requirements of these functions.
pub fn gamma_(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.gamma_: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "gamma_", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.gamma_(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "gamma_", fn (*O, X) void, &.{ *O, X }))
                return O.gamma_(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "gamma_", fn (*O, X) void, &.{ *O, X }))
            return X.gamma_(o, x);
    }

    return numeric.set(o, numeric.gamma(x));
}
