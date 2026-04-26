const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs in-place computation of the hyperbolic tangent `tanh(x)` of a
/// numeric `x` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.tanh_(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the hyperbolic tangent of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `tanh_` method. The expected
/// signature and behavior of `tanh_` are as follows:
/// * `fn tanh_(*O, X) void`: Computes the hyperbolic tangent of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `tanh_` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.tanh`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn tanh_(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.tanh_: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "tanh_", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.tanh_(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "tanh_", fn (*O, X) void, &.{ *O, X }))
                return O.tanh_(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "tanh_", fn (*O, X) void, &.{ *O, X }))
            return X.tanh_(o, x);
    }

    return numeric.set(o, numeric.tanh(x));
}
