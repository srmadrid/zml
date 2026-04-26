const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs in-place computation of the 1-norm of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.abs1_(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the 1-norm of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `abs1_` method. The expected
/// signature and behavior of `abs1_` are as follows:
/// * `fn abs1_(*O, X) void`: Computes the 1-norm of `x` and stores it in `o`.
///
/// If neither `O` nor `X` implement the required `abs1_` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.abs1`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn abs1_(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.abs1_: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "abs1_", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.abs1_(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "abs1_", fn (*O, X) void, &.{ *O, X }))
                return O.abs1_(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "abs1_", fn (*O, X) void, &.{ *O, X }))
            return X.abs1_(o, x);
    }

    return numeric.set(o, numeric.abs1(x));
}
