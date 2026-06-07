const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs computation of the exponential `eˣ` of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.expInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the exponential of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `expInto` method. The expected
/// signature and behavior of `expInto` are as follows:
/// * `fn expInto(*O, X) void`: Computes the exponential of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `expInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.exp`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn expInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.expInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "expInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.expInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "expInto", fn (*O, X) void, &.{ *O, X }))
                return O.expInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "expInto", fn (*O, X) void, &.{ *O, X }))
            return X.expInto(o, x);
    }

    return numeric.set(o, numeric.exp(x));
}
