const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs computation of the absolute value of a numeric `x` into a numeric
/// `o`.
///
/// ## Signature
/// ```zig
/// numeric.absInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the absolute value of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `absInto` method. The expected
/// signature and behavior of `absInto` are as follows:
/// * `fn absInto(*O, X) void`: Computes the absolute value of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `absInto` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.abs`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn absInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.absInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "absInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.absInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "absInto", fn (*O, X) void, &.{ *O, X }))
                return O.absInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "absInto", fn (*O, X) void, &.{ *O, X }))
            return X.absInto(o, x);
    }

    return numeric.set(o, numeric.abs(x));
}
