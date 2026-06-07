const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs computation of the complex conjugate of a numeric `x` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.conjInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the complex conjugate of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `conjInto` method. The expected
/// signature and behavior of `conjInto` are as follows:
/// * `fn conjInto(*O, X) void`: Computes the complex conjugate of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `conjInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.conj`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn conjInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.conjInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "conjInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.conjInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "conjInto", fn (*O, X) void, &.{ *O, X }))
                return O.conjInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "conjInto", fn (*O, X) void, &.{ *O, X }))
            return X.conjInto(o, x);
    }

    return numeric.set(o, numeric.conj(x));
}
