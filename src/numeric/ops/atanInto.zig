const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs computation of the arctangent `tan⁻¹(x)` of a numeric `x` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.atanInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the arctangent of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `atanInto` method. The expected
/// signature and behavior of `atanInto` are as follows:
/// * `fn atanInto(*O, X) void`: Computes the arctangent of `x` and stores it
///   in `o`.
///
/// If neither `O` nor `X` implement the required `atanInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.atan`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn atanInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.atanInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "atanInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.atanInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "atanInto", fn (*O, X) void, &.{ *O, X }))
                return O.atanInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "atanInto", fn (*O, X) void, &.{ *O, X }))
            return X.atanInto(o, x);
    }

    return numeric.set(o, numeric.atan(x));
}
