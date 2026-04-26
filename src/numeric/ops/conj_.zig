const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs in-place computation of the complex conjugate of a numeric `x` into
/// a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.conj_(o: *O, x: X) void
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
/// `O` or `X` should implement the required `conj_` method. The expected
/// signature and behavior of `conj_` are as follows:
/// * `fn conj_(*O, X) void`: Computes the complex conjugate of `x` and stores
///   it in `o`.
///
/// If neither `O` nor `X` implement the required `conj_` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.conj`,
/// potentially resulting in a less efficient implementation. In this case, `O`
/// and `X` must adhere to the requirements of these functions.
pub fn conj_(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.conj_: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "conj_", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.conj_(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "conj_", fn (*O, X) void, &.{ *O, X }))
                return O.conj_(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "conj_", fn (*O, X) void, &.{ *O, X }))
            return X.conj_(o, x);
    }

    return numeric.set(o, numeric.conj(x));
}
