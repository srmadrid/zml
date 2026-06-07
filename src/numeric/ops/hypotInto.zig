const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs computation of the hypotenuse `√(x² + y²)` of two numerics `x` and
/// `y` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.hypotInto(o: *O, x: X, y: Y) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O`, `X` or `Y` should implement the required `hypotInto` method. The
/// expected signature and behavior of `hypotInto` are as follows:
/// * `fn hypotInto(*O, X, Y) void`: Computes the hypotenuse of `x` and `y` and
///   stores it in `o`.
///
/// If none of `O`, `X` and `Y` implement the required `hypotInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.hypot`, potentially resulting in a less efficient implementation. In
/// this case, `O`, `X` and `Y` must adhere to the requirements of these
/// functions.
pub fn hypotInto(o: anytype, x: anytype, y: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X) or
        !meta.isNumeric(Y))
        @compileError("zsl.numeric.hypotInto: o must be a mutable one-item pointer to a numeric, and x and y must be numerics, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // O, X and Y all custom
                if (comptime meta.anyHasMethod(&.{ O, X, Y }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            } else { // only O and X custom
                if (comptime meta.anyHasMethod(&.{ O, X }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            }
        } else {
            if (comptime meta.isCustomNumeric(Y)) { // only O and Y custom
                if (comptime meta.anyHasMethod(&.{ O, Y }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            } else { // only O custom
                if (comptime meta.hasMethod(O, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return O.hypotInto(o, x, y);
            }
        }
    } else {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) { // only X and Y custom
                if (comptime meta.anyHasMethod(&.{ X, Y }, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y })) |Impl|
                    return Impl.hypotInto(o, x, y);
            } else { // only X custom
                if (comptime meta.hasMethod(X, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                    return X.hypotInto(o, x, y);
            }
        } else if (comptime meta.isCustomNumeric(Y)) { // only Y custom
            if (comptime meta.hasMethod(Y, "hypotInto", fn (*O, X, Y) void, &.{ *O, X, Y }))
                return Y.hypotInto(o, x, y);
        }
    }

    return numeric.set(o, numeric.hypot(x, y));
}
