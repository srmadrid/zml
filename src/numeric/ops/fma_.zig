const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

/// Performs in-place computation of the fused multiplication and addition of
/// three numerics `x`, `y` and `z` into a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.fma_(o: *O, x: X, y: Y, z: Z) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left multiplication operand.
/// * `y` (`anytype`): The right multiplication operand.
/// * `z` (`anytype`): The addition operand.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O`, `X`, `Y` or `Z` should implement the required `fma_` method. The
/// expected signature and behavior of `fma_` are as follows:
/// * `fn fma_(*O, X, Y, Z) void`: Computes the fused multiplication and
///   addition of `x`, `y` and `z` and stores it in `o`.
///
/// If none of `O`, `X`, `Y` and `Z` implement the required `fma_` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.fma`, potentially resulting in a less efficient implementation. In
/// this case, `O`, `X`, `Y` and `Z` must adhere to the requirements of these
/// functions.
pub fn fma_(o: anytype, x: anytype, y: anytype, z: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Z: type = @TypeOf(z);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X) or
        !meta.isNumeric(Y) or
        !meta.isNumeric(Z))
        @compileError("zsl.numeric.fma_: o must be a mutable one-item pointer to a numeric, and x, y and < must be numerics, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tz: " ++ @typeName(Z) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) {
            if (comptime meta.isCustomNumeric(Y)) {
                if (comptime meta.isCustomNumeric(Z)) { // O, X, Y and Z all custom
                    if (comptime meta.anyHasMethod(&.{ O, X, Y, Z }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fma_(o, x, y, z);
                } else { // only O, X and Y custom
                    if (comptime meta.anyHasMethod(&.{ O, X, Y }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fma_(o, x, y, z);
                }
            } else {
                if (comptime meta.isCustomNumeric(Z)) { // only O, X and Z custom
                    if (comptime meta.anyHasMethod(&.{ O, X, Z }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fma_(o, x, y, z);
                } else { // only O and X custom
                    if (comptime meta.anyHasMethod(&.{ O, X }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fma_(o, x, y, z);
                }
            }
        } else {
            if (comptime meta.isCustomNumeric(Y)) {
                if (comptime meta.isCustomNumeric(Z)) { // only O, Y and Z custom
                    if (comptime meta.anyHasMethod(&.{ O, Y, Z }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fma_(o, x, y, z);
                } else { // only O and Y custom
                    if (comptime meta.anyHasMethod(&.{ O, Y }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fma_(o, x, y, z);
                }
            } else {
                if (comptime meta.isCustomNumeric(Z)) { // only O and Z custom
                    if (comptime meta.anyHasMethod(&.{ O, Z }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                        return Impl.fma_(o, x, y, z);
                } else { // only O custom
                    if (comptime meta.hasMethod(O, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
                        return O.fma_(o, x, y, z);
                }
            }
        }
    } else if (comptime meta.isCustomNumeric(X)) {
        if (comptime meta.isCustomNumeric(Y)) {
            if (comptime meta.isCustomNumeric(Z)) { // only X, Y and Z custom
                if (comptime meta.anyHasMethod(&.{ X, Y, Z }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                    return Impl.fma_(o, x, y, z);
            } else { // only X and Y custom
                if (comptime meta.anyHasMethod(&.{ X, Y }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                    return Impl.fma_(o, x, y, z);
            }
        } else {
            if (comptime meta.isCustomNumeric(Z)) { // only X and Z custom
                if (comptime meta.anyHasMethod(&.{ X, Z }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                    return Impl.fma_(o, x, y, z);
            } else { // only X custom
                if (comptime meta.hasMethod(X, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
                    return X.fma_(o, x, y, z);
            }
        }
    } else if (comptime meta.isCustomNumeric(Y)) {
        if (comptime meta.isCustomNumeric(Z)) { // only Y and Z custom
            if (comptime meta.anyHasMethod(&.{ Y, Z }, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z })) |Impl|
                return Impl.fma_(o, x, y, z);
        } else { // only Y custom
            if (comptime meta.hasMethod(Y, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
                return Y.fma_(o, x, y, z);
        }
    } else if (comptime meta.isCustomNumeric(Z)) { // only Z custom
        if (comptime meta.hasMethod(Z, "fma_", fn (*O, X, Y, Z) void, &.{ *O, X, Y, Z }))
            return Z.fma_(o, x, y, z);
    }

    return numeric.set(o, numeric.fma(x, y, z));
}
