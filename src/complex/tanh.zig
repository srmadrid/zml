const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

/// Returns the hyperbolic tangent `tanh(z)` of a complex.
///
/// ## Signature
/// ```zig
/// complex.tanh(z: Z) Z
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The value to get the hyperbolic tangent of.
///
/// ## Returns
/// `@TypeOf(z)`: The hyperbolic tangent of `z`.
pub fn tanh(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.tanh: z must be a complex, got \n\tz: " ++ @typeName(Z) ++ "\n");

    if (numeric.ge(z.re, 0)) {
        const w = numeric.exp(numeric.mul(z, -2));

        return numeric.div(
            numeric.sub(1, w),
            numeric.add(1, w),
        );
    } else {
        const w = numeric.exp(numeric.mul(z, 2));

        return numeric.div(
            numeric.sub(w, 1),
            numeric.add(w, 1),
        );
    }
}
