const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

/// Returns the tangent `tan(z)` of a complex.
///
/// ## Signature
/// ```zig
/// complex.tan(z: Z) Z
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The value to get the tangent of.
///
/// ## Returns
/// `@TypeOf(z)`: The tangent of `z`.
pub fn tan(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.tan: z must be a complex, got \n\tz: " ++ @typeName(Z) ++ "\n");

    const two_re = numeric.add(z.re, z.re);
    const two_im = numeric.add(z.im, z.im);

    const E = numeric.exp(numeric.neg(numeric.abs(two_im)));
    const E_sq = numeric.mul(E, E);
    const two_E = numeric.add(E, E);

    const sin_two_re = numeric.sin(two_re);
    const cos_two_re = numeric.cos(two_re);

    const den = numeric.add(1, numeric.fma(two_E, cos_two_re, E_sq));
    const num_re = numeric.mul(two_E, sin_two_re);
    const num_im = numeric.mul(numeric.sign(z.im), numeric.sub(1, E_sq));

    return .{
        .re = numeric.div(num_re, den),
        .im = numeric.div(num_im, den),
    };
}
