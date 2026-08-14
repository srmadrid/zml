const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

/// Returns the hyperbolic arctangent `tanh⁻¹(z)` of a complex.
///
/// ## Signature
/// ```zig
/// complex.atanh(z: Z) Z
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The value to get the hyperbolic arctangent of.
///
/// ## Returns
/// `@TypeOf(z)`: The hyperbolic arctangnet of `z`.
pub fn atanh(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.atanh: z must be a complex, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return numeric.div(
        numeric.sub(
            numeric.ln(numeric.add(1, z)),
            numeric.ln(numeric.sub(1, z)),
        ),
        2,
    );
}
