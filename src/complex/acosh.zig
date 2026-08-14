const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

/// Returns the hyperbolic arccosine `cosh⁻¹(z)` of a complex.
///
/// ## Signature
/// ```zig
/// complex.acosh(z: Z) Z
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The value to get the hyperbolic arccosine of.
///
/// ## Returns
/// `@TypeOf(z)`: The hyperbolic arccosine of `z`.
pub fn acosh(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.acosh: z must be a complex, got \n\tz: " ++ @typeName(Z) ++ "\n");

    return numeric.mul(
        numeric.ln(
            numeric.add(
                numeric.sqrt(numeric.div(numeric.add(z, 1), 2)),
                numeric.sqrt(numeric.div(numeric.sub(z, 1), 2)),
            ),
        ),
        2,
    );
}
