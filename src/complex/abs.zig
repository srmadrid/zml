const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

pub fn Abs(comptime Z: type) type {
    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.Abs: Z must be a complex type, got \n\tZ = " ++ @typeName(Z) ++ "\n");

    return meta.Scalar(Z);
}

/// Returns the absolute value of a complex `z`.
///
/// ## Signature
/// ```zig
/// complex.abs(z: Z) complex.Abs(Z)
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The complex value to get the absolute value of.
///
/// ## Returns
/// `complex.Abs(@TypeOf(z))`: The absolute value of `z`.
pub fn abs(z: anytype) complex.Abs(@TypeOf(z)) {
    return numeric.hypot(z.re, z.im);
}
