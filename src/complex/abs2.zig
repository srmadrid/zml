const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

pub fn Abs2(comptime Z: type) type {
    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.Abs2: Z must be a complex type, got \n\tZ = " ++ @typeName(Z) ++ "\n");

    return meta.Scalar(Z);
}

/// Returns the squared absolute value of a complex `z`.
///
/// ## Signature
/// ```zig
/// complex.abs2(z: Z) complex.Abs2(Z)
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The complex value to get the squared absolute value of.
///
/// ## Returns
/// `complex.Abs2(@TypeOf(z))`: The squared absolute value of `z`.
pub fn abs2(z: anytype) complex.Abs2(@TypeOf(z)) {
    return numeric.fma(z.re, z.re, numeric.mul(z.im, z.im));
}
