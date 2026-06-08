const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

pub fn Abs1(comptime Z: type) type {
    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.Abs1: Z must be a complex type, got \n\tZ = " ++ @typeName(Z) ++ "\n");

    return meta.Scalar(Z);
}

/// Returns the 1-norm of a complex `z`.
///
/// ## Signature
/// ```zig
/// complex.abs1(z: Z) complex.Abs1(Z)
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The complex value to get the 1-norm of.
///
/// ## Returns
/// `complex.Abs1(@TypeOf(z))`: The 1-norm of `z`.
pub fn abs1(z: anytype) complex.Abs1(@TypeOf(z)) {
    return numeric.add(numeric.abs(z.re), numeric.abs(z.im));
}
