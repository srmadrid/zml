const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

pub fn Arg(comptime Z: type) type {
    comptime if (!meta.isNumeric(Z) or meta.numericType(Z) != .complex)
        @compileError("zsl.complex.Arg: z must be a complex type, got \n\tZ = " ++ @typeName(Z) ++ "\n");

    return meta.Scalar(Z);
}

/// Returns the argument of a complex `z`.
///
/// ## Signature
/// ```zig
/// complex.arg(z: Z) complex.Arg(Z)
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The complex value to get the argument of.
///
/// ## Returns
/// `complex.Arg(@TypeOf(z))`: The argument of `z`.
pub fn arg(z: anytype) complex.Arg(@TypeOf(z)) {
    return numeric.atan2(z.im, z.re);
}
