const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or !meta.isNumeric(Z) or
        !meta.numericType(X).le(.complex) or !meta.numericType(Y).le(.complex) or !meta.numericType(Z).le(.complex) or
        (meta.numericType(X) != .complex and meta.numericType(Y) != .complex and meta.numericType(Z) != .complex))
        @compileError("zsl.complex.Fma: at least one of X, Y or Z must be a complex type, the others must be bool, int, float, dyadic or complex types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tZ = " ++ @typeName(Z) ++ "\n");

    return complex.Coerce(X, numeric.Coerce(Y, Z));
}

/// Performs fused multiplication and addition (x * y + z) between three
/// operands of complex, dyadic, float, int or bool types, where at least one
/// operand must be of complex type. The result type is determined by coercing
/// the operand types, and the operation is performed by casting all three
/// operands to the result type, then performing the fused operation.
///
/// ## Signature
/// ```zig
/// complex.fma(x: X, y: Y, z: Z) complex.Fma(X, Y, Z)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left multiplication operand.
/// * `y` (`anytype`): The right multiplication operand.
/// * `z` (`anytype`): The addition operand.
///
/// ## Returns
/// `complex.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z))`: The result of the fused
/// multiplication and addition.
pub fn fma(x: anytype, y: anytype, z: anytype) complex.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Z: type = @TypeOf(z);
    const R: type = complex.Fma(X, Y, Z);

    if (comptime !meta.isComplex(X)) {
        if (comptime !meta.isComplex(Y)) {
            // x, y are real, z is complex
            return .{
                .re = numeric.fma(x, y, z.re),
                .im = numeric.cast(meta.Scalar(R), z.im),
            };
        } else {
            if (comptime !meta.isComplex(Z)) {
                // x, z are real, y is complex
                return .{
                    .re = numeric.fma(x, y.re, z),
                    .im = numeric.mul(x, y.im),
                };
            } else {
                // x is real, y, z are complex
                return .{
                    .re = numeric.fma(x, y.re, z.re),
                    .im = numeric.fma(x, y.im, z.im),
                };
            }
        }
    } else {
        if (comptime !meta.isComplex(Y)) {
            if (comptime !meta.isComplex(Z)) {
                // y, z are real, x is complex
                return .{
                    .re = numeric.fma(x.re, y, z),
                    .im = numeric.mul(x.im, y),
                };
            } else {
                // y is real, x, z are complex
                return .{
                    .re = numeric.fma(x.re, y, z.re),
                    .im = numeric.fma(x.im, y, z.im),
                };
            }
        } else {
            if (comptime !meta.isComplex(Z)) {
                // z is real, x, y are complex
                return .{
                    .re = numeric.fma(x.re, y.re, numeric.fma(numeric.neg(x.im), y.im, z)),
                    .im = numeric.fma(x.re, y.im, numeric.mul(x.im, y.re)),
                };
            } else {
                // x, y, z are complex
                return .{
                    .re = numeric.fma(x.re, y.re, numeric.fma(numeric.neg(x.im), y.im, z.re)),
                    .im = numeric.fma(x.re, y.im, numeric.fma(x.im, y.re, z.im)),
                };
            }
        }
    }
}
