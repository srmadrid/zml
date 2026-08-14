const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const complex = @import("../complex.zig");

/// Returns the square root `√z` of a complex.
///
/// ## Signature
/// ```zig
/// complex.sqrt(z: Z) Z
/// ```
///
/// ## Arguments
/// * `z` (`anytype`): The value to get the square root of.
///
/// ## Returns
/// `@TypeOf(z)`: The square root of `z`.
pub fn sqrt(z: anytype) @TypeOf(z) {
    if (numeric.eq(z.re, 0) and
        numeric.eq(z.im, 0))
        return z;

    const a = numeric.abs(z.re);
    const b = numeric.abs(z.im);
    const r = numeric.abs(z);

    const w = numeric.sqrt(numeric.div(numeric.add(r, a), 2));

    if (numeric.ge(z.re, 0)) {
        return .{
            .re = w,
            .im = numeric.div(z.im, numeric.mul(2, w)),
        };
    } else {
        return .{
            .re = numeric.div(b, numeric.mul(2, w)),
            .im = if (numeric.ge(z.im, 0)) w else numeric.neg(w),
        };
    }
}
