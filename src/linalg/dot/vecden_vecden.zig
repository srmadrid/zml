const linalg = @import("../../linalg.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn dotUnchecked(x: anytype, y: anytype) linalg.Dot(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = linalg.Dot(X, Y);

    return numeric.cast(
        R,
        linalg.blas.dotc(
            x.len,
            x.data,
            x.inc,
            y.data,
            y.inc,
        ) catch unreachable,
    );
}
