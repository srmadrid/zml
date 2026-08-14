const linalg = @import("../../linalg.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn dotUnchecked(x: anytype, y: anytype) linalg.Dot(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = linalg.Dot(X, Y);

    var sum = numeric.cast(meta.Accumulator(R), 0);

    inline for (0..X.len) |i| {
        // sum += conj(x[i]) * y[i]
        numeric.fmaInto(
            &sum,
            numeric.conj(x.data[i]),
            y.data[i],
            sum,
        );
    }

    return numeric.cast(R, sum);
}
