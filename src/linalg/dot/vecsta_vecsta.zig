const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

const linalg = @import("../../linalg.zig");

pub fn dot(x: anytype, y: anytype) linalg.Dot(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = linalg.Dot(X, Y);

    var sum = numeric.zero(meta.Accumulator(R));

    inline for (0..X.len) |i| {
        // sum += x[i] * conj(y[i])
        numeric.fmaInto(
            &sum,
            x.data[i],
            numeric.conj(y.data[i]),
            sum,
        );
    }

    return numeric.cast(R, sum);
}
