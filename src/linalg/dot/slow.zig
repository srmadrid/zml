const linalg = @import("../../linalg.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn dotUnchecked(x: anytype, y: anytype) linalg.Dot(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = linalg.Dot(X, Y);

    var sum = numeric.cast(meta.Accumulator(R), 0);

    const len = if (comptime meta.isStaticVector(X)) X.len else x.len;

    for (0..len) |i| {
        // sum += x[i] * conj(y[i])
        numeric.fmaInto(
            &sum,
            x.getAssumeInBounds(i),
            numeric.conj(y.getAssumeInBounds(i)),
            sum,
        );
    }

    return numeric.cast(R, sum);
}
