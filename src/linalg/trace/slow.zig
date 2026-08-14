const linalg = @import("../../linalg.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn traceUnchecked(x: anytype) linalg.Trace(@TypeOf(x)) {
    const X = @TypeOf(x);
    const R = linalg.Trace(X);

    var sum = numeric.cast(meta.Accumulator(R), 0);

    const len = if (comptime meta.isStaticMatrix(X)) X.rows else x.rows;

    for (0..len) |i| {
        // sum += x[i, i]
        numeric.addInto(
            &sum,
            x.get(i, i) catch unreachable,
            sum,
        );
    }

    return numeric.cast(R, sum);
}
