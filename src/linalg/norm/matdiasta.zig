const int = @import("../../int.zig");
const linalg = @import("../../linalg.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Norm(@TypeOf(x), order) {
    const X = @TypeOf(x);
    const R = linalg.Norm(X, order);

    switch (comptime order) {
        .l1, .l2, .inf, .p => {
            var max_val = numeric.zero(R);

            inline for (0..int.min(X.rows, X.cols)) |i| {
                numeric.maxInto(
                    &max_val,
                    max_val,
                    numeric.abs(x.data[i]),
                );
            }

            return numeric.cast(R, max_val);
        },
        .frobenius => {
            var sum = numeric.zero(meta.Accumulator(R));

            inline for (0..int.min(X.rows, X.cols)) |i| {
                // sum += abs(x[i])²
                numeric.addInto(
                    &sum,
                    sum,
                    numeric.abs2(x.data[i]),
                );
            }

            return numeric.cast(R, numeric.sqrt(sum));
        },
    }
}
