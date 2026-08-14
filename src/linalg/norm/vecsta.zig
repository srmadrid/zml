const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

const linalg = @import("../../linalg.zig");

pub fn norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Norm(@TypeOf(x), order) {
    const X = @TypeOf(x);
    const R = linalg.Norm(X, order);

    switch (comptime order) {
        .l1 => {
            var sum = numeric.cast(meta.Accumulator(R), 0);

            inline for (0..X.len) |i| {
                // sum += abs(x[i])
                numeric.addInto(
                    &sum,
                    sum,
                    numeric.abs(x.data[i]),
                );
            }

            return numeric.cast(R, sum);
        },
        .l2, .frobenius => {
            var sum = numeric.cast(meta.Accumulator(R), 0);

            inline for (0..X.len) |i| {
                // sum += abs(x[i])²
                numeric.addInto(
                    &sum,
                    sum,
                    numeric.abs2(x.data[i]),
                );
            }

            return numeric.cast(R, numeric.sqrt(sum));
        },
        .inf => {
            var max_val = numeric.cast(R, 0);

            inline for (0..X.len) |i| {
                numeric.maxInto(
                    &max_val,
                    max_val,
                    numeric.abs(x.data[i]),
                );
            }

            return numeric.cast(R, max_val);
        },
        .p => |p_| {
            var sum = numeric.cast(meta.Accumulator(R), 0);

            inline for (0..X.len) |i| {
                // sum += abs(x[i])ᵖ
                numeric.addInto(
                    &sum,
                    sum,
                    numeric.pow(numeric.abs(x.data[i]), p_),
                );
            }

            return numeric.cast(R, numeric.pow(sum, numeric.div(1, p_)));
        },
    }
}
