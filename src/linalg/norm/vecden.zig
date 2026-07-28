const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

const linalg = @import("../../linalg.zig");

pub fn norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Norm(@TypeOf(x), order) {
    return switch (x.flags.noconj) {
        true => k_norm(x, order, true),
        false => k_norm(x, order, false),
    };
}

fn k_norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x))), comptime noconj_x: bool) linalg.Norm(@TypeOf(x), order) {
    const X = @TypeOf(x);
    const R = linalg.Norm(X, order);

    switch (comptime order) {
        .l1 => {
            var sum = numeric.zero(meta.Accumulator(R));

            if (x.inc == 1) {
                for (0..x.len) |i| {
                    // sum += abs(x[i])
                    numeric.addInto(
                        &sum,
                        sum,
                        numeric.abs(
                            if (comptime noconj_x)
                                x.data[i]
                            else
                                numeric.conj(x.data[i]),
                        ),
                    );
                }
            } else {
                var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
                for (0..x.len) |_| {
                    // sum += abs(x[ix])
                    numeric.addInto(
                        &sum,
                        sum,
                        numeric.abs(
                            if (comptime noconj_x)
                                x.data[numeric.cast(usize, ix)]
                            else
                                numeric.conj(x.data[numeric.cast(usize, ix)]),
                        ),
                    );

                    ix += x.inc;
                }
            }

            return numeric.cast(R, sum);
        },
        .l2, .frobenius => {
            var sum = numeric.zero(meta.Accumulator(R));

            if (x.inc == 1) {
                for (0..x.len) |i| {
                    // sum += abs(x[i])²
                    numeric.addInto(
                        &sum,
                        sum,
                        numeric.abs2(
                            if (comptime noconj_x)
                                x.data[i]
                            else
                                numeric.conj(x.data[i]),
                        ),
                    );
                }
            } else {
                var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
                for (0..x.len) |_| {
                    // sum += abs(x[ix])²
                    numeric.addInto(
                        &sum,
                        sum,
                        numeric.abs2(
                            if (comptime noconj_x)
                                x.data[numeric.cast(usize, ix)]
                            else
                                numeric.conj(x.data[numeric.cast(usize, ix)]),
                        ),
                    );

                    ix += x.inc;
                }
            }

            return numeric.cast(R, numeric.sqrt(sum));
        },
        .inf => {
            var max_val = numeric.zero(R);

            if (x.inc == 1) {
                for (0..x.len) |i| {
                    numeric.maxInto(
                        &max_val,
                        max_val,
                        numeric.abs(
                            if (comptime noconj_x)
                                x.data[i]
                            else
                                numeric.conj(x.data[i]),
                        ),
                    );
                }
            } else {
                var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
                for (0..x.len) |_| {
                    numeric.maxInto(
                        &max_val,
                        max_val,
                        numeric.abs(
                            if (comptime noconj_x)
                                x.data[numeric.cast(usize, ix)]
                            else
                                numeric.conj(x.data[numeric.cast(usize, ix)]),
                        ),
                    );

                    ix += x.inc;
                }
            }

            return numeric.cast(R, max_val);
        },
        .p => |p_| {
            var sum = numeric.zero(meta.Accumulator(R));

            if (x.inc == 1) {
                for (0..x.len) |i| {
                    // sum += abs(x[i])ᵖ
                    numeric.addInto(
                        &sum,
                        sum,
                        numeric.pow(
                            numeric.abs(
                                if (comptime noconj_x)
                                    x.data[i]
                                else
                                    numeric.conj(x.data[i]),
                            ),
                            p_,
                        ),
                    );
                }
            } else {
                var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
                for (0..x.len) |_| {
                    // sum += abs(x[i])ᵖ
                    numeric.addInto(
                        &sum,
                        sum,
                        numeric.pow(
                            numeric.abs(
                                if (comptime noconj_x)
                                    x.data[numeric.cast(usize, ix)]
                                else
                                    numeric.conj(x.data[numeric.cast(usize, ix)]),
                            ),
                            p_,
                        ),
                    );

                    ix += x.inc;
                }
            }

            return numeric.cast(R, numeric.pow(sum, numeric.div(1, p_)));
        },
    }
}
