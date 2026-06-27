const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => switch (y.flags.noconj) {
            true => k_apply2IntoUnchecked(o, x, y, opInto, true, true),
            false => k_apply2IntoUnchecked(o, x, y, opInto, true, false),
        },
        false => switch (y.flags.noconj) {
            true => k_apply2IntoUnchecked(o, x, y, opInto, false, true),
            false => k_apply2IntoUnchecked(o, x, y, opInto, false, false),
        },
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_o: bool, comptime noconj_y: bool) void {
    const x_eff = if (comptime noconj_o) x else numeric.conj(x);

    o.nnz = y.nnz;

    var i: usize = 0;
    while (i < o.nnz) : (i += 1) {
        opInto(
            &o.data[i],
            x_eff,
            if (comptime noconj_o == noconj_y)
                y.data[i]
            else
                numeric.conj(y.data[i]),
        );

        o.idx[i] = y.idx[i];
    }
}
