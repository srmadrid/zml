const meta = @import("../../../meta.zig");
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

    const O: type = meta.Child(@TypeOf(o));

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    x_eff,
                    if (comptime noconj_o == noconj_y)
                        y.data[y._index(i, j)]
                    else
                        numeric.conj(y.data[y._index(i, j)]),
                );
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    x_eff,
                    if (comptime noconj_o == noconj_y)
                        y.data[y._index(i, j)]
                    else
                        numeric.conj(y.data[y._index(i, j)]),
                );
            }
        }
    }
}
