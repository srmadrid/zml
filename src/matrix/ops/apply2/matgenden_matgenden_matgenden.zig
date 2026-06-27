const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => switch (x.flags.noconj) {
            true => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, true, true, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, true, true, false),
            },
            false => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, true, false, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, true, false, false),
            },
        },
        false => switch (x.flags.noconj) {
            true => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, false, true, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, false, true, false),
            },
            false => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, false, false, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, false, false, false),
            },
        },
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_o: bool, comptime noconj_x: bool, comptime noconj_y: bool) void {
    if (comptime meta.layoutOf(meta.Child(@TypeOf(o))) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        x.data[x._index(i, j)]
                    else
                        numeric.conj(x.data[x._index(i, j)]),
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
                    if (comptime noconj_o == noconj_x)
                        x.data[x._index(i, j)]
                    else
                        numeric.conj(x.data[x._index(i, j)]),
                    if (comptime noconj_o == noconj_y)
                        y.data[y._index(i, j)]
                    else
                        numeric.conj(y.data[y._index(i, j)]),
                );
            }
        }
    }
}
