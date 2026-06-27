const int = @import("../../../int.zig");
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
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < int.min(j, o.rows)) : (i += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    x_eff,
                    numeric.zero(meta.Numeric(Y)),
                );
            }

            if (j < o.rows) {
                opInto(
                    &o.data[o._index(j, j)],
                    x_eff,
                    if (comptime noconj_o == noconj_y)
                        y.data[j]
                    else
                        numeric.conj(y.data[j]),
                );
            }

            i = int.min(j + 1, o.rows);
            while (i < o.rows) : (i += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    x_eff,
                    numeric.zero(meta.Numeric(Y)),
                );
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < int.min(i, o.cols)) : (j += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    x_eff,
                    numeric.zero(meta.Numeric(Y)),
                );
            }

            if (i < o.cols) {
                opInto(
                    &o.data[o._index(i, i)],
                    x_eff,
                    if (comptime noconj_o == noconj_y)
                        y.data[i]
                    else
                        numeric.conj(y.data[i]),
                );
            }

            j = int.min(i + 1, o.cols);
            while (j < o.cols) : (j += 1) {
                opInto(
                    &o.data[o._index(i, j)],
                    x_eff,
                    numeric.zero(meta.Numeric(Y)),
                );
            }
        }
    }
}
