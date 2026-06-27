const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => switch (x.flags.noconj) {
            true => k_apply2IntoUnchecked(o, x, y, opInto, true, true),
            false => k_apply2IntoUnchecked(o, x, y, opInto, true, false),
        },
        false => switch (x.flags.noconj) {
            true => k_apply2IntoUnchecked(o, x, y, opInto, false, true),
            false => k_apply2IntoUnchecked(o, x, y, opInto, false, false),
        },
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_o: bool, comptime noconj_x: bool) void {
    const y_eff = if (comptime noconj_o) y else numeric.conj(y);

    const X: type = @TypeOf(x);

    if (comptime meta.layoutOf(meta.Child(@TypeOf(o))) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < j) : (i += 1) {
                const tx = if (comptime meta.uploOf(X) == .upper) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    y_eff,
                );
            }

            opInto(
                &o.data[o._index(j, j)],
                if (comptime noconj_o == noconj_x)
                    x.data[x._index(j, j)]
                else
                    numeric.conj(x.data[x._index(j, j)]),
                y_eff,
            );

            i = j + 1;
            while (i < o.rows) : (i += 1) {
                const tx = if (comptime meta.uploOf(X) == .lower) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    y_eff,
                );
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < i) : (j += 1) {
                const tx = if (comptime meta.uploOf(X) == .lower) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    y_eff,
                );
            }

            opInto(
                &o.data[o._index(i, i)],
                if (comptime noconj_o == noconj_x)
                    x.data[x._index(i, i)]
                else
                    numeric.conj(x.data[x._index(i, i)]),
                y_eff,
            );

            j = i + 1;
            while (j < o.cols) : (j += 1) {
                const tx = if (comptime meta.uploOf(X) == .upper) x.data[x._index(i, j)] else x.data[x._index(j, i)];
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    y_eff,
                );
            }
        }
    }
}
