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
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(meta.Child(@TypeOf(o))) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            while (i < j) : (i += 1) {
                const tx = if (comptime meta.uploOf(X) == .upper) x.data[x._index(i, j)] else numeric.conj(x.data[x._index(j, i)]);
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    numeric.cast(meta.Numeric(Y), 0),
                );
            }

            opInto(
                &o.data[o._index(j, j)],
                if (comptime noconj_o == noconj_x)
                    x.data[x._index(j, j)]
                else
                    numeric.conj(x.data[x._index(j, j)]),
                if (comptime noconj_o == noconj_y)
                    y.data[j]
                else
                    numeric.conj(y.data[j]),
            );

            i = j + 1;
            while (i < o.rows) : (i += 1) {
                const tx = if (comptime meta.uploOf(X) == .lower) x.data[x._index(i, j)] else numeric.conj(x.data[x._index(j, i)]);
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    numeric.cast(meta.Numeric(Y), 0),
                );
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            while (j < i) : (j += 1) {
                const tx = if (comptime meta.uploOf(X) == .lower) x.data[x._index(i, j)] else numeric.conj(x.data[x._index(j, i)]);
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    numeric.cast(meta.Numeric(Y), 0),
                );
            }

            opInto(
                &o.data[o._index(i, i)],
                if (comptime noconj_o == noconj_x)
                    x.data[x._index(i, i)]
                else
                    numeric.conj(x.data[x._index(i, i)]),
                if (comptime noconj_o == noconj_y)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );

            j = i + 1;
            while (j < o.cols) : (j += 1) {
                const tx = if (comptime meta.uploOf(X) == .upper) x.data[x._index(i, j)] else numeric.conj(x.data[x._index(j, i)]);
                opInto(
                    &o.data[o._index(i, j)],
                    if (comptime noconj_o == noconj_x)
                        tx
                    else
                        numeric.conj(tx),
                    numeric.cast(meta.Numeric(Y), 0),
                );
            }
        }
    }
}
