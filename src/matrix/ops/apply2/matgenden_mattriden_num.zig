const int = @import("../../../int.zig");
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

    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            if (comptime meta.uploOf(X) == .upper) {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        if (comptime noconj_o == noconj_x)
                            x.data[x._index(i, j)]
                        else
                            numeric.conj(x.data[x._index(i, j)]),
                        y_eff,
                    );
                }
            } else {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        numeric.cast(meta.Numeric(X), 0),
                        y_eff,
                    );
                }
            }

            if (j < o.rows) {
                if (comptime meta.diagOf(X) == .unit)
                    opInto(
                        &o.data[o._index(j, j)],
                        numeric.cast(meta.Numeric(X), 1),
                        y_eff,
                    )
                else
                    opInto(
                        &o.data[o._index(j, j)],
                        if (comptime noconj_o == noconj_x)
                            x.data[x._index(j, j)]
                        else
                            numeric.conj(x.data[x._index(j, j)]),
                        y_eff,
                    );
            }

            i = int.min(j + 1, o.rows);
            if (comptime meta.uploOf(X) == .lower) {
                while (i < o.rows) : (i += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        if (comptime noconj_o == noconj_x)
                            x.data[x._index(i, j)]
                        else
                            numeric.conj(x.data[x._index(i, j)]),
                        y_eff,
                    );
                }
            } else {
                while (i < o.rows) : (i += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        numeric.cast(meta.Numeric(X), 0),
                        y_eff,
                    );
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            if (comptime meta.uploOf(X) == .lower) {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        if (comptime noconj_o == noconj_x)
                            x.data[x._index(i, j)]
                        else
                            numeric.conj(x.data[x._index(i, j)]),
                        y_eff,
                    );
                }
            } else {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        numeric.cast(meta.Numeric(X), 0),
                        y_eff,
                    );
                }
            }

            if (i < o.cols) {
                if (comptime meta.diagOf(X) == .unit)
                    opInto(
                        &o.data[o._index(i, i)],
                        numeric.cast(meta.Numeric(X), 1),
                        y_eff,
                    )
                else
                    opInto(
                        &o.data[o._index(i, i)],
                        if (comptime noconj_o == noconj_x)
                            x.data[x._index(i, i)]
                        else
                            numeric.conj(x.data[x._index(i, i)]),
                        y_eff,
                    );
            }

            j = int.min(i + 1, o.cols);
            if (comptime meta.uploOf(X) == .upper) {
                while (j < o.cols) : (j += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        if (comptime noconj_o == noconj_x)
                            x.data[x._index(i, j)]
                        else
                            numeric.conj(x.data[x._index(i, j)]),
                        y_eff,
                    );
                }
            } else {
                while (j < o.cols) : (j += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        numeric.cast(meta.Numeric(X), 0),
                        y_eff,
                    );
                }
            }
        }
    }
}
