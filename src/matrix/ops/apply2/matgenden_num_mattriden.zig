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
            if (comptime meta.uploOf(Y) == .upper) {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        x_eff,
                        if (comptime noconj_o == noconj_y)
                            y.data[y._index(i, j)]
                        else
                            numeric.conj(y.data[y._index(i, j)]),
                    );
                }
            } else {
                while (i < int.min(j, o.rows)) : (i += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        x_eff,
                        numeric.cast(meta.Numeric(Y), 0),
                    );
                }
            }

            if (j < o.rows) {
                if (comptime meta.diagOf(Y) == .unit)
                    opInto(
                        &o.data[o._index(j, j)],
                        x_eff,
                        numeric.cast(meta.Numeric(Y), 1),
                    )
                else
                    opInto(
                        &o.data[o._index(j, j)],
                        x_eff,
                        if (comptime noconj_o == noconj_y)
                            y.data[y._index(j, j)]
                        else
                            numeric.conj(y.data[y._index(j, j)]),
                    );
            }

            i = int.min(j + 1, o.rows);
            if (comptime meta.uploOf(Y) == .lower) {
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
            } else {
                while (i < o.rows) : (i += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        x_eff,
                        numeric.cast(meta.Numeric(Y), 0),
                    );
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            if (comptime meta.uploOf(Y) == .lower) {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        x_eff,
                        if (comptime noconj_o == noconj_y)
                            y.data[y._index(i, j)]
                        else
                            numeric.conj(y.data[y._index(i, j)]),
                    );
                }
            } else {
                while (j < int.min(i, o.cols)) : (j += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        x_eff,
                        numeric.cast(meta.Numeric(Y), 0),
                    );
                }
            }

            if (i < o.cols) {
                if (comptime meta.diagOf(Y) == .unit)
                    opInto(
                        &o.data[o._index(i, i)],
                        x_eff,
                        numeric.cast(meta.Numeric(Y), 1),
                    )
                else
                    opInto(
                        &o.data[o._index(i, i)],
                        x_eff,
                        if (comptime noconj_o == noconj_y)
                            y.data[y._index(i, i)]
                        else
                            numeric.conj(y.data[y._index(i, i)]),
                    );
            }

            j = int.min(i + 1, o.cols);
            if (comptime meta.uploOf(Y) == .upper) {
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
            } else {
                while (j < o.cols) : (j += 1) {
                    opInto(
                        &o.data[o._index(i, j)],
                        x_eff,
                        numeric.cast(meta.Numeric(Y), 0),
                    );
                }
            }
        }
    }
}
