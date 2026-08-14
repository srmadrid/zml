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
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    if (o.inc == 1) {
        var ix: usize = 0;
        var iy: usize = 0;
        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            if (ix < x.nnz and x.idx[ix] == i) {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(
                        &o.data[i],
                        if (comptime noconj_o == noconj_x)
                            x.data[ix]
                        else
                            numeric.conj(x.data[ix]),
                        if (comptime noconj_o == noconj_y)
                            y.data[iy]
                        else
                            numeric.conj(y.data[iy]),
                    );

                    ix += 1;
                    iy += 1;
                } else {
                    opInto(
                        &o.data[i],
                        if (comptime noconj_o == noconj_x)
                            x.data[ix]
                        else
                            numeric.conj(x.data[ix]),
                        numeric.cast(meta.Numeric(Y), 0),
                    );

                    ix += 1;
                }
            } else {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(
                        &o.data[i],
                        numeric.cast(meta.Numeric(X), 0),
                        if (comptime noconj_o == noconj_y)
                            y.data[iy]
                        else
                            numeric.conj(y.data[iy]),
                    );

                    iy += 1;
                } else {
                    opInto(&o.data[i], numeric.cast(meta.Numeric(X), 0), numeric.cast(meta.Numeric(Y), 0));
                }
            }
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: usize = 0;
        var iy: usize = 0;
        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            if (ix < x.nnz and x.idx[ix] == i) {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(
                        &o.data[numeric.cast(usize, io)],
                        if (comptime noconj_o == noconj_x)
                            x.data[ix]
                        else
                            numeric.conj(x.data[ix]),
                        if (comptime noconj_o == noconj_y)
                            y.data[iy]
                        else
                            numeric.conj(y.data[iy]),
                    );

                    ix += 1;
                    iy += 1;
                } else {
                    opInto(
                        &o.data[numeric.cast(usize, io)],
                        if (comptime noconj_o == noconj_x)
                            x.data[ix]
                        else
                            numeric.conj(x.data[ix]),
                        numeric.cast(meta.Numeric(Y), 0),
                    );

                    ix += 1;
                }
            } else {
                if (iy < y.nnz and y.idx[iy] == i) {
                    opInto(
                        &o.data[numeric.cast(usize, io)],
                        numeric.cast(meta.Numeric(X), 0),
                        if (comptime noconj_o == noconj_y)
                            y.data[iy]
                        else
                            numeric.conj(y.data[iy]),
                    );

                    iy += 1;
                } else {
                    opInto(&o.data[numeric.cast(usize, io)], numeric.cast(meta.Numeric(X), 0), numeric.cast(meta.Numeric(Y), 0));
                }
            }

            io += o.inc;
        }
    }
}
