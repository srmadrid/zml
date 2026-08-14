const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (x.flags.noconj) {
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

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_x: bool, comptime noconj_y: bool) void {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    var i: usize = 0;
    var ix: usize = 0;
    var iy: usize = 0;
    while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
        if (ix < x.nnz and x.idx[ix] == i) {
            if (iy < y.nnz and y.idx[iy] == i) {
                opInto(
                    &o.data[i],
                    if (comptime noconj_x)
                        x.data[ix]
                    else
                        numeric.conj(x.data[ix]),
                    if (comptime noconj_y)
                        y.data[iy]
                    else
                        numeric.conj(y.data[iy]),
                );

                ix += 1;
                iy += 1;
            } else {
                opInto(
                    &o.data[i],
                    if (comptime noconj_x)
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
                    if (comptime noconj_y)
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
}
