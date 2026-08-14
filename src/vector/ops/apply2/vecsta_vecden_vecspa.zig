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
    const Y = @TypeOf(y);

    if (x.inc == 1) {
        var iy: usize = 0;
        var i: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(
                    &o.data[i],
                    if (comptime noconj_x)
                        x.data[i]
                    else
                        numeric.conj(x.data[i]),
                    numeric.cast(meta.Numeric(Y), 0),
                );
            }

            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                if (comptime noconj_y)
                    y.data[iy]
                else
                    numeric.conj(y.data[iy]),
            );

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                numeric.cast(meta.Numeric(Y), 0),
            );
        }
    } else {
        var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
        var iy: usize = 0;
        var i: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(
                    &o.data[i],
                    if (comptime noconj_x)
                        x.data[numeric.cast(usize, ix)]
                    else
                        numeric.conj(x.data[numeric.cast(usize, ix)]),
                    numeric.cast(meta.Numeric(Y), 0),
                );

                ix += x.inc;
            }

            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[numeric.cast(usize, ix)]
                else
                    numeric.conj(x.data[numeric.cast(usize, ix)]),
                if (comptime noconj_y)
                    y.data[iy]
                else
                    numeric.conj(y.data[iy]),
            );

            i += 1;
            ix += x.inc;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[numeric.cast(usize, ix)]
                else
                    numeric.conj(x.data[numeric.cast(usize, ix)]),
                numeric.cast(meta.Numeric(Y), 0),
            );

            ix += x.inc;
        }
    }
}
