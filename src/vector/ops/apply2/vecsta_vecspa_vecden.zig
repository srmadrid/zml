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

    if (y.inc == 1) {
        var ix: usize = 0;
        var i: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(
                    &o.data[i],
                    numeric.cast(meta.Numeric(X), 0),
                    if (comptime noconj_y)
                        y.data[i]
                    else
                        numeric.conj(y.data[i]),
                );
            }

            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[ix]
                else
                    numeric.conj(x.data[ix]),
                if (comptime noconj_y)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            opInto(
                &o.data[i],
                numeric.cast(meta.Numeric(X), 0),
                if (comptime noconj_y)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    } else {
        var ix: usize = 0;
        var iy: isize = if (y.inc < 0) (-numeric.cast(isize, y.len) + 1) * y.inc else 0;
        var i: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(
                    &o.data[i],
                    numeric.cast(meta.Numeric(X), 0),
                    if (comptime noconj_y)
                        y.data[numeric.cast(usize, iy)]
                    else
                        numeric.conj(y.data[numeric.cast(usize, iy)]),
                );

                iy += y.inc;
            }

            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[ix]
                else
                    numeric.conj(x.data[ix]),
                if (comptime noconj_y)
                    y.data[numeric.cast(usize, iy)]
                else
                    numeric.conj(y.data[numeric.cast(usize, iy)]),
            );

            i += 1;
            iy += y.inc;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            opInto(
                &o.data[i],
                numeric.cast(meta.Numeric(X), 0),
                if (comptime noconj_y)
                    y.data[numeric.cast(usize, iy)]
                else
                    numeric.conj(y.data[numeric.cast(usize, iy)]),
            );

            iy += y.inc;
        }
    }
}
