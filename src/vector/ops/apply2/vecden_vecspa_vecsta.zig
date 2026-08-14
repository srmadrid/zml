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
    const X = @TypeOf(x);

    if (o.inc == 1) {
        var ix: usize = 0;
        var i: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(
                    &o.data[i],
                    numeric.cast(meta.Numeric(X), 0),
                    if (comptime noconj_o)
                        y.data[i]
                    else
                        numeric.conj(y.data[i]),
                );
            }

            opInto(
                &o.data[i],
                if (comptime noconj_o == noconj_x)
                    x.data[ix]
                else
                    numeric.conj(x.data[ix]),
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );

            i += 1;
        }

        while (i < o.len) : (i += 1) {
            opInto(
                &o.data[i],
                numeric.cast(meta.Numeric(X), 0),
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: usize = 0;
        var i: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                opInto(
                    &o.data[numeric.cast(usize, io)],
                    numeric.cast(meta.Numeric(X), 0),
                    if (comptime noconj_o)
                        y.data[i]
                    else
                        numeric.conj(y.data[i]),
                );

                io += o.inc;
            }

            opInto(
                &o.data[numeric.cast(usize, io)],
                if (comptime noconj_o == noconj_x)
                    x.data[ix]
                else
                    numeric.conj(x.data[ix]),
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );

            i += 1;
            io += o.inc;
        }

        while (i < o.len) : (i += 1) {
            opInto(
                &o.data[numeric.cast(usize, io)],
                numeric.cast(meta.Numeric(X), 0),
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );

            io += o.inc;
        }
    }
}
