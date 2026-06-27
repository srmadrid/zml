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
    const Y = @TypeOf(y);

    if (o.inc == 1) {
        var iy: usize = 0;
        var i: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(
                    &o.data[i],
                    if (comptime noconj_o)
                        x.data[i]
                    else
                        numeric.conj(x.data[i]),
                    numeric.zero(meta.Numeric(Y)),
                );
            }

            opInto(
                &o.data[i],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                if (comptime noconj_o == noconj_y)
                    y.data[iy]
                else
                    numeric.conj(y.data[iy]),
            );

            i += 1;
        }

        while (i < o.len) : (i += 1) {
            opInto(
                &o.data[i],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                numeric.zero(meta.Numeric(Y)),
            );
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var iy: usize = 0;
        var i: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(
                    &o.data[numeric.cast(usize, io)],
                    if (comptime noconj_o)
                        x.data[i]
                    else
                        numeric.conj(x.data[i]),
                    numeric.zero(meta.Numeric(Y)),
                );

                io += o.inc;
            }

            opInto(
                &o.data[numeric.cast(usize, io)],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                if (comptime noconj_o == noconj_y)
                    y.data[iy]
                else
                    numeric.conj(y.data[iy]),
            );

            i += 1;
            io += o.inc;
        }

        while (i < o.len) : (i += 1) {
            opInto(
                &o.data[numeric.cast(usize, io)],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                numeric.zero(meta.Numeric(Y)),
            );

            io += o.inc;
        }
    }
}
