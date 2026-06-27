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
    if (o.inc == 1 and x.inc == 1 and y.inc == 1) {
        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            opInto(
                &o.data[i],
                if (comptime noconj_o == noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                if (comptime noconj_o == noconj_y)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    } else {
        var io: isize = if (o.inc < 0) (-numeric.cast(isize, o.len) + 1) * o.inc else 0;
        var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
        var iy: isize = if (y.inc < 0) (-numeric.cast(isize, y.len) + 1) * y.inc else 0;
        var i: usize = 0;
        while (i < o.len) : (i += 1) {
            opInto(
                &o.data[numeric.cast(usize, io)],
                if (comptime noconj_o == noconj_x)
                    x.data[numeric.cast(usize, ix)]
                else
                    numeric.conj(x.data[numeric.cast(usize, ix)]),
                if (comptime noconj_o == noconj_y)
                    y.data[numeric.cast(usize, iy)]
                else
                    numeric.conj(y.data[numeric.cast(usize, iy)]),
            );

            io += o.inc;
            ix += x.inc;
            iy += y.inc;
        }
    }
}
