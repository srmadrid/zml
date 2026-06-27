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
    if (o.inc == 1 and x.inc == 1) {
        inline for (0..@TypeOf(y).len) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_o == noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    } else {
        const io: isize = if (o.inc < 0) (-numeric.cast(isize, @TypeOf(y).len) + 1) * o.inc else 0;
        const ix: isize = if (x.inc < 0) (-numeric.cast(isize, @TypeOf(y).len) + 1) * x.inc else 0;
        inline for (0..@TypeOf(y).len) |i| {
            opInto(
                &o.data[numeric.cast(usize, io + numeric.cast(isize, i) * o.inc)],
                if (comptime noconj_o == noconj_x)
                    x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)]
                else
                    numeric.conj(x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)]),
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    }
}
