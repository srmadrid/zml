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
    if (o.inc == 1 and y.inc == 1) {
        inline for (0..@TypeOf(x).len) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_o)
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
        const io: isize = if (o.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * o.inc else 0;
        const iy: isize = if (y.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * y.inc else 0;
        inline for (0..@TypeOf(x).len) |i| {
            opInto(
                &o.data[numeric.cast(usize, io + numeric.cast(isize, i) * o.inc)],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                if (comptime noconj_o == noconj_y)
                    y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]
                else
                    numeric.conj(y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]),
            );
        }
    }
}
