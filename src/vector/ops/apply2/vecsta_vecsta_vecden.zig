const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (y.flags.noconj) {
        true => k_apply2IntoUnchecked(o, x, y, opInto, true),
        false => k_apply2IntoUnchecked(o, x, y, opInto, false),
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_y: bool) void {
    if (y.inc == 1) {
        inline for (0..@TypeOf(x).len) |i| {
            opInto(
                &o.data[i],
                x.data[i],
                if (comptime noconj_y)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    } else {
        const iy: isize = if (y.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * y.inc else 0;
        inline for (0..@TypeOf(x).len) |i| {
            opInto(
                &o.data[i],
                x.data[i],
                if (comptime noconj_y)
                    y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]
                else
                    numeric.conj(y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]),
            );
        }
    }
}
