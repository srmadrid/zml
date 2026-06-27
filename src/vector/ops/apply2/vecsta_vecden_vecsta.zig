const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (x.flags.noconj) {
        true => k_apply2IntoUnchecked(o, x, y, opInto, true),
        false => k_apply2IntoUnchecked(o, x, y, opInto, false),
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_x: bool) void {
    if (x.inc == 1) {
        inline for (0..@TypeOf(y).len) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                y.data[i],
            );
        }
    } else {
        const ix: isize = if (x.inc < 0) (-numeric.cast(isize, @TypeOf(y).len) + 1) * x.inc else 0;
        inline for (0..@TypeOf(y).len) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)]
                else
                    numeric.conj(x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)]),
                y.data[i],
            );
        }
    }
}
