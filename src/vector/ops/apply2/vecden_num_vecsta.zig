const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => k_apply2IntoUnchecked(o, x, y, opInto, true),
        false => k_apply2IntoUnchecked(o, x, y, opInto, false),
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_o: bool) void {
    const x_eff = if (comptime noconj_o) x else numeric.conj(x);

    if (o.inc == 1) {
        inline for (0..@TypeOf(y).len) |i| {
            opInto(
                &o.data[i],
                x_eff,
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    } else {
        const io: isize = if (o.inc < 0) (-numeric.cast(isize, @TypeOf(y).len) + 1) * o.inc else 0;
        inline for (0..@TypeOf(y).len) |i| {
            opInto(
                &o.data[numeric.cast(usize, io + numeric.cast(isize, i) * o.inc)],
                x_eff,
                if (comptime noconj_o)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    }
}
