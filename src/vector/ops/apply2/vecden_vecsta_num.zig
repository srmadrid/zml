const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => k_apply2IntoUnchecked(o, x, y, opInto, true),
        false => k_apply2IntoUnchecked(o, x, y, opInto, false),
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_o: bool) void {
    const y_eff = if (comptime noconj_o) y else numeric.conj(y);

    if (o.inc == 1) {
        inline for (0..@TypeOf(x).len) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                y_eff,
            );
        }
    } else {
        const io: isize = if (o.inc < 0) (-numeric.cast(isize, @TypeOf(x).len) + 1) * o.inc else 0;
        inline for (0..@TypeOf(x).len) |i| {
            opInto(
                &o.data[numeric.cast(usize, io + numeric.cast(isize, i) * o.inc)],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                y_eff,
            );
        }
    }
}
