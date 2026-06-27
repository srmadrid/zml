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
    if (x.inc == 1 and y.inc == 1) {
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                if (comptime noconj_y)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }
    } else {
        const ix: isize = if (x.inc < 0) (-numeric.cast(isize, meta.Child(@TypeOf(o)).len) + 1) * x.inc else 0;
        const iy: isize = if (y.inc < 0) (-numeric.cast(isize, meta.Child(@TypeOf(o)).len) + 1) * y.inc else 0;
        inline for (0..meta.Child(@TypeOf(o)).len) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)]
                else
                    numeric.conj(x.data[numeric.cast(usize, ix + numeric.cast(isize, i) * x.inc)]),
                if (comptime noconj_y)
                    y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]
                else
                    numeric.conj(y.data[numeric.cast(usize, iy + numeric.cast(isize, i) * y.inc)]),
            );
        }
    }
}
