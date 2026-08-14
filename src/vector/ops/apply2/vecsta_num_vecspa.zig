const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (y.flags.noconj) {
        true => k_apply2IntoUnchecked(o, x, y, opInto, true),
        false => k_apply2IntoUnchecked(o, x, y, opInto, false),
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_y: bool) void {
    const O = meta.Child(@TypeOf(o));
    const Y = @TypeOf(y);

    var i: usize = 0;
    var iy: usize = 0;
    while (iy < y.nnz) : (iy += 1) {
        while (i < y.idx[iy]) : (i += 1) {
            opInto(&o.data[i], x, numeric.cast(meta.Numeric(Y), 0));
        }

        opInto(
            &o.data[i],
            x,
            if (comptime noconj_y)
                y.data[iy]
            else
                numeric.conj(y.data[iy]),
        );

        i += 1;
    }

    while (i < O.len) : (i += 1) {
        opInto(&o.data[i], x, numeric.cast(meta.Numeric(Y), 0));
    }
}
