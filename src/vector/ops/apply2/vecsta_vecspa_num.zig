const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (x.flags.noconj) {
        true => k_apply2IntoUnchecked(o, x, y, opInto, true),
        false => k_apply2IntoUnchecked(o, x, y, opInto, false),
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_x: bool) void {
    const X = @TypeOf(x);

    var i: usize = 0;
    var ix: usize = 0;
    while (ix < x.nnz) : (ix += 1) {
        while (i < x.idx[ix]) : (i += 1) {
            opInto(&o.data[i], numeric.cast(meta.Numeric(X), 0), y);
        }

        opInto(
            &o.data[i],
            if (comptime noconj_x)
                x.data[ix]
            else
                numeric.conj(x.data[ix]),
            y,
        );

        i += 1;
    }

    while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
        opInto(&o.data[i], numeric.cast(meta.Numeric(X), 0), y);
    }
}
