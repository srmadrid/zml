const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const Y = @TypeOf(y);

    if (x.inc == 1) {
        var iy: usize = 0;
        var i: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(&o.data[i], x.data[i], numeric.zero(meta.Numeric(Y)));
            }

            opInto(&o.data[i], x.data[i], y.data[iy]);

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            opInto(&o.data[i], x.data[i], numeric.zero(meta.Numeric(Y)));
        }
    } else {
        var ix: isize = if (x.inc < 0) (-numeric.cast(isize, x.len) + 1) * x.inc else 0;
        var iy: usize = 0;
        var i: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                opInto(&o.data[i], x.data[numeric.cast(usize, ix)], numeric.zero(meta.Numeric(Y)));

                ix += x.inc;
            }

            opInto(&o.data[i], x.data[numeric.cast(usize, ix)], y.data[iy]);

            i += 1;
            ix += x.inc;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            opInto(&o.data[i], x.data[numeric.cast(usize, ix)], numeric.zero(meta.Numeric(Y)));

            ix += x.inc;
        }
    }
}
