const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));

    var i: usize = 0;
    var iy: usize = 0;
    while (iy < y.nnz) : (iy += 1) {
        while (i < y.idx[iy]) : (i += 1) {
            o.data[i] = numeric.zero(meta.Numeric(O));
        }

        opInto(&o.data[i], x, y.data[iy]);

        i += 1;
    }

    while (i < O.len) : (i += 1) {
        o.data[i] = numeric.zero(meta.Numeric(O));
    }
}
