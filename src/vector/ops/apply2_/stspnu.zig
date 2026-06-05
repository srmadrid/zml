const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));

    var i: usize = 0;
    var ix: usize = 0;
    while (ix < x.nnz) : (ix += 1) {
        while (i < x.idx[ix]) : (i += 1) {
            o.data[i] = numeric.zero(meta.Numeric(O));
        }

        op_(&o.data[i], x.data[ix], y);

        i += 1;
    }

    while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
        o.data[i] = numeric.zero(meta.Numeric(O));
    }
}
