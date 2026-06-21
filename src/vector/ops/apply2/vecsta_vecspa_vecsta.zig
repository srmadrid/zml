const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X = @TypeOf(x);

    var i: usize = 0;
    var ix: usize = 0;
    while (ix < x.nnz) : (ix += 1) {
        while (i < x.idx[ix]) : (i += 1) {
            opInto(&o.data[i], numeric.zero(meta.Numeric(X)), y.data[i]);
        }

        opInto(&o.data[i], x.data[ix], y.data[i]);

        i += 1;
    }

    while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
        opInto(&o.data[i], numeric.zero(meta.Numeric(X)), y.data[i]);
    }
}
