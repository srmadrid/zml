const std = @import("std");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const aliased = (comptime O == Y) and std.meta.eql(o.*, y);

    if ((comptime opInto == numeric.addInto) and aliased) {
        var ix: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            opInto(&o.data[x.idx[ix]], x.data[ix], y.data[x.idx[ix]]);
        }
    } else {
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
}
