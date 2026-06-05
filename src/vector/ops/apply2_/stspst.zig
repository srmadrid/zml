const std = @import("std");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const Y: type = @TypeOf(y);

    const aliased = (comptime O == Y) and std.meta.eql(o.*, y);

    if ((comptime op_ == numeric.add_) and aliased) {
        var ix: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            op_(&o.data[x.idx[ix]], x.data[ix], y.data[x.idx[ix]]);
        }
    } else {
        var i: usize = 0;
        var ix: usize = 0;
        while (ix < x.nnz) : (ix += 1) {
            while (i < x.idx[ix]) : (i += 1) {
                if (comptime op_ == numeric.add_)
                    numeric.set(&o.data[i], y.data[i])
                else
                    numeric.set(&o.data[i], numeric.neg(y.data[i]));
            }

            op_(&o.data[i], x.data[ix], y.data[i]);

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            if (comptime op_ == numeric.add_)
                numeric.set(&o.data[i], y.data[i])
            else
                numeric.set(&o.data[i], numeric.neg(y.data[i]));
        }
    }
}
