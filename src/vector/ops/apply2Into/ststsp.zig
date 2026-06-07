const std = @import("std");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);

    const aliased = (comptime O == X) and std.meta.eql(o.*, x);

    if (aliased) {
        var iy: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            opInto(&o.data[y.idx[iy]], x.data[y.idx[iy]], y.data[iy]);
        }
    } else {
        var i: usize = 0;
        var iy: usize = 0;
        while (iy < y.nnz) : (iy += 1) {
            while (i < y.idx[iy]) : (i += 1) {
                numeric.set(&o.data[i], x.data[i]);
            }

            opInto(&o.data[i], x.data[i], y.data[iy]);

            i += 1;
        }

        while (i < meta.Child(@TypeOf(o)).len) : (i += 1) {
            numeric.set(&o.data[i], x.data[i]);
        }
    }
}
