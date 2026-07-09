const std = @import("std");

const meta = @import("../../../meta.zig");

pub fn apply1IntoUnchecked(o: anytype, x: anytype, comptime opInto: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    if (comptime O.storage_order == X.storage_order and std.mem.eql(usize, O.shape[0..O.ndim], X.shape[0..X.ndim])) {
        for (0..O.size) |i| {
            opInto(&o.data[i], x.data[i]);
        }

        return;
    }

    k_apply1IntoUnchecked(o, x, opInto, 0, 0, 0);
}

fn k_apply1IntoUnchecked(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    o_offset: usize,
    x_offset: usize,
    comptime level: usize,
) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    if (comptime level == O.ndim) {
        opInto(&o.data[o_offset], x.data[x_offset]);

        return;
    }

    const dim = if (comptime O.storage_order == .c) level else O.ndim - 1 - level;

    var curr_o = o_offset;
    var curr_x = x_offset;
    for (0..O.shape[dim]) |_| {
        k_apply1IntoUnchecked(o, x, opInto, curr_o, curr_x, level + 1);

        curr_o += O.strides[dim];
        curr_x += X.strides[dim];
    }
}
