const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const array = @import("../../../array.zig");

pub fn apply1IntoUnchecked(o: anytype, x: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => k_apply1IntoUnchecked(o, x, opInto, true),
        false => k_apply1IntoUnchecked(o, x, opInto, false),
    };
}

fn k_apply1IntoUnchecked(o: anytype, x: anytype, comptime opInto: anytype, comptime noconj_x: bool) void {
    const O: type = meta.Child(@TypeOf(o));

    if (x.ndim == 0) {
        if (x.nnz > 0) {
            var out_val: meta.Numeric(O) = undefined;
            opInto(&out_val, if (comptime noconj_x) x.data[0] else numeric.conj(x.data[0]));
            o.appendAssumeCapacity(&[_]usize{}, out_val);
        }

        return;
    }

    var current_idx: [array.max_dimensions]usize = undefined;
    broadcastSparseLevel(o, x, opInto, noconj_x, 0, 0, &current_idx);
}

fn broadcastSparseLevel(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_x: bool,
    level: usize,
    node: usize,
    current_idx: *[array.max_dimensions]usize,
) void {
    const O: type = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    const dim = if (comptime X.storage_order == .c) level else x.ndim - 1 - level;

    const child_start = if (level == 0) 0 else x.ptr[level][node];
    const child_end = if (level == 0) x.ptr[0][1] else x.ptr[level][node + 1];

    if (o.shape[dim] > 1 and x.shape[dim] == 1) {
        if (child_start < child_end) {
            const child_node = child_start;

            for (0..o.shape[dim]) |new_coord| {
                current_idx[dim] = new_coord;

                if (level == x.ndim - 1) {
                    var out_val: meta.Numeric(O) = undefined;
                    opInto(&out_val, if (comptime noconj_x) x.data[child_node] else numeric.conj(x.data[child_node]));
                    o.appendAssumeCapacity(current_idx[0..o.ndim], out_val);
                } else {
                    broadcastSparseLevel(o, x, opInto, noconj_x, level + 1, child_node, current_idx);
                }
            }
        }
    } else {
        var child = child_start;
        while (child < child_end) : (child += 1) {
            current_idx[dim] = x.idx[level][child];

            if (level == x.ndim - 1) {
                var out_val: meta.Numeric(O) = undefined;
                opInto(&out_val, if (comptime noconj_x) x.data[child] else numeric.conj(x.data[child]));
                o.appendAssumeCapacity(current_idx[0..o.ndim], out_val);
            } else {
                broadcastSparseLevel(o, x, opInto, noconj_x, level + 1, child, current_idx);
            }
        }
    }
}
