const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const array = @import("../../../array.zig");

pub fn apply1IntoUnchecked(o: anytype, x: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => switch (x.flags.noconj) {
            true => k_apply1IntoUnchecked(o, x, opInto, true, true),
            false => k_apply1IntoUnchecked(o, x, opInto, true, false),
        },
        false => switch (x.flags.noconj) {
            true => k_apply1IntoUnchecked(o, x, opInto, false, true),
            false => k_apply1IntoUnchecked(o, x, opInto, false, false),
        },
    };
}

fn k_apply1IntoUnchecked(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
    comptime noconj_x: bool,
) void {
    if (x.ndim == 0) {
        o.nnz = x.nnz;

        if (x.nnz > 0) {
            opInto(&o.data[0], if (comptime noconj_o == noconj_x) x.data[0] else numeric.conj(x.data[0]));
        }

        return;
    }

    if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
        o.nnz = x.nnz;

        if (o.data != x.data) {
            var x_nodes_at_level: usize = 1;
            for (0..x.ndim) |level| {
                const parent_nodes = x_nodes_at_level;
                x_nodes_at_level = if (level == x.ndim - 1) x.nnz else x.ptr[level][parent_nodes];

                @memcpy(o.idx[level][0..x_nodes_at_level], x.idx[level][0..x_nodes_at_level]);

                if (level > 0) {
                    @memcpy(o.ptr[level][0 .. parent_nodes + 1], x.ptr[level][0 .. parent_nodes + 1]);
                } else {
                    @memcpy(o.ptr[0][0..2], x.ptr[0][0..2]);
                }
            }
        }

        for (0..x.nnz) |i| {
            opInto(&o.data[i], if (comptime noconj_o == noconj_x) x.data[i] else numeric.conj(x.data[i]));
        }

        return;
    }

    var state = BroadcastState{
        .idx_pos = .{0} ** array.max_dimensions,
        .data_pos = 0,
    };

    broadcastSparseLevel(o, x, opInto, noconj_o, noconj_x, 0, 0, 0, &state);

    o.nnz = state.data_pos;
}

const BroadcastState = struct {
    idx_pos: [array.max_dimensions]usize,
    data_pos: usize,
};

fn broadcastSparseLevel(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
    comptime noconj_x: bool,
    level: usize,
    x_parent_node: usize,
    o_parent_node: usize,
    state: *BroadcastState,
) void {
    const X = @TypeOf(x);

    const dim = if (comptime X.storage_order == .c) level else x.ndim - 1 - level;

    const x_child_start = if (level == 0) 0 else x.ptr[level][x_parent_node];
    const x_child_end = if (level == 0) x.ptr[0][1] else x.ptr[level][x_parent_node + 1];

    if (level > 0)
        o.ptr[level][o_parent_node] = state.idx_pos[level]
    else
        o.ptr[0][0] = 0;

    if ((o.shape[dim] > 1 and x.shape[dim] == 1)) {
        if (x_child_start < x_child_end) {
            const x_child_node = x_child_start;

            for (0..o.shape[dim]) |new_coord| {
                const current_o_node = state.idx_pos[level];

                o.idx[level][current_o_node] = new_coord;
                state.idx_pos[level] += 1;

                if (level == x.ndim - 1) {
                    opInto(&o.data[state.data_pos], if (comptime noconj_o == noconj_x) x.data[x_child_node] else numeric.conj(x.data[x_child_node]));

                    state.data_pos += 1;
                } else {
                    broadcastSparseLevel(o, x, opInto, noconj_o, noconj_x, level + 1, x_child_node, current_o_node, state);
                }
            }
        }
    } else {
        var curr_x_child = x_child_start;
        while (curr_x_child < x_child_end) : (curr_x_child += 1) {
            const current_o_node = state.idx_pos[level];

            o.idx[level][current_o_node] = x.idx[level][curr_x_child];
            state.idx_pos[level] += 1;

            if (level == x.ndim - 1) {
                opInto(&o.data[state.data_pos], if (comptime noconj_o == noconj_x) x.data[curr_x_child] else numeric.conj(x.data[curr_x_child]));

                state.data_pos += 1;
            } else {
                broadcastSparseLevel(o, x, opInto, noconj_o, noconj_x, level + 1, curr_x_child, current_o_node, state);
            }
        }
    }

    if (level > 0)
        o.ptr[level][o_parent_node + 1] = state.idx_pos[level]
    else
        o.ptr[0][1] = state.idx_pos[0];
}
