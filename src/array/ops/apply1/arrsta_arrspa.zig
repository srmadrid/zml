const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply1IntoUnchecked(o: anytype, x: anytype, comptime opInto: anytype) void {
    return switch (x.flags.noconj) {
        true => k_apply1IntoUnchecked(o, x, opInto, true),
        false => k_apply1IntoUnchecked(o, x, opInto, false),
    };
}

fn k_apply1IntoUnchecked(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_x: bool,
) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    for (0..O.size) |i| {
        opInto(&o.data[i], numeric.cast(meta.Numeric(X), 0));
    }

    if (x.nnz == 0)
        return;

    if (x.ndim == 0) {
        opInto(&o.data[0], if (comptime noconj_x) x.data[0] else numeric.conj(x.data[0]));

        return;
    }

    scatterSparse(o, x, opInto, noconj_x, 0, 0, x._ilen[0], 0);
}

fn scatterSparse(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_x: bool,
    o_offset: usize,
    range_start: usize,
    range_end: usize,
    comptime level: usize,
) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    const dim = if (comptime X.storage_order == .c) level else x.ndim - 1 - level;

    var i: usize = range_start;
    while (i < range_end) : (i += 1) {
        if ((comptime O.shape[dim] > 1) and x.shape[dim] == 1) {
            var curr_o = o_offset;
            for (0..O.shape[dim]) |_| {
                if (comptime level == x.ndim - 1)
                    opInto(&o.data[curr_o], if (comptime noconj_x) x.data[i] else numeric.conj(x.data[i]))
                else
                    scatterSparse(o, x, opInto, noconj_x, curr_o, x.ptr[level][i], x.ptr[level][i + 1], level + 1);

                curr_o += O.strides[dim];
            }
        } else {
            const curr_o = o_offset + x.idx[level][i] * O.strides[dim];

            if (level == x.ndim - 1)
                opInto(&o.data[curr_o], if (comptime noconj_x) x.data[i] else numeric.conj(x.data[i]))
            else
                scatterSparse(o, x, opInto, noconj_x, curr_o, x.ptr[level][i], x.ptr[level][i + 1], level + 1);
        }
    }
}
