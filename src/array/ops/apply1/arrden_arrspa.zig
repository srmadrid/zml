const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

const arrutils = @import("../../utils.zig");

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
    const X = @TypeOf(x);

    if (o.isContiguous(.c) or o.isContiguous(.f)) {
        for (0..o._size()) |i| {
            opInto(&o.data[i], numeric.cast(meta.Numeric(X), 0));
        }
    } else {
        if (numeric.abs(o.strides[o.ndim - 1]) <= numeric.abs(o.strides[0]))
            fillZeroStrided(o, opInto, numeric.cast(meta.Numeric(X), 0), true, 0, 0)
        else
            fillZeroStrided(o, opInto, numeric.cast(meta.Numeric(X), 0), false, 0, 0);
    }

    if (x.nnz == 0)
        return;

    if (x.ndim == 0) {
        opInto(&o.data[0], if (comptime noconj_o == noconj_x) x.data[0] else numeric.conj(x.data[0]));

        return;
    }

    scatterSparse(o, x, opInto, noconj_o, noconj_x, 0, 0, x._ilen[0], 0);
}

fn fillZeroStrided(
    o: anytype,
    comptime opInto: anytype,
    zero_x: anytype,
    comptime is_c_order: bool,
    o_offset: isize,
    level: usize,
) void {
    if (level == o.ndim) {
        const o_ptr = arrutils.getPtr(o.data, o_offset);

        opInto(&o_ptr[0], zero_x);

        return;
    }

    const dim = if (comptime is_c_order) level else o.ndim - 1 - level;

    var curr_o = o_offset;
    for (0..o.shape[dim]) |_| {
        fillZeroStrided(o, opInto, zero_x, is_c_order, curr_o, level + 1);

        curr_o += o.strides[dim];
    }
}

fn scatterSparse(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
    comptime noconj_x: bool,
    o_offset: isize,
    range_start: usize,
    range_end: usize,
    level: usize,
) void {
    const X = @TypeOf(x);

    const dim = if (comptime X.storage_order == .c) level else x.ndim - 1 - level;

    var i: usize = range_start;
    while (i < range_end) : (i += 1) {
        if (o.shape[dim] > 1 and x.shape[dim] == 1) {
            var curr_o = o_offset;

            for (0..o.shape[dim]) |_| {
                if (level == x.ndim - 1) {
                    const o_ptr = arrutils.getPtr(o.data, curr_o);

                    opInto(&o_ptr[0], if (comptime noconj_o == noconj_x) x.data[i] else numeric.conj(x.data[i]));
                } else {
                    scatterSparse(o, x, opInto, noconj_o, noconj_x, curr_o, x.ptr[level][i], x.ptr[level][i + 1], level + 1);
                }

                curr_o += o.strides[dim];
            }
        } else {
            const curr_o = o_offset + numeric.cast(isize, x.idx[level][i]) * o.strides[dim];

            if (level == x.ndim - 1) {
                const o_ptr = arrutils.getPtr(o.data, curr_o);

                opInto(&o_ptr[0], if (comptime noconj_o == noconj_x) x.data[i] else numeric.conj(x.data[i]));
            } else {
                scatterSparse(o, x, opInto, noconj_o, noconj_x, curr_o, x.ptr[level][i], x.ptr[level][i + 1], level + 1);
            }
        }
    }
}
