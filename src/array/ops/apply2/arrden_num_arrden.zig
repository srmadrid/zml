const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

const arrutils = @import("../../utils.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => switch (y.flags.noconj) {
            true => k_apply2IntoUnchecked(o, x, y, opInto, true, true),
            false => k_apply2IntoUnchecked(o, x, y, opInto, true, false),
        },
        false => switch (y.flags.noconj) {
            true => k_apply2IntoUnchecked(o, x, y, opInto, false, true),
            false => k_apply2IntoUnchecked(o, x, y, opInto, false, false),
        },
    };
}

fn k_apply2IntoUnchecked(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
    comptime noconj_y: bool,
) void {
    if ((o.isContiguous(.f) and y.isContiguous(.f)) or (o.isContiguous(.c) and y.isContiguous(.c))) {
        for (0..o._size()) |i| {
            opInto(
                &o.data[i],
                x,
                if (comptime noconj_o == noconj_y)
                    y.data[i]
                else
                    numeric.conj(y.data[i]),
            );
        }

        return;
    }

    if (numeric.abs(o.strides[o.ndim - 1]) <= numeric.abs(o.strides[0]))
        stridedLoop(o, x, y, opInto, noconj_o, noconj_y, true, 0, 0, 0)
    else
        stridedLoop(o, x, y, opInto, noconj_o, noconj_y, false, 0, 0, 0);
}

fn stridedLoop(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
    comptime noconj_y: bool,
    comptime is_c_order: bool,
    o_offset: isize,
    y_offset: isize,
    level: usize,
) void {
    if (level == o.ndim) {
        const o_ptr = arrutils.getPtr(o.data, o_offset);
        const y_ptr = arrutils.getPtr(y.data, y_offset);

        opInto(
            &o_ptr[0],
            x,
            if (comptime noconj_o == noconj_y)
                y_ptr[0]
            else
                numeric.conj(y_ptr[0]),
        );

        return;
    }

    const dim = if (comptime is_c_order) level else o.ndim - 1 - level;

    var curr_o = o_offset;
    var curr_y = y_offset;
    for (0..o.shape[dim]) |_| {
        stridedLoop(o, x, y, opInto, noconj_o, noconj_y, is_c_order, curr_o, curr_y, level + 1);

        curr_o += o.strides[dim];
        curr_y += y.strides[dim];
    }
}
