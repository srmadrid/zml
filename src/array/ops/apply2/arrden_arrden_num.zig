const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

const arrutils = @import("../../utils.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => switch (x.flags.noconj) {
            true => k_apply2IntoUnchecked(o, x, y, opInto, true, true),
            false => k_apply2IntoUnchecked(o, x, y, opInto, true, false),
        },
        false => switch (x.flags.noconj) {
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
    comptime noconj_x: bool,
) void {
    if ((o.isContiguous(.f) and x.isContiguous(.f)) or (o.isContiguous(.c) and x.isContiguous(.c))) {
        for (0..o._size()) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_o == noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
                y,
            );
        }

        return;
    }

    if (numeric.abs(o.strides[o.ndim - 1]) <= numeric.abs(o.strides[0]))
        stridedLoop(o, x, y, opInto, noconj_o, noconj_x, true, 0, 0, 0)
    else
        stridedLoop(o, x, y, opInto, noconj_o, noconj_x, false, 0, 0, 0);
}

fn stridedLoop(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
    comptime noconj_x: bool,
    comptime is_c_order: bool,
    o_offset: isize,
    x_offset: isize,
    level: usize,
) void {
    if (level == o.ndim) {
        const o_ptr = arrutils.getPtr(o.data, o_offset);
        const x_ptr = arrutils.getPtr(x.data, x_offset);

        opInto(
            &o_ptr[0],
            if (comptime noconj_o == noconj_x)
                x_ptr[0]
            else
                numeric.conj(x_ptr[0]),
            y,
        );

        return;
    }

    const dim = if (comptime is_c_order) level else o.ndim - 1 - level;

    var curr_o = o_offset;
    var curr_x = x_offset;
    for (0..o.shape[dim]) |_| {
        stridedLoop(o, x, y, opInto, noconj_o, noconj_x, is_c_order, curr_o, curr_x, level + 1);

        curr_o += o.strides[dim];
        curr_x += x.strides[dim];
    }
}
