const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

const arrutils = @import("../../utils.zig");

pub fn apply1IntoUnchecked(o: anytype, x: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => k_apply1IntoUnchecked(o, x, opInto, true),
        false => k_apply1IntoUnchecked(o, x, opInto, false),
    };
}

fn k_apply1IntoUnchecked(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
) void {
    const X = @TypeOf(x);

    if (((comptime X.storage_order == .f) and o.isContiguous(.f)) or ((comptime X.storage_order == .c) and o.isContiguous(.c))) {
        for (0..o._size()) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_o)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
            );
        }

        return;
    }

    if (numeric.abs(o.strides[o.ndim - 1]) <= numeric.abs(o.strides[0]))
        stridedLoop(o, x, opInto, noconj_o, true, 0, 0, 0)
    else
        stridedLoop(o, x, opInto, noconj_o, false, 0, 0, 0);
}

fn stridedLoop(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_o: bool,
    comptime is_c_order: bool,
    o_offset: isize,
    x_offset: usize,
    level: usize,
) void {
    const X = @TypeOf(x);

    if (level == o.ndim) {
        const o_ptr = arrutils.getPtr(o.data, o_offset);

        opInto(
            &o_ptr[0],
            if (comptime noconj_o)
                x.data[x_offset]
            else
                numeric.conj(x.data[x_offset]),
        );

        return;
    }

    const dim = if (comptime is_c_order) level else o.ndim - 1 - level;

    var curr_o = o_offset;
    var curr_x = x_offset;
    for (0..o.shape[dim]) |_| {
        stridedLoop(o, x, opInto, noconj_o, is_c_order, curr_o, curr_x, level + 1);

        curr_o += o.strides[dim];
        curr_x += X.strides[dim];
    }
}
