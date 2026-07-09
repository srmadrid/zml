const std = @import("std");

const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

const arrutils = @import("../../utils.zig");

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

    if (((comptime O.storage_order == .f) and x.isContiguous(.f)) or ((comptime O.storage_order == .c) and x.isContiguous(.c))) {
        for (0..O.size) |i| {
            opInto(
                &o.data[i],
                if (comptime noconj_x)
                    x.data[i]
                else
                    numeric.conj(x.data[i]),
            );
        }

        return;
    }

    stridedLoop(o, x, opInto, noconj_x, 0, 0, 0);
}

fn stridedLoop(
    o: anytype,
    x: anytype,
    comptime opInto: anytype,
    comptime noconj_x: bool,
    o_offset: usize,
    x_offset: isize,
    comptime level: usize,
) void {
    const O = meta.Child(@TypeOf(o));

    if (comptime level == O.ndim) {
        const x_ptr = arrutils.getPtr(x.data, x_offset);

        opInto(
            &o.data[o_offset],
            if (comptime noconj_x)
                x_ptr[0]
            else
                numeric.conj(x_ptr[0]),
        );

        return;
    }

    const dim = if (comptime O.storage_order == .c) level else O.ndim - 1 - level;

    var curr_o = o_offset;
    var curr_x = x_offset;
    for (0..O.shape[dim]) |_| {
        stridedLoop(o, x, opInto, noconj_x, curr_o, curr_x, level + 1);

        curr_o += O.strides[dim];
        curr_x += x.strides[dim];
    }
}
