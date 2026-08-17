const meta = @import("../../meta.zig");

const linalg = @import("../../linalg.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    linalg.blas.gemv(
        meta.layoutOf(@TypeOf(y)).?,
        if (y.flags.noconj) .trans else .conj_trans,
        x.len,
        o.len,
        1,
        y.data,
        y.ld,
        x.data,
        x.inc,
        0,
        o.data,
        o.inc,
    ) catch unreachable;
}
