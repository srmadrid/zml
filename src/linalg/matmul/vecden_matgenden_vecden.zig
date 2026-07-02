const meta = @import("../../meta.zig");

const linalg = @import("../../linalg.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    linalg.blas.gemv(
        meta.layoutOf(@TypeOf(x)).?,
        if (x.flags.noconj) .no_trans else .conj_no_trans,
        o.len,
        y.len,
        1,
        x.data,
        x.ld,
        y.data,
        y.inc,
        0,
        o.data,
        o.inc,
        .{ .num_threads = 1 },
    ) catch unreachable;
}
