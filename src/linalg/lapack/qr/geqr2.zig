const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Side = linalg.Side;

const utils = @import("../utils.zig");

pub fn geqr2(
    order: Order,
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    tau: anytype,
    work: anytype,
    ctx: anytype,
) !void {
    if (m < 0 or n < 0 or lda < int.max(1, if (order == .col_major) m else n))
        return lapack.Error.InvalidArgument;

    var i: i32 = 0;
    while (i < int.min(m, n)) : (i += 1) {
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        try lapack.larfg(
            m - i,
            &a[utils.index(order, i, i, lda)],
            a + utils.index(order, int.min(i + 1, m - 1), i, lda),
            utils.col_ld(order, lda),
            &tau[types.scast(u32, i)],
            ctx,
        );

        if (i < n - 1) {
            // Apply H(i) to A(i:m,i+1:n) from the left
            try lapack.larf1f(
                order,
                .left,
                m - i,
                n - i - 1,
                a + utils.index(order, i, i, lda),
                utils.col_ld(order, lda),
                try ops.conj(tau[types.scast(u32, i)], ctx),
                a + utils.index(order, i, i + 1, lda),
                lda,
                work,
                ctx,
            );
        }
    }
}
