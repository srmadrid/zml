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
const Direction = lapack.Direction;
const Storage = lapack.Storage;

const utils = @import("../utils.zig");

pub fn org2r(
    order: Order,
    m: i32,
    n: i32,
    k: i32,
    a: anytype,
    lda: i32,
    tau: anytype,
    work: anytype,
    ctx: anytype,
) !void {
    if (m < 0 or n < 0 or n > m or k < 0 or k > n or lda < int.max(1, if (order == .col_major) m else n))
        return lapack.Error.InvalidArgument;

    if (n == 0)
        return;

    // Initialise columns k+1:n to columns of the unit matrix
    var j: i32 = k;
    while (j < n) : (j += 1) {
        var l: i32 = 0;
        while (l < m) : (l += 1) {
            try ops.set(
                &a[utils.index(order, l, j, lda)],
                0,
                ctx,
            );
        }

        try ops.set(
            &a[utils.index(order, j, j, lda)],
            1,
            ctx,
        );
    }

    var i: i32 = k - 1;
    while (i >= 0) : (i -= 1) {
        // Apply h(i) to a(i:m,i:n) from the left
        if (i < n - 1) {
            try lapack.larf1f(
                order,
                .left,
                m - i,
                n - i - 1,
                a + utils.index(order, i, i, lda),
                utils.col_ld(order, lda),
                tau[types.scast(u32, i)],
                a + utils.index(order, i, i + 1, lda),
                lda,
                work,
                ctx,
            );
        }

        if (i < m - 1) {
            try blas.scal(
                m - i - 1,
                try ops.neg(tau[types.scast(u32, i)], ctx),
                a + utils.index(order, i + 1, i, lda),
                utils.col_ld(order, lda),
                ctx,
            );
        }

        try ops.sub_(
            &a[utils.index(order, i, i, lda)],
            1,
            tau[types.scast(u32, i)],
            ctx,
        );

        // set a(1:i-1,i) to zero
        var l: i32 = 0;
        while (l < i) : (l += 1) {
            try ops.set(
                &a[utils.index(order, l, i, lda)],
                0,
                ctx,
            );
        }
    }

    return;
}
