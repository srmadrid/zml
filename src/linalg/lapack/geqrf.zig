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

pub fn geqrf(
    order: Order,
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    tau: anytype,
    work: anytype,
    lwork: i32,
    ctx: anytype,
) !void {
    if (m < 0 or n < 0 or lda < int.max(1, if (order == .col_major) m else n) or
        (lwork < 1 and (lwork != -1 or (m > 0 and lwork < int.max(1, n)))))
        return lapack.Error.InvalidArgument;

    const k: i32 = int.min(m, n);
    var nb = lapack.ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
    if (lwork == -1) {
        const lwkopt: i32 = if (k == 0) 1 else n * nb;

        try ops.set(
            &work[0],
            lwkopt,
            ctx,
        );

        return;
    }

    if (k == 0) {
        try ops.set(
            &work[0],
            1,
            ctx,
        );

        return;
    }

    var nbmin: i32 = 2;
    var nx: i32 = 0;
    var iws: i32 = n;
    var ldwork: i32 = 0;
    if (nb > 1 and nb < k) {
        // Determine when to cross over from blocked to unblocked code.
        nx = int.max(0, lapack.ilaenv(3, "DGEQRF", " ", m, n, -1, -1));
        if (nx < k) {
            // Determine if workspace is large enough for blocked code.
            ldwork = if (order == .col_major) n else nb;
            iws = ldwork * nb;
            if (lwork < iws) {
                // Not enough workspace to use optimal nb:  reduce nb and
                // determine the minimum value of nb.
                nb = int.div(lwork, ldwork);
                nbmin = int.max(2, lapack.ilaenv(2, "DGEQRF", " ", m, n, -1, -1));
            }
        }
    }

    var i: i32 = 0;
    if (nb >= nbmin and nb < k and nx < k) {
        // use blocked code initially
        while (i < k - nx) : (i += nb) {
            const ib: i32 = int.min(k - i, nb);

            // Compute the qr factorization of the current block
            // a(i:m,i:i+ib-1)
            std.debug.print("Calling geqr2 with i = {d}, ib = {d}\n", .{ i, ib });
            try lapack.geqr2(
                order,
                m - i,
                ib,
                a + utils.index(order, i, i, lda),
                lda,
                tau + types.scast(u32, i),
                work,
                ctx,
            );

            if (i + ib < n) {
                // Form the triangular factor of the block reflector
                // h = h(i) h(i+1) . . . h(i+ib-1)
                std.debug.print("Calling larft with i = {d}, ib = {d}\n", .{ i, ib });
                try lapack.larft(
                    order,
                    .forward,
                    .columnwise,
                    m - i,
                    ib,
                    a + utils.index(order, i, i, lda),
                    lda,
                    tau + types.scast(u32, i),
                    work,
                    ldwork,
                    ctx,
                );

                // Apply h**h to a(i:m,i+ib:n) from the left
                std.debug.print("Calling larfb with i = {d}, ib = {d}\n", .{ i, ib });
                try lapack.larfb(
                    order,
                    .left,
                    .conj_trans,
                    .forward,
                    .columnwise,
                    m - i,
                    n - i - ib,
                    ib,
                    a + utils.index(order, i, i, lda),
                    lda,
                    work,
                    ldwork,
                    a + utils.index(order, i, i + ib, lda),
                    lda,
                    work + utils.index(order, ib, 0, ldwork),
                    ldwork,
                    ctx,
                );
            }
        }
    }

    // Use unblocked code to factor the last or only block.
    if (i < k) {
        std.debug.print("Calling geqr2 with i = {d}, ib = {d}\n", .{ i, n - i });
        try lapack.geqr2(
            order,
            m - i,
            n - i,
            a + utils.index(order, i, i, lda),
            lda,
            tau + types.scast(u32, i),
            work,
            ctx,
        );
    }

    try ops.set(
        &work[0],
        iws,
        ctx,
    );

    return;
}
