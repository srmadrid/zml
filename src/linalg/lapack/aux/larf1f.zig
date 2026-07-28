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

pub fn larf1f(
    order: Order,
    side: Side,
    m: i32,
    n: i32,
    v: anytype,
    incv: i32,
    tau: anytype,
    c: anytype,
    ldc: i32,
    work: anytype,
    ctx: anytype,
) !void {
    const V = types.Child(@TypeOf(v));
    const Ta = @TypeOf(tau);
    const C = types.Child(@TypeOf(c));
    const W = types.Child(@TypeOf(work));
    const CC = types.Coerce(V, types.Coerce(Ta, types.Coerce(C, W)));

    var lastv: i32 = 1;
    var lastc: i32 = 0;
    if (try ops.ne(tau, 0, ctx)) {
        // Set up variables for scanning v. lastv begins pointing to the end
        // of v up to v[0]
        if (side == .left) {
            lastv = m;
        } else {
            lastv = n;
        }

        // Look for the last non-zero row in v
        var i: i32 = if (incv > 0) (lastv - 1) * incv else 0;
        while (lastv > 1 and try ops.eq(v[types.scast(u32, i)], 0, ctx)) {
            lastv -= 1;
            i -= incv;
        }

        if (side == .left) {
            // Scan for the last non-zero column in c[0:lastv,:]
            lastc = try lapack.ilalc(order, lastv, n, c, ldc);
        } else {
            // Scan for the last non-zero row in c[:,0:lastv]
            lastc = try lapack.ilalr(order, m, lastv, c, ldc);
        }
    }

    if (lastc == 0)
        return;

    if (side == .left) {
        // Form h * c
        if (lastv == 1) {
            // c[0,0:lastc-1] := (1 - tau) * c[0,0:lastc-1]
            try blas.scal(
                lastc,
                try ops.sub(1, tau, ctx),
                c,
                utils.row_ld(order, ldc),
                ctx,
            );
        } else {
            // w[0:lastc-1,0] = c[1:lastv-1,0:lastc-1]**h * v[1:lastv-1,0]
            try blas.gemv(
                order,
                if (comptime !types.isComplex(CC)) .trans else .conj_trans,
                lastv - 1,
                lastc,
                1,
                c + utils.index(order, 1, 0, ldc),
                ldc,
                v + types.scast(u32, incv),
                incv,
                0,
                work,
                1,
                ctx,
            );

            // w[0:lastc-1,0] += v[0,0] * c[0,0:lastc-1]**h
            if (comptime !types.isComplex(CC)) {
                try blas.axpy(
                    lastc,
                    1,
                    c,
                    utils.row_ld(order, ldc),
                    work,
                    1,
                    ctx,
                );
            } else {
                var i: i32 = 0;
                while (i < lastc) : (i += 1) {
                    try ops.add_(
                        &work[types.scast(u32, i)],
                        work[types.scast(u32, i)],
                        try ops.conj(c[utils.index(order, 0, i, ldc)], ctx),
                        ctx,
                    );
                }
            }

            // c[0, 0:lastc-1] -= tau * w[0:lastc-1,0]**h
            if (comptime !types.isComplex(CC)) {
                try blas.axpy(
                    lastc,
                    try ops.neg(tau, ctx),
                    work,
                    1,
                    c,
                    utils.row_ld(order, ldc),
                    ctx,
                );
            } else {
                var i: i32 = 0;
                while (i < lastc) : (i += 1) {
                    try ops.sub_(
                        &c[utils.index(order, 0, i, ldc)],
                        c[utils.index(order, 0, i, ldc)],
                        try ops.mul(
                            tau,
                            try ops.conj(work[types.scast(u32, i)], ctx),
                            ctx,
                        ),
                        ctx,
                    );
                }
            }

            // c[1:lastv-1,0:lastc-1] += -tau * v[1:lastv-1,0] * w[0:lastc-1,0]**h
            if (comptime !types.isComplex(CC)) {
                try blas.ger(
                    order,
                    lastv - 1,
                    lastc,
                    try ops.neg(tau, ctx),
                    v + types.scast(u32, incv),
                    incv,
                    work,
                    1,
                    c + utils.index(order, 1, 0, ldc),
                    ldc,
                    ctx,
                );
            } else {
                try blas.gerc(
                    order,
                    lastv - 1,
                    lastc,
                    try ops.neg(tau, ctx),
                    v + types.scast(u32, incv),
                    incv,
                    work,
                    1,
                    c + utils.index(order, 1, 0, ldc),
                    ldc,
                    ctx,
                );
            }
        }
    } else {
        // Form c * h
        if (lastv == 1) {
            // c[0:lastc-1,0] := (1 - tau) * c[0:lastc-1,0]
            try blas.scal(
                lastc,
                try ops.sub(1, tau, ctx),
                c,
                utils.col_ld(order, ldc),
                ctx,
            );
        } else {
            // w[0:lastc-1,0] := c[0:lastc-1,1:lastv-1] * v[1:lastv-1,0]
            try blas.gemv(
                order,
                .no_trans,
                lastc,
                lastv - 1,
                1,
                c + utils.index(order, 0, 1, ldc),
                ldc,
                v + types.scast(u32, incv),
                incv,
                0,
                work,
                1,
                ctx,
            );

            // w[0:lastc-1,0] += v[0,0] * c[0:lastc-1,0]
            try blas.axpy(
                lastc,
                1,
                c,
                utils.col_ld(order, ldc),
                work,
                1,
                ctx,
            );

            // c[0:lastc-1,0] += -tau * v[0,0] * w[0:lastc-1,0]
            try blas.axpy(
                lastc,
                try ops.neg(tau, ctx),
                work,
                1,
                c,
                utils.col_ld(order, ldc),
                ctx,
            );

            // c[0:lastc-1,1:lastv-1] += -tau * w[0:lastc-1,0] * v[1:lastv-1]**h
            if (comptime !types.isComplex(CC)) {
                try blas.ger(
                    order,
                    lastc,
                    lastv - 1,
                    try ops.neg(tau, ctx),
                    work,
                    1,
                    v + types.scast(u32, incv),
                    incv,
                    c + utils.index(order, 0, 1, ldc),
                    ldc,
                    ctx,
                );
            } else {
                try blas.gerc(
                    order,
                    lastc,
                    lastv - 1,
                    try ops.neg(tau, ctx),
                    work,
                    1,
                    v + types.scast(u32, incv),
                    incv,
                    c + utils.index(order, 0, 1, ldc),
                    ldc,
                    ctx,
                );
            }
        }
    }
}
