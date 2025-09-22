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
const Uplo = types.Uplo;

const utils = @import("../utils.zig");

pub fn larfg(
    n: i32,
    alpha: anytype,
    x: anytype,
    incx: i32,
    tau: anytype,
    ctx: anytype,
) !void {
    const X: type = types.Child(@TypeOf(x));
    const Al: type = types.Child(@TypeOf(alpha));
    const C: type = types.Coerce(X, Al);

    if (n <= 0) {
        try ops.set(
            tau,
            0,
            ctx,
        );

        return;
    }

    var xnorm: types.Scalar(X) = try blas.nrm2(
        n - 1,
        x,
        incx,
        ctx,
    );

    if (try ops.eq(xnorm, 0, ctx) and try ops.eq(try ops.im(alpha.*, ctx), 0, ctx)) {
        // h  =  i
        try ops.set(
            tau,
            0,
            ctx,
        );
    } else {
        // General case
        var beta: types.Scalar(C) = try ops.neg(try ops.copysign(
            lapack.lapy3(
                try ops.re(alpha.*, ctx),
                try ops.im(alpha.*, ctx),
                xnorm,
                ctx,
            ),
            try ops.re(alpha.*, ctx),
            ctx,
        ), ctx);
        const safmin: types.Scalar(X) = try ops.div(
            lapack.lamch(types.Scalar(X), .sfmin),
            lapack.lamch(types.Scalar(X), .eps),
            ctx,
        );
        const rsafmn: types.Scalar(X) = try ops.div(
            1,
            safmin,
            ctx,
        );

        var knt: i32 = 0;
        if (try ops.lt(try ops.abs(beta, ctx), safmin, ctx)) {
            // xnorm, beta may be inaccurate; scale x and recompute them
            while ((try ops.lt(try ops.abs(beta, ctx), safmin, ctx)) and (knt < 20)) : (knt += 1) {
                try blas.scal(
                    n - 1,
                    rsafmn,
                    x,
                    incx,
                    ctx,
                );

                try ops.mul_(
                    &beta,
                    beta,
                    rsafmn,
                    ctx,
                );
                try ops.mul_(
                    alpha,
                    alpha.*,
                    rsafmn,
                    ctx,
                );
            }

            // new beta is at most 1, at least safmin
            xnorm = try blas.nrm2(
                n - 1,
                x,
                incx,
                ctx,
            );

            beta = try ops.neg(try ops.copysign(
                lapack.lapy3(
                    try ops.re(alpha.*, ctx),
                    try ops.im(alpha.*, ctx),
                    xnorm,
                    ctx,
                ),
                try ops.re(alpha.*, ctx),
                ctx,
            ), ctx);
        }

        try ops.set(
            tau,
            try ops.div(
                try ops.sub(beta, alpha.*, ctx),
                beta,
                ctx,
            ),
            ctx,
        );

        try ops.div_(
            alpha,
            1,
            try ops.sub(
                alpha.*,
                beta,
                ctx,
            ),
            ctx,
        );

        try blas.scal(
            n - 1,
            alpha.*,
            x,
            incx,
            ctx,
        );

        // if alpha is subnormal, it may lose relative accuracy
        var j: i32 = 0;
        while (j < knt) : (j += 1) {
            try ops.mul_(
                &beta,
                beta,
                safmin,
                ctx,
            );
        }

        try ops.set(
            alpha,
            beta,
            ctx,
        );
    }

    return;
}
