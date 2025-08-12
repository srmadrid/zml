const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;

pub inline fn gesv(
    order: Order,
    n: i32,
    nrhs: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !i32 {
    if (order == .col_major) {
        return k_gesv_c(n, nrhs, a, lda, ipiv, b, ldb, ctx);
    } else {
        return k_gesv_r(n, nrhs, a, lda, ipiv, b, ldb, ctx);
    }
}

fn k_gesv_c(
    n: i32,
    nrhs: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Coerce(A, B);

    if (n < 0 or nrhs < 0 or lda < int.max(1, n) or ldb < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(C)) {
        // Compute the LU factorization of A.
        info = lapack.getrf(
            .col_major,
            n,
            n,
            a,
            lda,
            ipiv,
            ctx,
        ) catch unreachable;

        if (info == 0) {
            // Solve the system A * X = B, overwriting B with X.
            lapack.getrs(
                .col_major,
                .no_trans,
                n,
                nrhs,
                a,
                lda,
                ipiv,
                b,
                ldb,
                ctx,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.gesv not implemented for arbitrary precision types yet");
    }

    return info;
}

fn k_gesv_r(
    n: i32,
    nrhs: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Coerce(A, B);

    if (n < 0 or nrhs < 0 or lda < int.max(1, n) or ldb < int.max(1, nrhs))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(C)) {
        // Compute the LU factorization of A.
        info = lapack.getrf(
            .row_major,
            n,
            n,
            a,
            lda,
            ipiv,
            ctx,
        ) catch unreachable;

        if (info == 0) {
            // Solve the system A * X = B, overwriting B with X.
            lapack.getrs(
                .row_major,
                .no_trans,
                n,
                nrhs,
                a,
                lda,
                ipiv,
                b,
                ldb,
                ctx,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.gesv not implemented for arbitrary precision types yet");
    }

    return info;
}
