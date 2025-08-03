const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = linalg.Order;
const Uplo = linalg.Uplo;

pub inline fn posv(
    order: Order,
    uplo: Uplo,
    n: isize,
    nrhs: isize,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    ctx: anytype,
) !isize {
    if (order == .col_major) {
        return k_posv_c(uplo, n, nrhs, a, lda, b, ldb, ctx);
    } else {
        return k_posv_r(uplo, n, nrhs, a, lda, b, ldb, ctx);
    }
}

fn k_posv_c(
    uplo: Uplo,
    n: isize,
    nrhs: isize,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    ctx: anytype,
) !isize {
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Coerce(A, B);

    if (n < 0 or nrhs < 0 or lda < int.max(1, n) or ldb < int.max(1, n))
        return lapack.Error.InvalidArgument;

    var info: isize = 0;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(C)) {
        // Compute the Cholesky factorization A = U^T * U, A = U^H * U, A = L * L^T, or A = L * L^H.
        info = lapack.potrf(
            .col_major,
            uplo,
            n,
            a,
            lda,
            ctx,
        ) catch unreachable;

        if (info == 0) {
            // Solve the system A * X = B, overwriting B with X.
            lapack.potrs(
                .col_major,
                uplo,
                n,
                nrhs,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.posv not implemented for arbitrary precision types yet");
    }

    return info;
}

fn k_posv_r(
    uplo: Uplo,
    n: isize,
    nrhs: isize,
    a: anytype,
    lda: isize,
    b: anytype,
    ldb: isize,
    ctx: anytype,
) !isize {
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Coerce(A, B);

    if (n < 0 or nrhs < 0 or lda < int.max(1, n) or ldb < int.max(1, nrhs))
        return lapack.Error.InvalidArgument;

    var info: isize = 0;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(C)) {
        // Compute the Cholesky factorization A = U^T * U, A = U^H * U, A = L * L^T, or A = L * L^H.
        info = lapack.potrf(
            .row_major,
            uplo,
            n,
            a,
            lda,
            ctx,
        ) catch unreachable;

        if (info == 0) {
            // Solve the system A * X = B, overwriting B with X.
            lapack.potrs(
                .row_major,
                uplo,
                n,
                nrhs,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.posv not implemented for arbitrary precision types yet");
    }

    return info;
}
