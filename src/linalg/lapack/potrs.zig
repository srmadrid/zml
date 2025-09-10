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
const Uplo = types.Uplo;

pub inline fn potrs(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Coerce(A, B);

    if (n < 0 or nrhs < 0 or lda < int.max(1, n) or ldb < int.max(1, n))
        return lapack.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return;

    if (comptime !types.isArbitraryPrecision(C)) {
        if (uplo == .upper) {
            // Solve A * X = B, where A = U^T * U or A = U^H * U.

            // Solve U^T * X = B or U^H * X = B, overwriting B with X.
            blas.trsm(
                order,
                .left,
                .upper,
                .conj_trans,
                .non_unit,
                n,
                nrhs,
                1,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;

            // Solve U * X = B, overwriting B with X.
            blas.trsm(
                order,
                .left,
                .upper,
                .no_trans,
                .non_unit,
                n,
                nrhs,
                1,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;
        } else {
            // Solve A * X = B, where A = L * L^T or A = L * L^H.

            // Solve L * X = B, overwriting B with X.
            blas.trsm(
                order,
                .left,
                .lower,
                .no_trans,
                .non_unit,
                n,
                nrhs,
                1,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;

            // Solve L^T * X = B or L^H * X = B, overwriting B with X.
            blas.trsm(
                order,
                .left,
                .lower,
                .conj_trans,
                .non_unit,
                n,
                nrhs,
                1,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.potrs not implemented for arbitrary precision types yet");
    }

    return;
}
