const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;

pub fn getrs(
    order: Order,
    transa: Transpose,
    n: i32,
    nrhs: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Coerce(A, B);

    const nota: bool = transa == .no_trans or transa == .conj_no_trans;

    if (n < 0 or nrhs < 0 or lda < int.max(1, n) or ldb < int.max(1, if (order == .col_major) n else nrhs))
        return lapack.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return;

    if (comptime !types.isArbitraryPrecision(C)) {
        if (nota) {
            // Solve A * X = B or conj(A) * X = B.

            // Apply row interchanges to the right hand sides.
            lapack.laswp(
                order,
                nrhs,
                b,
                ldb,
                1,
                n,
                ipiv,
                1,
            ) catch unreachable;

            // Solve L * X = B or conj(L) * X = B, overwriting B with X.
            blas.trsm(
                order,
                .left,
                .lower,
                transa,
                .unit,
                n,
                nrhs,
                1,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;

            // Solve U * X = B or conj(U) * X = B, overwriting B with X.
            blas.trsm(
                order,
                .left,
                .upper,
                transa,
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
            // Solve A^T * X = B or A^H * X = B.

            // Solve U^T * X = B or U^H * X = B, overwriting B with X.
            blas.trsm(
                order,
                .left,
                .upper,
                transa,
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
                transa,
                .unit,
                n,
                nrhs,
                1,
                a,
                lda,
                b,
                ldb,
                ctx,
            ) catch unreachable;

            // Apply row interchanges to the solution vectors.
            lapack.laswp(
                order,
                nrhs,
                b,
                ldb,
                1,
                n,
                ipiv,
                -1,
            ) catch unreachable;
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.getrs not implemented for arbitrary precision types yet");
    }

    return;
}
