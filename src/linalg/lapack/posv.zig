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

pub fn posv(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: anytype,
    lda: i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !i32 {
    const A: type = types.Child(@TypeOf(a));
    const B: type = types.Child(@TypeOf(b));
    const C: type = types.Coerce(A, B);

    if (n < 0 or nrhs < 0 or lda < int.max(1, n) or ldb < int.max(1, if (order == .col_major) n else nrhs))
        return lapack.Error.InvalidArgument;

    var info: i32 = 0;

    // Quick return if possible.
    if (n == 0 or nrhs == 0)
        return info;

    if (comptime !types.isArbitraryPrecision(C)) {
        // Compute the Cholesky factorization A = U^T * U, A = U^H * U, A = L * L^T, or A = L * L^H.
        info = lapack.potrf(
            order,
            uplo,
            n,
            a,
            lda,
            ctx,
        ) catch unreachable;

        if (info == 0) {
            // Solve the system A * X = B, overwriting B with X.
            lapack.potrs(
                order,
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
