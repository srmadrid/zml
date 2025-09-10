const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;

pub fn gesv(
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
        // Compute the LU factorization of A.
        info = lapack.getrf(
            order,
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
                order,
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
