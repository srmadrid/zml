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

/// Solves a system of linear equations with a Cholesky-factored symmetric or
/// Hermitian positive-definite coefficient matrix.
///
/// The `potrs` routine solves for `X` the system of linear equations
/// `A * X = B` with a symmetric positive-definite or, for complex data,
/// Hermitian positive-definite matrix `A`, given the Cholesky factorization of
/// `A`:
///
/// ```zig
///     A = U^T * U,
/// ```
///
/// or
///
/// ```zig
///     A = L * L^T,
/// ```
///
/// or
///
/// ```zig
///     A = U^H * U,
/// ```
///
/// or
///
/// ```zig
///     A = L * L^H,
/// ```
///
/// where `U` is an upper triangular matrix and `L` is lower triangular. The
/// system is solved with multiple right-hand sides stored in the columns of the
/// matrix `B`. Before calling this routine, you must call `potrf` to compute
/// the Cholesky factorization of `A`.
///
/// Signature
/// ---------
/// ```zig
/// fn potrs(order: Order, uplo: Uplo, n: i32, nrhs: i32, a: [*]const A, lda: i32, b: [*]B, ldb: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies how the matrix `A` has been factored:
/// - If `uplo = upper`, then the factorization is `A = U^T * U` or
/// `A = U^H * U`.
/// - If `uplo = lower`, then the factorization is `A = L * L^T` or
/// `A = L * L^H`.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (many-item pointer to `bool`, `int`, `float`, `cfloat`, `integer`,
/// `rational`, `real`, `complex` or `expression`): Array, size at least
/// `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `ldb * nrhs` if `order = col_major`, or `ldb * n` if
/// `order = row_major`. On return, contains the solution matrix `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Errors
/// ------
/// `linalg.lapack.Error.InvalidArgument`: If `n` or `nrhs` is less than 0, if
/// `lda` is less than `max(1, n)`, or if `ldb` is less than `max(1, n)` or
/// `max(1, nrhs)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding LAPACKE function, if available. In that case, no errors will
/// be raised even if the arguments are invalid.
pub fn potrs(
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
