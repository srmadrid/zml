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

/// Computes the solution to the system of linear equations with a symmetric or
/// Hermitian positive-definite coefficient matrix `A` and multiple right-hand
/// sides.
///
/// The `posv` routine solves for `X` the real or complex system of linear
/// equations `A * X = B`, where `A` is an `n`-by-`n` symmetric or Hermitian
/// positive-definite matrix, the columns of matrix `B` are individual
/// right-hand sides, and the columns of `X` are the corresponding solutions.
/// The Cholesky decomposition is used to factor `A` as:
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
/// factored form of `A` is then used to solve the system of equations
/// `A * X = B`.
///
/// Signature
/// ---------
/// ```zig
/// fn posv(order: Order, uplo: Uplo, n: i32, nrhs: i32, a: [*]A, lda: i32, b: [*]B, ldb: i32, ctx: anytype) !i32
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies which part of the matrix `A` is stored, and
/// which factorization is computed:
/// - If `uplo = upper`, then the upper triangular part of `A` is stored, and
/// the factorization is `A = U^T * U` or `A = U^H * U` is computed.
/// - If `uplo = lower`, then the lower triangular part of `A` is stored, and
/// the factorization is `A = L * L^T` or `A = L * L^H` is computed.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `lda * n`. On return, contains the Cholesky factorization of the
/// matrix `A`.
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
/// `void`: The result is stored in `a` and `b`.
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
