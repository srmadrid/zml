const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;

/// Solves a system of linear equations with an LU-factored square coefficient
/// matrix, with multiple right-hand sides.
///
/// The `getrs` routine solves for `X` the following systems of linear
/// equations:
///
/// ```zig
///     A * X = B,
/// ```
///
/// or
///
/// ```zig
///     A^T * X = B,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * X = B,
/// ```
///
/// or
///
/// ```zig
///     A^H * X = B,
/// ```
///
/// where `A` is the LU factorization of a general `n`-by-`n` matrix `A`,
/// computed by `getrf`, `B` is an `n`-by-`nrhs` matrix of right-hand
/// sides, and `X` is an `n`-by-`nrhs` matrix of solutions.
///
/// Signature
/// ---------
/// ```zig
/// fn getrs(order: Order, transa: Transpose, n: i32, nrhs: i32, a: [*]const A, lda: i32, ipiv: [*]const i32, b: [*]B, ldb: i32, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `transa` (`Transpose`): Specifies the form of the system of equations to
/// solve:
/// - If `transa = .no_trans`, then `A * X = B`.
/// - If `transa = .trans`, then `A^T * X = B`.
/// - If `transa = .conj_no_trans`, then `conj(A) * X = B`.
/// - If `transa = .conj_trans`, then `A^H * X = B`.
///
/// `n` (`i32`): The order of the matrix `A` and the number of rows of the
/// matrix `B`. Must be greater than or equal to 0.
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
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. The pivot indices as
/// returned by `getrf`.
///
/// `b` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `ldb * nrhs` if `order = .col_major` or `ldb * n` if
/// `order = .row_major`. On return, contains the solution matrix `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = .col_major` or `max(1, nrhs)` if
/// `order = .row_major`.
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
