const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Transpose = linalg.Transpose;

/// Computes the solution to the system of linear equations with a square
/// coefficient matrix `A` and multiple right-hand sides.
///
/// The `gesv` routine solves for `X` the system of linear equations
///
/// ```zig
///     A * X = B,
/// ```
///
/// where `A` is an `n`-by-`n` matrix, the columns of matrix `B` are individual
/// right-hand sides, and the columns of `X` are the corresponding solutions.
///
/// The LU decomposition with partial pivoting and row interchanges is used to
/// factor `A` as `A = P * L * U`, where `P` is a permutation matrix, `L` is
/// unit lower triangular, and `U` is upper triangular. The factored form of `A`
/// is then used to solve the system of equations `A * X = B`.
///
/// Signature
/// ---------
/// ```zig
/// fn gesv(order: Order, n: i32, nrhs: i32, a: [*]A, lda: i32, ipiv: [*]i32, b: B[*], ldb: i32, ctx: anytype) !i32
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `n` (`i32`): The order of the matrix `A` and the number of rows of the
/// matrix `B`. Must be greater than or equal to 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `lda * n`. On return, contains the LU factorization of the matrix `A`
/// as computed by `getrf`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. On return contains the
/// pivot indices as returned by `getrf`. For `1 <= i <= n`, row `i` of the
/// matrix was interchanged with row `ipiv[i - 1]`.
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
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`, `ipiv`, and `b`.
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
