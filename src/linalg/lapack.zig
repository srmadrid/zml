const std = @import("std");
const options = @import("options");

const types = @import("../types.zig");
const scast = types.scast;
const validateContext = types.validateContext;
const Scalar = types.Scalar;
const Child = types.Child;
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const cfloat = @import("../cfloat.zig");
const cf32 = cfloat.cf32;
const cf64 = cfloat.cf64;
const ops = @import("../ops.zig");
const linalg = @import("../linalg.zig");

const ci = @import("../c.zig");

const Order = types.Order;
const Transpose = linalg.Transpose;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Side = linalg.Side;

pub const Mach = enum {
    eps,
    sfmin,
    base,
    prec,
    t,
    rnd,
    emin,
    rmin,
    emax,
    rmax,

    pub inline fn toChar(self: Mach) u8 {
        return switch (self) {
            .eps => 'E',
            .sfmin => 'S',
            .base => 'B',
            .prec => 'P',
            .t => 'N',
            .rnd => 'R',
            .emin => 'M',
            .rmin => 'U',
            .emax => 'L',
            .rmax => 'O',
        };
    }
};

//
pub const ilaenv = @import("lapack/ilaenv.zig").ilaenv;

pub inline fn lamch(
    comptime T: type,
    cmach: Mach,
) Scalar(T) {
    comptime if (!types.isNumeric(T))
        @compileError("zml.linalg.lapack.lamch requires T to be a numeric, got " ++ @typeName(T));

    if (comptime options.link_lapacke != null) {
        switch (comptime types.numericType(T)) {
            .float => {
                if (comptime T == f32) {
                    return ci.LAPACKE_slamch(cmach.toChar());
                } else if (comptime T == f64) {
                    return ci.LAPACKE_dlamch(cmach.toChar());
                }
            },
            .cfloat => {
                if (comptime Scalar(T) == f32) {
                    return ci.LAPACKE_clamch(cmach.toChar());
                } else if (comptime Scalar(T) == f64) {
                    return ci.LAPACKE_zlamch(cmach.toChar());
                }
            },
            else => {},
        }
    }

    return @import("lapack/lamch.zig").lamch(comptime T, cmach);
}

/// Performs a series of row interchanges on a general rectangular matrix.
///
/// The routine performs a series of row interchanges on the `m`-by-`n` matrix
/// `A`. One row interchange is initiated for each of rows `k1` through `k2` of
/// `A`.
///
/// Signature
/// ---------
/// ```zig
/// fn laswp(order: Order, n: i32, a: [*]A, lda: i32, k1: i32, k2: i32, ipiv: [*]const i32, incx: i32) !void
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size
/// at least `lda * n` if `order = .col_major` or `lda * mm` if
/// `order = .row_major`, where `mm` is not less than the maximum of
/// `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`i32`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`i32`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`i32`): The increment for the elements of `ipiv`. If `incx` is
/// negative, the pivots are applied in reverse order. Must be different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Errors
/// ------
/// `linalg.lapack.Error.InvalidArgument`: If `n` is less than or equal to 0, if
/// `lda` is less than `max(1, n)` (only checked for `order = .row_major`), or
/// if `incx` is 0.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding LAPACKE function, if available. In that case, no errors will
/// be raised even if the arguments are invalid.
pub inline fn laswp(
    order: Order,
    n: i32,
    a: anytype,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) !void {
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.laswp requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.laswp requires a's child type to be a numeric, got " ++ @typeName(A));

    if (comptime options.link_lapacke != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    const result = ci.LAPACKE_slaswp(order.toCInt(), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

                    if (result < 0) {
                        if (result == -1011)
                            return std.mem.Allocator.Error.OutOfMemory;

                        return Error.InvalidArgument;
                    }

                    return;
                } else if (comptime A == f64) {
                    const result = ci.LAPACKE_dlaswp(order.toCInt(), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

                    if (result < 0) {
                        if (result == -1011)
                            return std.mem.Allocator.Error.OutOfMemory;

                        return Error.InvalidArgument;
                    }

                    return;
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    const result = ci.LAPACKE_claswp(order.toCInt(), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

                    if (result < 0) {
                        if (result == -1011)
                            return std.mem.Allocator.Error.OutOfMemory;

                        return Error.InvalidArgument;
                    }

                    return;
                } else if (comptime Scalar(A) == f64) {
                    const result = ci.LAPACKE_zlaswp(order.toCInt(), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

                    if (result < 0) {
                        if (result == -1011)
                            return std.mem.Allocator.Error.OutOfMemory;

                        return Error.InvalidArgument;
                    }

                    return;
                }
            },
            else => {},
        }
    }

    return @import("lapack/laswp.zig").laswp(order, n, a, lda, k1, k2, ipiv, incx);
}

/// Performs a series of row interchanges on a general rectangular matrix.
///
/// The routine performs a series of row interchanges on the `m`-by-`n` matrix
/// `A`. One row interchange is initiated for each of rows `k1` through `k2` of
/// `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`i32`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`i32`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`i32`): The increment for the elements of `ipiv`. If `incx` is
/// negative, the pivots are applied in reverse order. Must be different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn slaswp(
    order: Order,
    n: i32,
    a: [*]f32,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) void {
    return laswp(order, n, a, lda, k1, k2, ipiv, incx) catch {};
}

/// Performs a series of row interchanges on a general rectangular matrix.
///
/// The routine performs a series of row interchanges on the `m`-by-`n` matrix
/// `A`. One row interchange is initiated for each of rows `k1` through `k2` of
/// `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`i32`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`i32`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`i32`): The increment for the elements of `ipiv`. If `incx` is
/// negative, the pivots are applied in reverse order. Must be different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dlaswp(
    order: Order,
    n: i32,
    a: [*]f64,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) void {
    return laswp(order, n, a, lda, k1, k2, ipiv, incx) catch {};
}

/// Performs a series of row interchanges on a general rectangular matrix.
///
/// The routine performs a series of row interchanges on the `m`-by-`n` matrix
/// `A`. One row interchange is initiated for each of rows `k1` through `k2` of
/// `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`i32`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`i32`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`i32`): The increment for the elements of `ipiv`. If `incx` is
/// negative, the pivots are applied in reverse order. Must be different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn claswp(
    order: Order,
    n: i32,
    a: [*]cf32,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) void {
    return laswp(order, n, a, lda, k1, k2, ipiv, incx) catch {};
}

/// Performs a series of row interchanges on a general rectangular matrix.
///
/// The routine performs a series of row interchanges on the `m`-by-`n` matrix
/// `A`. One row interchange is initiated for each of rows `k1` through `k2` of
/// `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`i32`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`i32`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`i32`): The increment for the elements of `ipiv`. If `incx` is
/// negative, the pivots are applied in reverse order. Must be different from 0.
///
/// Returns
/// -------
/// `void`: The result is stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zlaswp(
    order: Order,
    n: i32,
    a: [*]cf64,
    lda: i32,
    k1: i32,
    k2: i32,
    ipiv: [*]const i32,
    incx: i32,
) void {
    return laswp(order, n, a, lda, k1, k2, ipiv, incx) catch {};
}

/// Computes LU factorization using partial pivoting with row interchanges.
///
/// The `getrf2` routine computes an LU factorization of a general `m`-by-`n`
/// matrix `A` using partial pivoting with row interchanges. The factorization
/// has the form:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// This is the recursive version of the algorithm. It divides the matrix into
/// four submatrices:
///
/// ```zig
///         [ A11 A12 ]
///     A = [ A21 A22 ],
/// ```
///
/// where `A11` is `n1` by `n1` and `A22` is `n2` by `n2` with `n1 = min(m, n)`,
/// and `n2 = n - n1`. The subroutine calls itself to factor `[ A11 A12 ]`, does
/// the swaps on `[ A12 A22 ]^T`, solves `A12`, updates `A22`, and then it calls
/// itself to factor `A22` and do the swaps on `A21`.
///
/// Signature
/// ---------
/// ```zig
/// fn getrf2(order: Order, m: i32, n: i32, a: [*]A, lda: i32, ipiv: [*]i32, ctx: anytype) !i32
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size
/// at least `lda * n` if `order = .col_major` or `lda * m` if
/// `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Errors
/// ------
/// `linalg.lapack.Error.InvalidArgument`: If `m` or `n` is less than 0, or if
/// `lda` is less than `max(1, m)` or `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding LAPACKE function, if available. In that case, no errors will
/// be raised even if the arguments are invalid.
pub inline fn getrf2(
    order: Order,
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    ctx: anytype,
) !i32 {
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.getrf2 requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.getrf2 requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(A)) {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.getrf2 not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime options.link_lapacke != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return scast(i32, ci.LAPACKE_sgetrf2(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime A == f64) {
                    return scast(i32, ci.LAPACKE_dgetrf2(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return scast(i32, ci.LAPACKE_cgetrf2(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime Scalar(A) == f64) {
                    return scast(i32, ci.LAPACKE_zgetrf2(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                }
            },
            else => {},
        }
    }

    return @import("lapack/getrf2.zig").getrf2(order, m, n, a, lda, ipiv, ctx);
}

/// Computes LU factorization using partial pivoting with row interchanges.
///
/// The `sgetrf2` routine computes an LU factorization of a general `m`-by-`n`
/// matrix `A` using partial pivoting with row interchanges. The factorization
/// has the form:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// This is the recursive version of the algorithm. It divides the matrix into
/// four submatrices:
///
/// ```zig
///         [ A11 A12 ]
///     A = [ A21 A22 ],
/// ```
///
/// where `A11` is `n1` by `n1` and `A22` is `n2` by `n2` with `n1 = min(m, n)`,
/// and `n2 = n - n1`. The subroutine calls itself to factor `[ A11 A12 ]`, does
/// the swaps on `[ A12 A22 ]^T`, solves `A12`, updates `A22`, and then it calls
/// itself to factor `A22` and do the swaps on `A21`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn sgetrf2(
    order: Order,
    m: i32,
    n: i32,
    a: [*]f32,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf2(order, m, n, a, lda, ipiv, null) catch -1;
}

/// Computes LU factorization using partial pivoting with row interchanges.
///
/// The `dgetrf2` routine computes an LU factorization of a general `m`-by-`n`
/// matrix `A` using partial pivoting with row interchanges. The factorization
/// has the form:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// This is the recursive version of the algorithm. It divides the matrix into
/// four submatrices:
///
/// ```zig
///         [ A11 A12 ]
///     A = [ A21 A22 ],
/// ```
///
/// where `A11` is `n1` by `n1` and `A22` is `n2` by `n2` with `n1 = min(m, n)`,
/// and `n2 = n - n1`. The subroutine calls itself to factor `[ A11 A12 ]`, does
/// the swaps on `[ A12 A22 ]^T`, solves `A12`, updates `A22`, and then it calls
/// itself to factor `A22` and do the swaps on `A21`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dgetrf2(
    order: Order,
    m: i32,
    n: i32,
    a: [*]f64,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf2(order, m, n, a, lda, ipiv, null) catch -1;
}

/// Computes LU factorization using partial pivoting with row interchanges.
///
/// The `cgetrf2` routine computes an LU factorization of a general `m`-by-`n`
/// matrix `A` using partial pivoting with row interchanges. The factorization
/// has the form:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// This is the recursive version of the algorithm. It divides the matrix into
/// four submatrices:
///
/// ```zig
///         [ A11 A12 ]
///     A = [ A21 A22 ],
/// ```
///
/// where `A11` is `n1` by `n1` and `A22` is `n2` by `n2` with `n1 = min(m, n)`,
/// and `n2 = n - n1`. The subroutine calls itself to factor `[ A11 A12 ]`, does
/// the swaps on `[ A12 A22 ]^T`, solves `A12`, updates `A22`, and then it calls
/// itself to factor `A22` and do the swaps on `A21`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cgetrf2(
    order: Order,
    m: i32,
    n: i32,
    a: [*]cf32,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf2(order, m, n, a, lda, ipiv, null) catch -1;
}

/// Computes LU factorization using partial pivoting with row interchanges.
///
/// The `zgetrf2` routine computes an LU factorization of a general `m`-by-`n`
/// matrix `A` using partial pivoting with row interchanges. The factorization
/// has the form:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// This is the recursive version of the algorithm. It divides the matrix into
/// four submatrices:
///
/// ```zig
///         [ A11 A12 ]
///     A = [ A21 A22 ],
/// ```
///
/// where `A11` is `n1` by `n1` and `A22` is `n2` by `n2` with `n1 = min(m, n)`,
/// and `n2 = n - n1`. The subroutine calls itself to factor `[ A11 A12 ]`, does
/// the swaps on `[ A12 A22 ]^T`, solves `A12`, updates `A22`, and then it calls
/// itself to factor `A22` and do the swaps on `A21`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zgetrf2(
    order: Order,
    m: i32,
    n: i32,
    a: [*]cf64,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf2(order, m, n, a, lda, ipiv, null) catch -1;
}

/// Computes the LU factorization of a general `m`-by-`n` matrix.
///
/// The `getrf` routine computes the LU factorization of a general `m`-by-`n`
/// matrix `A` as:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// Signature
/// ---------
/// ```zig
/// fn getrf(order: Order, m: i32, n: i32, a: [*]A, lda: i32, ipiv: [*]i32, ctx: anytype) !i32
/// ```
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size
/// at least `lda * n` if `order = .col_major` or `lda * m` if
/// `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Errors
/// ------
/// `linalg.lapack.Error.InvalidArgument`: If `m` or `n` is less than 0, or if
/// `lda` is less than `max(1, m)` or `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding LAPACKE function, if available. In that case, no errors will
/// be raised even if the arguments are invalid.
pub inline fn getrf(
    order: Order,
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]i32,
    ctx: anytype,
) !i32 {
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.getrf requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.getrf requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(A)) {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.getrf not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime options.link_lapacke != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return scast(i32, ci.LAPACKE_sgetrf(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime A == f64) {
                    return scast(i32, ci.LAPACKE_dgetrf(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return scast(i32, ci.LAPACKE_cgetrf(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime Scalar(A) == f64) {
                    return scast(i32, ci.LAPACKE_zgetrf(
                        order.toCInt(),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                }
            },
            else => {},
        }
    }

    return @import("lapack/getrf.zig").getrf(order, m, n, a, lda, ipiv, ctx);
}

/// Computes the LU factorization of a general `m`-by-`n` matrix.
///
/// The `sgetrf` routine computes the LU factorization of a general `m`-by-`n`
/// matrix `A` as:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`f32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn sgetrf(
    order: Order,
    m: i32,
    n: i32,
    a: [*]f32,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf(order, m, n, a, lda, ipiv, null) catch -1;
}

/// Computes the LU factorization of a general `m`-by-`n` matrix.
///
/// The `dgetrf` routine computes the LU factorization of a general `m`-by-`n`
/// matrix `A` as:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`f64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dgetrf(
    order: Order,
    m: i32,
    n: i32,
    a: [*]f64,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf(order, m, n, a, lda, ipiv, null) catch -1;
}

/// Computes the LU factorization of a general `m`-by-`n` matrix.
///
/// The `cgetrf` routine computes the LU factorization of a general `m`-by-`n`
/// matrix `A` as:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`cf32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cgetrf(
    order: Order,
    m: i32,
    n: i32,
    a: [*]cf32,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf(order, m, n, a, lda, ipiv, null) catch -1;
}

/// Computes the LU factorization of a general `m`-by-`n` matrix.
///
/// The `zgetrf` routine computes the LU factorization of a general `m`-by-`n`
/// matrix `A` as:
///
/// ```zig
///     A = P * L * U,
/// ```
///
/// where `P` is a permutation matrix, `L` is lower triangular with unit
/// diagonal elements (lower trapezoidal if `m > n`), and `U` is upper
/// triangular (upper trapezoidal if `m < n`).
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`i32`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`i32`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`cf64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zgetrf(
    order: Order,
    m: i32,
    n: i32,
    a: [*]cf64,
    lda: i32,
    ipiv: [*]i32,
) i32 {
    return getrf(order, m, n, a, lda, ipiv, null) catch -1;
}

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
pub inline fn getrs(
    order: Order,
    transa: Transpose,
    n: i32,
    nrhs: i32,
    a: anytype,
    lda: i32,
    ipiv: [*]const i32,
    b: anytype,
    ldb: i32,
    ctx: anytype,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.lapack.getrs requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.getrs requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B) or types.isConstPointer(B))
        @compileError("zml.linalg.lapack.getrs requires b to be a mutable many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.lapack.getrs requires b's child type to be a numeric, got " ++ @typeName(B));

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B))
    {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.getrs not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (transa != .conj_no_trans and (comptime A == B and options.link_lapacke != null)) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    _ = ci.LAPACKE_sgetrs(
                        order.toCInt(),
                        transa.toChar(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    );

                    return;
                } else if (comptime A == f64) {
                    _ = ci.LAPACKE_dgetrs(
                        order.toCInt(),
                        transa.toChar(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    );

                    return;
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    _ = ci.LAPACKE_cgetrs(
                        order.toCInt(),
                        transa.toChar(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    );

                    return;
                } else if (comptime Scalar(A) == f64) {
                    _ = ci.LAPACKE_zgetrs(
                        order.toCInt(),
                        transa.toChar(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    );

                    return;
                }
            },
            else => {},
        }
    }

    return @import("lapack/getrs.zig").getrs(order, transa, n, nrhs, a, lda, ipiv, b, ldb, ctx);
}

/// Solves a system of linear equations with an LU-factored square coefficient
/// matrix, with multiple right-hand sides.
///
/// The `sgetrs` routine solves for `X` the following systems of linear
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
/// computed by `sgetrf`, `B` is an `n`-by-`nrhs` matrix of right-hand
/// sides, and `X` is an `n`-by-`nrhs` matrix of solutions.
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
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. The pivot indices as
/// returned by `getrf`.
///
/// `b` (`[*]f32`): Array, size at least `ldb * nrhs` if `order = .col_major`
/// or `ldb * n` if `order = .row_major`. On return, contains the solution
/// matrix `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = .col_major` or `max(1, nrhs)` if
/// `order = .row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn sgetrs(
    order: Order,
    transa: Transpose,
    n: i32,
    nrhs: i32,
    a: [*]const f32,
    lda: i32,
    ipiv: [*]const i32,
    b: [*]f32,
    ldb: i32,
) void {
    return getrs(order, transa, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch {};
}

/// Solves a system of linear equations with an LU-factored square coefficient
/// matrix, with multiple right-hand sides.
///
/// The `dgetrs` routine solves for `X` the following systems of linear
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
/// computed by `dgetrf`, `B` is an `n`-by-`nrhs` matrix of right-hand
/// sides, and `X` is an `n`-by-`nrhs` matrix of solutions.
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
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. The pivot indices as
/// returned by `getrf`.
///
/// `b` (`[*]f64`): Array, size at least `ldb * nrhs` if `order = .col_major`
/// or `ldb * n` if `order = .row_major`. On return, contains the solution
/// matrix `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = .col_major` or `max(1, nrhs)` if
/// `order = .row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dgetrs(
    order: Order,
    transa: Transpose,
    n: i32,
    nrhs: i32,
    a: [*]const f64,
    lda: i32,
    ipiv: [*]const i32,
    b: [*]f64,
    ldb: i32,
) void {
    return getrs(order, transa, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch {};
}

/// Solves a system of linear equations with an LU-factored square coefficient
/// matrix, with multiple right-hand sides.
///
/// The `cgetrs` routine solves for `X` the following systems of linear
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
/// computed by `cgetrf`, `B` is an `n`-by-`nrhs` matrix of right-hand
/// sides, and `X` is an `n`-by-`nrhs` matrix of solutions.
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
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. The pivot indices as
/// returned by `getrf`.
///
/// `b` (`[*]cf32`): Array, size at least `ldb * nrhs` if `order = .col_major`
/// or `ldb * n` if `order = .row_major`. On return, contains the solution
/// matrix `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = .col_major` or `max(1, nrhs)` if
/// `order = .row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cgetrs(
    order: Order,
    transa: Transpose,
    n: i32,
    nrhs: i32,
    a: [*]const cf32,
    lda: i32,
    ipiv: [*]const i32,
    b: [*]cf32,
    ldb: i32,
) void {
    return getrs(order, transa, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch {};
}

/// Solves a system of linear equations with an LU-factored square coefficient
/// matrix, with multiple right-hand sides.
///
/// The `zgetrs` routine solves for `X` the following systems of linear
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
/// computed by `zgetrf`, `B` is an `n`-by-`nrhs` matrix of right-hand
/// sides, and `X` is an `n`-by-`nrhs` matrix of solutions.
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
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. The pivot indices as
/// returned by `getrf`.
///
/// `b` (`[*]cf64`): Array, size at least `ldb * nrhs` if `order = .col_major`
/// or `ldb * n` if `order = .row_major`. On return, contains the solution
/// matrix `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = .col_major` or `max(1, nrhs)` if
/// `order = .row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zgetrs(
    order: Order,
    transa: Transpose,
    n: i32,
    nrhs: i32,
    a: [*]const cf64,
    lda: i32,
    ipiv: [*]const i32,
    b: [*]cf64,
    ldb: i32,
) void {
    return getrs(order, transa, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch {};
}

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
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.gesv requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.gesv requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B) or types.isConstPointer(B))
        @compileError("zml.linalg.lapack.gesv requires b to be a mutable many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.lapack.gesv requires b's child type to be a numeric, got " ++ @typeName(B));

    comptime if (types.isArbitraryPrecision(A) or
        types.isArbitraryPrecision(B))
    {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.gesv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (comptime A == B and options.link_lapacke != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return scast(i32, ci.LAPACKE_sgesv(
                        order.toCInt(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    ));
                } else if (comptime A == f64) {
                    return scast(i32, ci.LAPACKE_dgesv(
                        order.toCInt(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    ));
                }
            },
            .cfloat => {
                if (comptime Scalar(A) == f32) {
                    return scast(i32, ci.LAPACKE_cgesv(
                        order.toCInt(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    ));
                } else if (comptime Scalar(A) == f64) {
                    return scast(i32, ci.LAPACKE_zgesv(
                        order.toCInt(),
                        scast(c_int, n),
                        scast(c_int, nrhs),
                        a,
                        scast(c_int, lda),
                        ipiv,
                        b,
                        scast(c_int, ldb),
                    ));
                }
            },
            else => {},
        }
    }

    return @import("lapack/gesv.zig").gesv(order, n, nrhs, a, lda, ipiv, b, ldb, ctx);
}

/// Computes the solution to the system of linear equations with a square
/// coefficient matrix `A` and multiple right-hand sides.
///
/// The `sgesv` routine solves for `X` the system of linear equations
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
/// `a` (`[*]f32`): Array, size at least `lda * n`. On return, contains the LU
/// factorization of the matrix `A` as computed by `sgetrf`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. On return contains the
/// pivot indices as returned by `sgetrf`. For `1 <= i <= n`, row `i` of the
/// matrix was interchanged with row `ipiv[i - 1]`.
///
/// `b` (`[*]f32`): Array, size at least `ldb * nrhs` if `order = .col_major` or
/// `ldb * n` if `order = .row_major`. On return, contains the solution matrix
/// `X`.
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
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn sgesv(
    order: Order,
    n: i32,
    nrhs: i32,
    a: [*]f32,
    lda: i32,
    ipiv: [*]i32,
    b: [*]f32,
    ldb: i32,
) i32 {
    return gesv(order, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch -1;
}

/// Computes the solution to the system of linear equations with a square
/// coefficient matrix `A` and multiple right-hand sides.
///
/// The `dgesv` routine solves for `X` the system of linear equations
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
/// `a` (`[*]f64`): Array, size at least `lda * n`. On return, contains the LU
/// factorization of the matrix `A` as computed by `dgetrf`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. On return contains the
/// pivot indices as returned by `dgetrf`. For `1 <= i <= n`, row `i` of the
/// matrix was interchanged with row `ipiv[i - 1]`.
///
/// `b` (`[*]f64`): Array, size at least `ldb * nrhs` if `order = .col_major` or
/// `ldb * n` if `order = .row_major`. On return, contains the solution matrix
/// `X`.
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
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dgesv(
    order: Order,
    n: i32,
    nrhs: i32,
    a: [*]f64,
    lda: i32,
    ipiv: [*]i32,
    b: [*]f64,
    ldb: i32,
) i32 {
    return gesv(order, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch -1;
}

/// Computes the solution to the system of linear equations with a square
/// coefficient matrix `A` and multiple right-hand sides.
///
/// The `cgesv` routine solves for `X` the system of linear equations
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
/// `a` (`[*]cf32`): Array, size at least `lda * n`. On return, contains the LU
/// factorization of the matrix `A` as computed by `cgetrf`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. On return contains the
/// pivot indices as returned by `cgetrf`. For `1 <= i <= n`, row `i` of the
/// matrix was interchanged with row `ipiv[i - 1]`.
///
/// `b` (`[*]cf32`): Array, size at least `ldb * nrhs` if `order = .col_major`
/// or `ldb * n` if `order = .row_major`. On return, contains the solution
/// matrix `X`.
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
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cgesv(
    order: Order,
    n: i32,
    nrhs: i32,
    a: [*]cf32,
    lda: i32,
    ipiv: [*]i32,
    b: [*]cf32,
    ldb: i32,
) i32 {
    return gesv(order, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch -1;
}

/// Computes the solution to the system of linear equations with a square
/// coefficient matrix `A` and multiple right-hand sides.
///
/// The `zgesv` routine solves for `X` the system of linear equations
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
/// `a` (`[*]cf64`): Array, size at least `lda * n`. On return, contains the LU
/// factorization of the matrix `A` as computed by `zgetrf`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, n)`. On return contains the
/// pivot indices as returned by `zgetrf`. For `1 <= i <= n`, row `i` of the
/// matrix was interchanged with row `ipiv[i - 1]`.
///
/// `b` (`[*]cf64`): Array, size at least `ldb * nrhs` if `order = .col_major`
/// or `ldb * n` if `order = .row_major`. On return, contains the solution
/// matrix `X`.
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
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zgesv(
    order: Order,
    n: i32,
    nrhs: i32,
    a: [*]cf64,
    lda: i32,
    ipiv: [*]i32,
    b: [*]cf64,
    ldb: i32,
) i32 {
    return gesv(order, n, nrhs, a, lda, ipiv, b, ldb, .{}) catch -1;
}

/// Computes Cholesky factorization using a recursive algorithm.
///
/// The `potrf2` routine computes the Cholesky factorization of a real symmetric
/// or complex Hermitian positive definite matrix `A` using the recursive
/// algorithm. The factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular. This
/// is the recursive version of the algorithm. It divides the matrix into four
/// submatrices:
///
/// ```zig
///         [ A11  A12 ]
///     A = [ A21  A22 ]
/// ```
///
/// where `A11` is `n1`-by-`n1` and `A22` is `n2`-by-`n2`, with `n1 = n / 2` and
/// `n2 = n - n1`. The subroutine calls itself to factor `A11`. Update and scale
/// `A21` or `A12`, update `A22` then call itself to factor `A22`.
///
/// Signature
/// ---------
/// ```zig
/// fn potrf2(order: Order, uplo: Uplo, n: i32, a: [*]A, lda: i32, ctx: anytype) !i32
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
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `lda * n`. On return, contains the Cholesky factorization of the
/// matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Errors
/// ------
/// `linalg.lapack.Error.InvalidArgument`: If `n` is less than 0, or if `lda` is
/// less than `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding LAPACKE function, if available. In that case, no errors will
/// be raised even if the arguments are invalid.
pub inline fn potrf2(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !i32 {
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.potrf2 requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.potrf2 requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(A)) {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.potrf2 not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    // if (comptime options.link_lapacke != null) {
    //     switch (comptime types.numericType(A)) {
    //         .float => {
    //             if (comptime A == f32) {
    //                 return scast(i32, ci.LAPACKE_spotrf2(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             } else if (comptime A == f64) {
    //                 return scast(i32, ci.LAPACKE_dpotrf2(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             }
    //         },
    //         .cfloat => {
    //             if (comptime Scalar(A) == f32) {
    //                 return scast(i32, ci.LAPACKE_cpotrf2(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             } else if (comptime Scalar(A) == f64) {
    //                 return scast(i32, ci.LAPACKE_zpotrf2(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             }
    //         },
    //         else => {},
    //     }
    // }

    return @import("lapack/potrf2.zig").potrf2(order, uplo, n, a, lda, ctx);
}

/// Computes Cholesky factorization using a recursive algorithm.
///
/// The `spotrf2` routine computes the Cholesky factorization of a real
/// symmetric positive definite matrix `A` using the recursive algorithm. The
/// factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular. This
/// is the recursive version of the algorithm. It divides the matrix into four
/// submatrices:
///
/// ```zig
///         [ A11  A12 ]
///     A = [ A21  A22 ]
/// ```
///
/// where `A11` is `n1`-by-`n1` and `A22` is `n2`-by-`n2`, with `n1 = n / 2` and
/// `n2 = n - n1`. The subroutine calls itself to factor `A11`. Update and scale
/// `A21` or `A12`, update `A22` then call itself to factor `A22`.
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
/// `a` (`[*]f32`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn spotrf2(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]f32,
    lda: i32,
) i32 {
    return potrf2(order, uplo, n, a, lda, .{}) catch -1;
}

/// Computes Cholesky factorization using a recursive algorithm.
///
/// The `dpotrf2` routine computes the Cholesky factorization of a real
/// symmetric positive definite matrix `A` using the recursive algorithm. The
/// factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular. This
/// is the recursive version of the algorithm. It divides the matrix into four
/// submatrices:
///
/// ```zig
///         [ A11  A12 ]
///     A = [ A21  A22 ]
/// ```
///
/// where `A11` is `n1`-by-`n1` and `A22` is `n2`-by-`n2`, with `n1 = n / 2` and
/// `n2 = n - n1`. The subroutine calls itself to factor `A11`. Update and scale
/// `A21` or `A12`, update `A22` then call itself to factor `A22`.
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
/// `a` (`[*]f64`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dpotrf2(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]f64,
    lda: i32,
) i32 {
    return potrf2(order, uplo, n, a, lda, .{}) catch -1;
}

/// Computes Cholesky factorization using a recursive algorithm.
///
/// The `cpotrf2` routine computes the Cholesky factorization of a complex
/// hermitian positive definite matrix `A` using the recursive algorithm. The
/// factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular. This
/// is the recursive version of the algorithm. It divides the matrix into four
/// submatrices:
///
/// ```zig
///         [ A11  A12 ]
///     A = [ A21  A22 ]
/// ```
///
/// where `A11` is `n1`-by-`n1` and `A22` is `n2`-by-`n2`, with `n1 = n / 2` and
/// `n2 = n - n1`. The subroutine calls itself to factor `A11`. Update and scale
/// `A21` or `A12`, update `A22` then call itself to factor `A22`.
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
/// `a` (`[*]cf32`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cpotrf2(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]cf32,
    lda: i32,
) i32 {
    return potrf2(order, uplo, n, a, lda, .{}) catch -1;
}

/// Computes Cholesky factorization using a recursive algorithm.
///
/// The `zpotrf2` routine computes the Cholesky factorization of a complex
/// hermitian positive definite matrix `A` using the recursive algorithm. The
/// factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular. This
/// is the recursive version of the algorithm. It divides the matrix into four
/// submatrices:
///
/// ```zig
///         [ A11  A12 ]
///     A = [ A21  A22 ]
/// ```
///
/// where `A11` is `n1`-by-`n1` and `A22` is `n2`-by-`n2`, with `n1 = n / 2` and
/// `n2 = n - n1`. The subroutine calls itself to factor `A11`. Update and scale
/// `A21` or `A12`, update `A22` then call itself to factor `A22`.
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
/// `a` (`[*]f32`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zpotrf2(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]cf64,
    lda: i32,
) i32 {
    return potrf2(order, uplo, n, a, lda, .{}) catch -1;
}

/// Computes the Cholesky factorization of a symmetric or Hermitian
/// positive-definite matrix.
///
/// The `potrf` routine computes the Cholesky factorization of a real symmetric
/// or complex hermitian positive definite matrix `A`. The factorization has the
/// form:
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
/// where `U` is an upper triangular matrix and `L` is lower triangular.
///
/// Signature
/// ---------
/// ```zig
/// fn potrf(order: Order, uplo: Uplo, n: i32, a: [*]A, lda: i32, ctx: anytype) !i32
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
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size at
/// least `lda * n`. On return, contains the Cholesky factorization of the
/// matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Errors
/// ------
/// `linalg.lapack.Error.InvalidArgument`: If `n` is less than 0, or if `lda` is
/// less than `max(1, n)`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding LAPACKE function, if available. In that case, no errors will
/// be raised even if the arguments are invalid.
pub inline fn potrf(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: anytype,
    lda: i32,
    ctx: anytype,
) !i32 {
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.potrf requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.potrf requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (types.isArbitraryPrecision(A)) {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.potrf not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    // if (comptime options.link_lapacke != null) {
    //     switch (comptime types.numericType(A)) {
    //         .float => {
    //             if (comptime A == f32) {
    //                 return scast(i32, ci.LAPACKE_spotrf(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             } else if (comptime A == f64) {
    //                 return scast(i32, ci.LAPACKE_dpotrf(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             }
    //         },
    //         .cfloat => {
    //             if (comptime Scalar(A) == f32) {
    //                 return scast(i32, ci.LAPACKE_cpotrf(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             } else if (comptime Scalar(A) == f64) {
    //                 return scast(i32, ci.LAPACKE_zpotrf(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     a,
    //                     scast(c_int, lda),
    //                 ));
    //             }
    //         },
    //         else => {},
    //     }
    // }

    return @import("lapack/potrf.zig").potrf(order, uplo, n, a, lda, ctx);
}

/// Computes the Cholesky factorization of a symmetric positive-definite matrix.
///
/// The `spotrf` routine computes the Cholesky factorization of a real symmetric
/// positive definite matrix `A`. The factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular.
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
/// `a` (`[*]f32`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn spotrf(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]f32,
    lda: i32,
) i32 {
    return potrf(order, uplo, n, a, lda, .{}) catch -1;
}

/// Computes the Cholesky factorization of a symmetric positive-definite matrix.
///
/// The `dpotrf` routine computes the Cholesky factorization of a real symmetric
/// positive definite matrix `A`. The factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular.
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
/// `a` (`[*]f64`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dpotrf(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]f64,
    lda: i32,
) i32 {
    return potrf(order, uplo, n, a, lda, .{}) catch -1;
}

/// Computes the Cholesky factorization of a Hermitian positive-definite matrix.
///
/// The `cpotrf` routine computes the Cholesky factorization of a complex
/// Hermitian positive definite matrix `A`. The factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular.
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
/// `a` (`[*]cf32`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cpotrf(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]cf32,
    lda: i32,
) i32 {
    return potrf(order, uplo, n, a, lda, .{}) catch -1;
}

/// Computes the Cholesky factorization of a Hermitian positive-definite matrix.
///
/// The `zpotrf` routine computes the Cholesky factorization of a complex
/// Hermitian positive definite matrix `A`. The factorization has the form:
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
///  where `U` is an upper triangular matrix and `L` is lower triangular.
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
/// `a` (`[*]cf64`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// Returns
/// -------
/// `i32`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zpotrf(
    order: Order,
    uplo: Uplo,
    n: i32,
    a: [*]cf64,
    lda: i32,
) i32 {
    return potrf(order, uplo, n, a, lda, .{}) catch -1;
}

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
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);

    comptime if (!types.isManyPointer(A))
        @compileError("zml.linalg.lapack.potrs requires a to be a many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.potrs requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B) or types.isConstPointer(B))
        @compileError("zml.linalg.lapack.potrs requires b to be a mutable many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.lapack.potrs requires b's child type to be a numeric, got " ++ @typeName(B));

    comptime if (types.isArbitraryPrecision(A) or types.isArbitraryPrecision(B)) {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.potrs not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    // if (comptime A == B and options.link_lapacke != null) {
    //     switch (comptime types.numericType(A)) {
    //         .float => {
    //             if (comptime A == f32) {
    //                 _ = ci.LAPACKE_spotrs(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 );

    //                 return;
    //             } else if (comptime A == f64) {
    //                 std.debug.print("Calling LAPACKE_dpotrs with order: {}, uplo: {}, n: {}, nrhs: {}, lda: {}, ldb: {}\n", .{ order.toCInt(), uplo.toChar(), n, nrhs, lda, ldb });
    //                 _ = ci.LAPACKE_dpotrs(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 );

    //                 return;
    //             }
    //         },
    //         .cfloat => {
    //             if (comptime Scalar(A) == f32) {
    //                 _ = ci.LAPACKE_cpotrs(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 );

    //                 return;
    //             } else if (comptime Scalar(A) == f64) {
    //                 _ = ci.LAPACKE_zpotrs(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 );

    //                 return;
    //             }
    //         },
    //         else => {},
    //     }
    // }

    return @import("lapack/potrs.zig").potrs(order, uplo, n, nrhs, a, lda, b, ldb, ctx);
}

/// Solves a system of linear equations with a Cholesky-factored symmetric
/// positive-definite coefficient matrix.
///
/// The `spotrs` routine solves for `X` the system of linear equations
/// `A * X = B` with a symmetric positive-definite matrix `A`, given the
/// Cholesky factorization of `A`:
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
/// where `U` is an upper triangular matrix and `L` is lower triangular. The
/// system is solved with multiple right-hand sides stored in the columns of the
/// matrix `B`. Before calling this routine, you must call `spotrf` to compute
/// the Cholesky factorization of `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies how the matrix `A` has been factored:
/// - If `uplo = upper`, then the factorization is `A = U^T * U`.
/// - If `uplo = lower`, then the factorization is `A = L * L^T`.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]const f32`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]f32`): Array, size at least `ldb * nrhs` if `order = col_major`, or
/// `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn spotrs(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]const f32,
    lda: i32,
    b: [*]f32,
    ldb: i32,
    ctx: anytype,
) i32 {
    return potrs(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

/// Solves a system of linear equations with a Cholesky-factored symmetric
/// positive-definite coefficient matrix.
///
/// The `dpotrs` routine solves for `X` the system of linear equations
/// `A * X = B` with a symmetric positive-definite matrix `A`, given the
/// Cholesky factorization of `A`:
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
/// where `U` is an upper triangular matrix and `L` is lower triangular. The
/// system is solved with multiple right-hand sides stored in the columns of the
/// matrix `B`. Before calling this routine, you must call `dpotrf` to compute
/// the Cholesky factorization of `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies how the matrix `A` has been factored:
/// - If `uplo = upper`, then the factorization is `A = U^T * U`.
/// - If `uplo = lower`, then the factorization is `A = L * L^T`.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]const f64`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]f64`): Array, size at least `ldb * nrhs` if `order = col_major`, or
/// `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dpotrs(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]const f64,
    lda: i32,
    b: [*]f64,
    ldb: i32,
    ctx: anytype,
) i32 {
    return potrs(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

/// Solves a system of linear equations with a Cholesky-factored Hermitian
/// positive-definite coefficient matrix.
///
/// The `cpotrs` routine solves for `X` the system of linear equations
/// `A * X = B` with a Hermitian positive-definite matrix `A`, given the
/// Cholesky factorization of `A`:
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
/// matrix `B`. Before calling this routine, you must call `cpotrf` to compute
/// the Cholesky factorization of `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies how the matrix `A` has been factored:
/// - If `uplo = upper`, then the factorization is `A = U^H * U`.
/// - If `uplo = lower`, then the factorization is `A = L * L^H`.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]const cf32`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]cf32`): Array, size at least `ldb * nrhs` if `order = col_major`,
/// or `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cpotrs(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]const cf32,
    lda: i32,
    b: [*]cf32,
    ldb: i32,
    ctx: anytype,
) i32 {
    return potrs(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

/// Solves a system of linear equations with a Cholesky-factored Hermitian
/// positive-definite coefficient matrix.
///
/// The `zpotrs` routine solves for `X` the system of linear equations
/// `A * X = B` with a Hermitian positive-definite matrix `A`, given the
/// Cholesky factorization of `A`:
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
/// matrix `B`. Before calling this routine, you must call `zpotrf` to compute
/// the Cholesky factorization of `A`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies how the matrix `A` has been factored:
/// - If `uplo = upper`, then the factorization is `A = U^H * U`.
/// - If `uplo = lower`, then the factorization is `A = L * L^H`.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]const cf64`): Array, size at least `lda * n`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]cf64`): Array, size at least `ldb * nrhs` if `order = col_major`,
/// or `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zpotrs(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]const cf64,
    lda: i32,
    b: [*]cf64,
    ldb: i32,
    ctx: anytype,
) i32 {
    return potrs(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

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
pub inline fn posv(
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
    comptime var A: type = @TypeOf(a);
    comptime var B: type = @TypeOf(b);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.posv requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.posv requires a's child type to be a numeric, got " ++ @typeName(A));

    comptime if (!types.isManyPointer(B) or types.isConstPointer(B))
        @compileError("zml.linalg.lapack.posv requires b to be a mutable many-item pointer, got " ++ @typeName(B));

    B = types.Child(B);

    comptime if (!types.isNumeric(B))
        @compileError("zml.linalg.lapack.posv requires b's child type to be a numeric, got " ++ @typeName(B));

    comptime if (types.isArbitraryPrecision(A) or types.isArbitraryPrecision(B)) {
        // When implemented, expand if
        @compileError("zml.linalg.lapack.posv not implemented for arbitrary precision types yet");
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    // if (comptime A == B and options.link_lapacke != null) {
    //     switch (comptime types.numericType(A)) {
    //         .float => {
    //             if (comptime A == f32) {
    //                 return scast(i32, ci.LAPACKE_sposv(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 ));
    //             } else if (comptime A == f64) {
    //                 return scast(i32, ci.LAPACKE_dposv(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 ));
    //             }
    //         },
    //         .cfloat => {
    //             if (comptime Scalar(A) == f32) {
    //                 return scast(i32, ci.LAPACKE_cposv(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 ));
    //             } else if (comptime Scalar(A) == f64) {
    //                 return scast(i32, ci.LAPACKE_zposv(
    //                     order.toCInt(),
    //                     uplo.toChar(),
    //                     scast(c_int, n),
    //                     scast(c_int, nrhs),
    //                     a,
    //                     scast(c_int, lda),
    //                     b,
    //                     scast(c_int, ldb),
    //                 ));
    //             }
    //         },
    //         else => {},
    //     }
    // }

    return @import("lapack/posv.zig").posv(order, uplo, n, nrhs, a, lda, b, ldb, ctx);
}

/// Computes the solution to the system of linear equations with a symmetric
/// positive-definite coefficient matrix `A` and multiple right-hand sides.
///
/// The `sposv` routine solves for `X` the real system of linear equations
/// `A * X = B`, where `A` is an `n`-by-`n` symmetric positive-definite matrix,
/// the columns of matrix `B` are individual right-hand sides, and the columns
/// of `X` are the corresponding solutions. The Cholesky decomposition is used
/// to factor `A` as:
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
/// where `U` is an upper triangular matrix and `L` is lower triangular. The
/// factored form of `A` is then used to solve the system of equations
/// `A * X = B`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies which part of the matrix `A` is stored, and
/// which factorization is computed:
/// - If `uplo = upper`, then the upper triangular part of `A` is stored, and
/// the factorization is `A = U^T * U` is computed.
/// - If `uplo = lower`, then the lower triangular part of `A` is stored, and
/// the factorization is `A = L * L^T` is computed.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]f32`): Array, size at least `ldb * nrhs` if `order = col_major`, or
/// `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a` and `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn sposv(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]f32,
    lda: i32,
    b: [*]f32,
    ldb: i32,
    ctx: anytype,
) i32 {
    return posv(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

/// Computes the solution to the system of linear equations with a symmetric
/// positive-definite coefficient matrix `A` and multiple right-hand sides.
///
/// The `dposv` routine solves for `X` the real system of linear equations
/// `A * X = B`, where `A` is an `n`-by-`n` symmetric positive-definite matrix,
/// the columns of matrix `B` are individual right-hand sides, and the columns
/// of `X` are the corresponding solutions. The Cholesky decomposition is used
/// to factor `A` as:
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
/// where `U` is an upper triangular matrix and `L` is lower triangular. The
/// factored form of `A` is then used to solve the system of equations
/// `A * X = B`.
///
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies which part of the matrix `A` is stored, and
/// which factorization is computed:
/// - If `uplo = upper`, then the upper triangular part of `A` is stored, and
/// the factorization is `A = U^T * U` is computed.
/// - If `uplo = lower`, then the lower triangular part of `A` is stored, and
/// the factorization is `A = L * L^T` is computed.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]f64`): Array, size at least `ldb * nrhs` if `order = col_major`, or
/// `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a` and `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dposv(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]f64,
    lda: i32,
    b: [*]f64,
    ldb: i32,
    ctx: anytype,
) i32 {
    return posv(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

/// Computes the solution to the system of linear equations with a Hermitian
/// positive-definite coefficient matrix `A` and multiple right-hand sides.
///
/// The `cposv` routine solves for `X` the complex system of linear equations
/// `A * X = B`, where `A` is an `n`-by-`n` Hermitian positive-definite matrix,
/// the columns of matrix `B` are individual right-hand sides, and the columns
/// of `X` are the corresponding solutions. The Cholesky decomposition is used
/// to factor `A` as:
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
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies which part of the matrix `A` is stored, and
/// which factorization is computed:
/// - If `uplo = upper`, then the upper triangular part of `A` is stored, and
/// the factorization is `A = U^H * U` is computed.
/// - If `uplo = lower`, then the lower triangular part of `A` is stored, and
/// the factorization is `A = L * L^H` is computed.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]cf32`): Array, size at least `ldb * nrhs` if `order = col_major`,
/// or `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a` and `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cposv(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]cf32,
    lda: i32,
    b: [*]cf32,
    ldb: i32,
    ctx: anytype,
) i32 {
    return posv(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

/// Computes the solution to the system of linear equations with a Hermitian
/// positive-definite coefficient matrix `A` and multiple right-hand sides.
///
/// The `zposv` routine solves for `X` the complex system of linear equations
/// `A * X = B`, where `A` is an `n`-by-`n` Hermitian positive-definite matrix,
/// the columns of matrix `B` are individual right-hand sides, and the columns
/// of `X` are the corresponding solutions. The Cholesky decomposition is used
/// to factor `A` as:
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
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `uplo` (`Uplo`): Specifies which part of the matrix `A` is stored, and
/// which factorization is computed:
/// - If `uplo = upper`, then the upper triangular part of `A` is stored, and
/// the factorization is `A = U^H * U` is computed.
/// - If `uplo = lower`, then the lower triangular part of `A` is stored, and
/// the factorization is `A = L * L^H` is computed.
///
/// `n` (`i32`): The order of the matrix `A`. Must be greater than or equal to
/// 0.
///
/// `nrhs` (`i32`): The number of right-hand sides, i.e., the number of
/// columns of the matrix `B`. Must be greater than or equal to 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * n`. On return, contains the
/// Cholesky factorization of the matrix `A`.
///
/// `lda` (`i32`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, n)`.
///
/// `b` (`[*]cf64`): Array, size at least `ldb * nrhs` if `order = col_major`,
/// or `ldb * n` if `order = row_major`. On return, contains the solution matrix
/// `X`.
///
/// `ldb` (`i32`): The leading dimension of the array `b`. Must be greater
/// than or equal to `max(1, n)` if `order = col_major`, or `max(1, nrhs)` if
/// `order = row_major`.
///
/// Returns
/// -------
/// `void`: The result is stored in `a` and `b`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zposv(
    order: Order,
    uplo: Uplo,
    n: i32,
    nrhs: i32,
    a: [*]cf64,
    lda: i32,
    b: [*]cf64,
    ldb: i32,
    ctx: anytype,
) i32 {
    return posv(order, uplo, n, nrhs, a, lda, b, ldb, ctx) catch {};
}

pub const Error = error{
    InvalidArgument,
};
