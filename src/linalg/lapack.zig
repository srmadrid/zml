const std = @import("std");
const opts = @import("options");

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

const Order = linalg.Order;
const Transpose = linalg.Transpose;
const Uplo = linalg.Uplo;
const Diag = linalg.Diag;
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

    if (comptime opts.link_lapacke != null) {
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
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size
/// at least `lda * n` if `order = .col_major` or `lda * mm` if
/// `order = .row_major`, where `mm` is not less than the maximum of
/// `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`isize`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`isize`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`isize`): The increment for the elements of `ipiv`. If `incx` is
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
    n: isize,
    a: anytype,
    lda: isize,
    k1: isize,
    k2: isize,
    ipiv: [*]const i32,
    incx: isize,
) !void {
    comptime var A: type = @TypeOf(a);

    comptime if (!types.isManyPointer(A) or types.isConstPointer(A))
        @compileError("zml.linalg.lapack.laswp requires a to be a mutable many-item pointer, got " ++ @typeName(A));

    A = types.Child(A);

    comptime if (!types.isNumeric(A))
        @compileError("zml.linalg.lapack.laswp requires a's child type to be a numeric, got " ++ @typeName(A));

    if (comptime opts.link_lapacke != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    const result = ci.LAPACKE_slaswp(@intFromEnum(order), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

                    if (result < 0) {
                        if (result == -1011)
                            return std.mem.Allocator.Error.OutOfMemory;

                        return Error.InvalidArgument;
                    }

                    return;
                } else if (comptime A == f64) {
                    const result = ci.LAPACKE_dlaswp(@intFromEnum(order), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

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
                    const result = ci.LAPACKE_claswp(@intFromEnum(order), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

                    if (result < 0) {
                        if (result == -1011)
                            return std.mem.Allocator.Error.OutOfMemory;

                        return Error.InvalidArgument;
                    }

                    return;
                } else if (comptime Scalar(A) == f64) {
                    const result = ci.LAPACKE_zlaswp(@intFromEnum(order), scast(c_int, n), a, scast(c_int, lda), scast(c_int, k1), scast(c_int, k2), ipiv, scast(c_int, incx));

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
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`isize`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`isize`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`isize`): The increment for the elements of `ipiv`. If `incx` is
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
    n: isize,
    a: [*]f32,
    lda: isize,
    k1: isize,
    k2: isize,
    ipiv: [*]const i32,
    incx: isize,
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
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`isize`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`isize`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`isize`): The increment for the elements of `ipiv`. If `incx` is
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
    n: isize,
    a: [*]f64,
    lda: isize,
    k1: isize,
    k2: isize,
    ipiv: [*]const i32,
    incx: isize,
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
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`isize`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`isize`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`isize`): The increment for the elements of `ipiv`. If `incx` is
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
    n: isize,
    a: [*]cf32,
    lda: isize,
    k1: isize,
    k2: isize,
    ipiv: [*]const i32,
    incx: isize,
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
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * mm` if `order = .row_major`, where `mm` is not less than the maximum
/// of `ipiv[k1 - 1 + j  *abs(incx)]`, for `0 <= j < k2 - k1`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `k1` (`isize`): The first element of `ipiv` for which a row interchange will
/// be done.
///
/// `k2` (`isize`): The last element of `ipiv` for which a row interchange will
/// be done.
///
/// `ipiv` (`[*]const i32`): Array, size at least `k1 + (k2 - k1) * abs(incx)`.
/// The vector of pivot indices. Only the elements in positions `k1` through
/// `k2` of ipiv are accessed. `ipiv[k] = l` implies rows `k` and `l` are to be
/// interchanged.
///
/// `incx` (`isize`): The increment for the elements of `ipiv`. If `incx` is
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
    n: isize,
    a: [*]cf64,
    lda: isize,
    k1: isize,
    k2: isize,
    ipiv: [*]const i32,
    incx: isize,
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
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size
/// at least `lda * n` if `order = .col_major` or `lda * m` if
/// `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
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
    m: isize,
    n: isize,
    a: anytype,
    lda: isize,
    ipiv: [*]i32,
    ctx: anytype,
) !isize {
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

    if (comptime opts.link_lapacke != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return scast(isize, ci.LAPACKE_sgetrf2(
                        @intFromEnum(order),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime A == f64) {
                    return scast(isize, ci.LAPACKE_dgetrf2(
                        @intFromEnum(order),
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
                    return scast(isize, ci.LAPACKE_cgetrf2(
                        @intFromEnum(order),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime Scalar(A) == f64) {
                    return scast(isize, ci.LAPACKE_zgetrf2(
                        @intFromEnum(order),
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn sgetrf2(
    order: Order,
    m: isize,
    n: isize,
    a: [*]f32,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf2(order, m, n, a, lda, ipiv, null) catch {};
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]f64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dgetrf2(
    order: Order,
    m: isize,
    n: isize,
    a: [*]f64,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf2(order, m, n, a, lda, ipiv, null) catch {};
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cgetrf2(
    order: Order,
    m: isize,
    n: isize,
    a: [*]cf32,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf2(order, m, n, a, lda, ipiv, null) catch {};
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`[*]cf64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]const i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zgetrf2(
    order: Order,
    m: isize,
    n: isize,
    a: [*]cf64,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf2(order, m, n, a, lda, ipiv, null) catch {};
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
/// Parameters
/// ----------
/// `order` (`Order`): Specifies whether two-dimensional array storage is
/// row-major or column-major.
///
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (mutable many-item pointer to `bool`, `int`, `float`, `cfloat`,
/// `integer`, `rational`, `real`, `complex` or `expression`): Array, size
/// at least `lda * n` if `order = .col_major` or `lda * m` if
/// `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
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
    m: isize,
    n: isize,
    a: anytype,
    lda: isize,
    ipiv: [*]i32,
    ctx: anytype,
) !isize {
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

    if (comptime opts.link_lapacke != null) {
        switch (comptime types.numericType(A)) {
            .float => {
                if (comptime A == f32) {
                    return scast(isize, ci.LAPACKE_sgetrf(
                        @intFromEnum(order),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime A == f64) {
                    return scast(isize, ci.LAPACKE_dgetrf(
                        @intFromEnum(order),
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
                    return scast(isize, ci.LAPACKE_cgetrf(
                        @intFromEnum(order),
                        scast(c_int, m),
                        scast(c_int, n),
                        a,
                        scast(c_int, lda),
                        ipiv,
                    ));
                } else if (comptime Scalar(A) == f64) {
                    return scast(isize, ci.LAPACKE_zgetrf(
                        @intFromEnum(order),
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`f32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn sgetrf(
    order: Order,
    m: isize,
    n: isize,
    a: [*]f32,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf(order, m, n, a, lda, ipiv, null) catch {};
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`f64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn dgetrf(
    order: Order,
    m: isize,
    n: isize,
    a: [*]f64,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf(order, m, n, a, lda, ipiv, null) catch {};
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`cf32`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn cgetrf(
    order: Order,
    m: isize,
    n: isize,
    a: [*]cf32,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf(order, m, n, a, lda, ipiv, null) catch {};
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
/// `m` (`isize`): The number of rows of the matrix `A`. Must be greater than or
/// equal to 0.
///
/// `n` (`isize`): The number of columns of the matrix `A`. Must be greater than
/// or equal to 0.
///
/// `a` (`cf64`): Array, size at least `lda * n` if `order = .col_major` or
/// `lda * m` if `order = .row_major`.
///
/// `lda` (`isize`): The leading dimension of the array `a`. Must be grater than
/// or equal to `max(1, m)` if `order = .col_major` or `max(1, n)` if
/// `order = .row_major`.
///
/// `ipiv` (`[*]i32`): Array, size at least `max(1, min(m, n))`. On return
/// contains the pivot indices; for `1 <= i <= min(m, n)`, row `i` of the matrix
/// was interchanged with row `ipiv[i - 1]`.
///
/// Returns
/// -------
/// `isize`: 0 if successful, or `i` if `u11` is exactly zero. The result is
/// stored in `a` and `ipiv`.
///
/// Notes
/// -----
/// If the `link_cblas` option is not `null`, the function will call the
/// corresponding LAPACKE function.
pub inline fn zgetrf(
    order: Order,
    m: isize,
    n: isize,
    a: [*]cf64,
    lda: isize,
    ipiv: [*]i32,
) !isize {
    return getrf(order, m, n, a, lda, ipiv, null) catch {};
}

pub const Error = error{
    InvalidArgument,
};
