const std = @import("std");

const options = @import("options");

const float = @import("../../../float.zig");
const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

/// Computes a matrix-vector product with a Hermitian matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are `n`-element vectors,
/// and `A` is an `n`-by-`n` Hermitian matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.hemv(layout: matrix.Layout, uplo: matrix.Uplo, n: usize, alpha: Al, a: [*]const A, lda: usize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the Hermitian matrix `A` is used.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `y` need not be set on input.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`. On return, contains the result of the
///   operation.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, n)`, or
///   if `incx` or `incy` is 0.
pub fn hemv(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.hemv: alpha and beta must be numerics, a and x must be many-item pointers to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);
    Y = meta.Child(Y);

    if (lda < int.max(1, n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const noconj = layout == .col_major;

    // Quick return if possible.
    if (n == 0)
        return;

    if (numeric.eq(alpha, 0)) {
        if (numeric.ne(beta, 1))
            linalg.blas.scal(n, beta, y, incy) catch unreachable;

        return;
    }

    if (noconj)
        k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, true)
    else
        k_hemv(eff_uplo, n, alpha, a, lda, x, incx, beta, y, incy, false);
}

fn k_hemv(uplo: matrix.Uplo, n: usize, alpha: anytype, a: anytype, lda: usize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Set up the start points in x and y.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;

    // First form  y = beta * y.
    if (numeric.ne(beta, 1))
        linalg.blas.scal(n, beta, y, incy) catch unreachable;

    if (numeric.eq(alpha, 0))
        return;

    if (uplo == .upper) {
        // Form  y  when upper triangle of A is stored.
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
            std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X))) orelse 2,
        );
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(Y) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;
            var local_y: [tile_size]Y = undefined;

            const px = if (incx == 1)
                x + tile_i
            else blk: {
                linalg.blas.copy(
                    b_len,
                    x + numeric.cast(usize, if (incx > 0)
                        numeric.cast(isize, tile_i) * incx
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incx),
                    incx,
                    @as([*]X, &local_x),
                    1,
                ) catch unreachable;

                break :blk @as([*]const X, &local_x);
            };

            const py = if (incy == 1)
                y + tile_i
            else blk: {
                linalg.blas.copy(
                    b_len,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                    @as([*]Y, &local_y),
                    1,
                ) catch unreachable;

                break :blk @as([*]Y, &local_y);
            };

            var j: usize = tile_i;
            var jx: isize = kx + numeric.cast(isize, tile_i) * incx;
            var jy: isize = ky + numeric.cast(isize, tile_i) * incy;
            while (j < n) : (j += 1) {
                if (j < tile_i + b_len) {
                    // temp1 = alpha * px[j - tile_i]
                    const temp1 = numeric.mul(alpha, px[j - tile_i]);

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = 0;
                    while (i < ((j - tile_i) / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fmaInto(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fmaInto(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.addInto(&temp2, temp2, sums[u]);
                    }

                    while (i < j - tile_i) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fmaInto(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fmaInto(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // py[j - tile_i] += temp1 * re(a[j + j * lda]) + alpha * temp2
                    numeric.fmaInto(
                        &py[j - tile_i],
                        temp1,
                        numeric.re(a[j + j * lda]),
                        numeric.fma(
                            alpha,
                            temp2,
                            py[j - tile_i],
                        ),
                    );
                } else {
                    // temp1 = alpha * x[jx]
                    const temp1 = numeric.mul(alpha, x[numeric.cast(usize, jx)]);

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = 0;
                    while (i + unroll <= b_len) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fmaInto(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fmaInto(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.addInto(&temp2, temp2, sums[u]);
                    }

                    while (i < b_len) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fmaInto(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) + px[i]
                        numeric.fmaInto(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // y[jy] += alpha * temp2
                    numeric.fmaInto(
                        &y[numeric.cast(usize, jy)],
                        alpha,
                        temp2,
                        y[numeric.cast(usize, jy)],
                    );
                }

                jx += incx;
                jy += incy;
            }

            if (incy != 1) {
                linalg.blas.copy(
                    b_len,
                    @as([*]Y, &local_y),
                    1,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                ) catch unreachable;
            }
        }
    } else {
        // Form  y  when lower triangle of A is stored.
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2,
            std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X))) orelse 2,
        );
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(Y) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;
            var local_y: [tile_size]Y = undefined;

            const px = if (incx == 1)
                x + tile_i
            else blk: {
                linalg.blas.copy(
                    b_len,
                    x + numeric.cast(usize, if (incx > 0)
                        numeric.cast(isize, tile_i) * incx
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incx),
                    incx,
                    @as([*]X, &local_x),
                    1,
                ) catch unreachable;

                break :blk @as([*]const X, &local_x);
            };

            const py = if (incy == 1)
                y + tile_i
            else blk: {
                linalg.blas.copy(
                    b_len,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                    @as([*]Y, &local_y),
                    1,
                ) catch unreachable;

                break :blk @as([*]Y, &local_y);
            };

            var j: usize = 0;
            var jx: isize = kx;
            var jy: isize = ky;
            while (j < tile_i + b_len) : (j += 1) {
                if (j >= tile_i) {
                    // temp1 = alpha * px[j - tile_i]
                    const temp1 = numeric.mul(
                        alpha,
                        px[j - tile_i],
                    );

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    // py[j - tile_i] += temp1 * re(a[j + j * lda])
                    numeric.fmaInto(
                        &py[j - tile_i],
                        temp1,
                        numeric.re(a[j + j * lda]),
                        py[j - tile_i],
                    );

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = j - tile_i + 1;
                    while (i < b_len and i % unroll != 0) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fmaInto(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fmaInto(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    while (i < (b_len / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fmaInto(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fmaInto(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.addInto(&temp2, temp2, sums[u]);
                    }

                    while (i < b_len) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fmaInto(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fmaInto(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // py[j - tile_i] += alpha * temp2
                    numeric.fmaInto(
                        &py[j - tile_i],
                        alpha,
                        temp2,
                        py[j - tile_i],
                    );
                } else {
                    // temp1 = alpha * x[jx]
                    const temp1 = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    var temp2 = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)));

                    var sums: [unroll]meta.Accumulator(numeric.Mul(if (noconj) A else numeric.Conj(A), X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) A else numeric.Conj(A), X)))} ** unroll;

                    var i: usize = 0;
                    while (i < (b_len / unroll) * unroll) : (i += unroll) {
                        inline for (0..unroll) |u| {
                            // py[i + u] += temp1 * a[tile_i + i + u + j * lda]
                            numeric.fmaInto(
                                &py[i + u],
                                temp1,
                                if (comptime noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                py[i + u],
                            );

                            // sums[u] += conj(a[tile_i + i + u + j * lda]) * px[i + u]
                            numeric.fmaInto(
                                &sums[u],
                                if (comptime !noconj)
                                    a[tile_i + i + u + j * lda]
                                else
                                    numeric.conj(a[tile_i + i + u + j * lda]),
                                px[i + u],
                                sums[u],
                            );
                        }
                    }

                    inline for (0..unroll) |u| {
                        numeric.addInto(&temp2, temp2, sums[u]);
                    }

                    while (i < b_len) : (i += 1) {
                        // py[i] += temp1 * a[tile_i + i + j * lda]
                        numeric.fmaInto(
                            &py[i],
                            temp1,
                            if (comptime noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            py[i],
                        );

                        // temp2 += conj(a[tile_i + i + j * lda]) * px[i]
                        numeric.fmaInto(
                            &temp2,
                            if (comptime !noconj)
                                a[tile_i + i + j * lda]
                            else
                                numeric.conj(a[tile_i + i + j * lda]),
                            px[i],
                            temp2,
                        );
                    }

                    // y[jy] += alpha * temp2
                    numeric.fmaInto(
                        &y[numeric.cast(usize, jy)],
                        alpha,
                        temp2,
                        y[numeric.cast(usize, jy)],
                    );
                }

                jx += incx;
                jy += incy;
            }

            if (incy != 1) {
                linalg.blas.copy(
                    b_len,
                    @as([*]Y, &local_y),
                    1,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                ) catch unreachable;
            }
        }
    }

    return;
}
