const std = @import("std");

const options = @import("options");

const float = @import("../../../float.zig");
const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

/// Performs a rank-1 update of a Hermitian matrix defined as:
///
/// ```zig
/// A = αx xᴴ + A,
/// ```
///
/// where `α` is a real numeric, `x` is an `n`-element vector, and `A` is an
/// `n × n` Hermitian matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.her(layout: matrix.Layout, uplo: matrix.Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, a: [*]A, lda: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the Hermitian matrix `A` is used.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `α`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * n`. On
///   return, contains the result of the operation.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is 0, or if `lda` is less
///   than `max(1, n)`.
pub fn her(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    a: anytype,
    lda: usize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var A: type = @TypeOf(a);

    comptime if (!meta.isNumeric(Al) or meta.isComplex(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.blas.her: alpha must be a real numeric, x must be a many-item pointer to numerics, and a must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ta: " ++ @typeName(A) ++ "\n");

    X = meta.Child(X);
    A = meta.Child(A);

    if (lda < int.max(1, n) or incx == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const noconj = layout == .col_major;

    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    if (noconj)
        k_her(eff_uplo, n, alpha, x, incx, a, lda, true)
    else
        k_her(eff_uplo, n, alpha, x, incx, a, lda, false);
}

fn k_her(uplo: matrix.Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, a: anytype, lda: usize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Set up the start point in `x`.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (uplo == .upper) {
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(if (comptime noconj) X else numeric.Conj(X), numeric.Mul(Al, if (comptime noconj) numeric.Conj(X) else X), A)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;

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

            var j: usize = tile_i;
            var jx: isize = kx + numeric.cast(isize, tile_i) * incx;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    const temp = numeric.mul(
                        alpha,
                        if (comptime noconj)
                            numeric.conj(x[numeric.cast(usize, jx)])
                        else
                            x[numeric.cast(usize, jx)],
                    );

                    if (j < tile_i + b_len) {
                        var i: usize = 0;
                        while (i < ((j - tile_i) / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    if (comptime noconj)
                                        px[i + u]
                                    else
                                        numeric.conj(px[i + u]),
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < j - tile_i) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                if (comptime noconj)
                                    px[i]
                                else
                                    numeric.conj(px[i]),
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }

                        // a[j + j * lda] = re(a[j + j * lda]) + re(px[j - tile_i] * temp)
                        numeric.fmaInto(
                            &a[j + j * lda],
                            if (comptime noconj)
                                numeric.neg(numeric.im(px[j - tile_i]))
                            else
                                numeric.im(px[j - tile_i]),
                            numeric.im(temp),
                            numeric.fma(
                                numeric.re(px[j - tile_i]),
                                numeric.re(temp),
                                numeric.re(a[j + j * lda]),
                            ),
                        );
                    } else {
                        var i: usize = 0;
                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    if (comptime noconj)
                                        px[i + u]
                                    else
                                        numeric.conj(px[i + u]),
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                if (comptime noconj)
                                    px[i]
                                else
                                    numeric.conj(px[i]),
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }
                    }
                } else {
                    if (j < tile_i + b_len) {
                        // a[j + j * lda] = re(a[j + j * lda])
                        numeric.set(
                            &a[j + j * lda],
                            numeric.re(a[j + j * lda]),
                        );
                    }
                }

                jx += incx;
            }
        }
    } else {
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(if (comptime noconj) X else numeric.Conj(X), numeric.Mul(Al, if (comptime noconj) numeric.Conj(X) else X), A)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < n) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, n - tile_i);
            var local_x: [tile_size]X = undefined;

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

            var j: usize = 0;
            var jx: isize = kx;
            while (j < tile_i + b_len) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                    const temp = numeric.mul(
                        alpha,
                        if (comptime noconj)
                            numeric.conj(x[numeric.cast(usize, jx)])
                        else
                            x[numeric.cast(usize, jx)],
                    );

                    if (j >= tile_i) {
                        // a[j + j * lda] = re(a[j + j * lda]) + re(conj(px[j - tile_i]) * temp)
                        numeric.fmaInto(
                            &a[j + j * lda],
                            if (comptime !noconj)
                                numeric.im(px[j - tile_i])
                            else
                                numeric.neg(numeric.im(px[j - tile_i])),
                            numeric.im(temp),
                            numeric.fma(
                                numeric.re(px[j - tile_i]),
                                numeric.re(temp),
                                numeric.re(a[j + j * lda]),
                            ),
                        );

                        var i: usize = j - tile_i + 1;
                        while (i < b_len and i % unroll != 0) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                if (comptime noconj)
                                    px[i]
                                else
                                    numeric.conj(px[i]),
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }

                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    if (comptime noconj)
                                        px[i + u]
                                    else
                                        numeric.conj(px[i + u]),
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                if (comptime noconj)
                                    px[i]
                                else
                                    numeric.conj(px[i]),
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }
                    } else {
                        var i: usize = 0;
                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += px[i + u] * temp
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    if (comptime noconj)
                                        px[i + u]
                                    else
                                        numeric.conj(px[i + u]),
                                    temp,
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += px[i] * temp
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                if (comptime noconj)
                                    px[i]
                                else
                                    numeric.conj(px[i]),
                                temp,
                                a[tile_i + i + j * lda],
                            );
                        }
                    }
                } else {
                    if (j >= tile_i) {
                        // a[j + j * lda] = re(a[j + j * lda])
                        numeric.set(
                            &a[j + j * lda],
                            numeric.re(a[j + j * lda]),
                        );
                    }
                }

                jx += incx;
            }
        }
    }
}
