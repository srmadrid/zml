const std = @import("std");

const options = @import("options");

const float = @import("../../../float.zig");
const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

/// Performs a rank-2 update of a symmetric matrix defined as:
///
/// ```zig
/// A = αx yᵀ + αy xᵀ + A,
/// ```
///
/// where `α` is a numeric, `x` and `y` are `n`-element vectors, and `A` is an
/// `n × n` symmetric matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.syr2(layout: matrix.Layout, uplo: matrix.Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, a: [*]A, lda: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the upper or lower triangular
///   part of the symmetric matrix `A` is used.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * n`. On
///   return, contains the result of the operation.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0, or if `lda`
///   is less than `max(1, n)`.
pub fn syr2(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: usize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or !meta.isNumeric(meta.Child(Y)) or
        !meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.blas.syr2: alpha must be a real numeric, x and y must be many-item pointers to numerics, and a must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\ta: " ++ @typeName(A) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);
    A = meta.Child(A);

    if (lda < int.max(1, n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();

    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    k_syr2(eff_uplo, n, alpha, x, incx, y, incy, a, lda);
}

fn k_syr2(uplo: matrix.Uplo, n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize, a: anytype, lda: usize) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Set up the start points in `x` and `y`.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;

    if (uplo == .upper) {
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, Y), X, A)) orelse 2,
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), Y, A)) orelse 2,
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

                break :blk @as([*]const Y, &local_y);
            };

            var j: usize = tile_i;
            var jx: isize = kx + numeric.cast(isize, tile_i) * incx;
            var jy: isize = ky + numeric.cast(isize, tile_i) * incy;
            while (j < n) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0) or numeric.ne(y[numeric.cast(usize, jy)], 0)) {
                    // temp1 = alpha * y[jy]
                    const temp1 = numeric.mul(
                        alpha,
                        y[numeric.cast(usize, jy)],
                    );

                    // temp2 = alpha * x[jx]
                    const temp2 = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    if (j < tile_i + b_len) {
                        var i: usize = 0;
                        while (i < ((j - tile_i) / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += temp1 * px[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp1,
                                    px[i + u],
                                    a[tile_i + i + u + j * lda],
                                );

                                // a[tile_i + i + u + j * lda] += temp2 * py[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp2,
                                    py[i + u],
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < j - tile_i) : (i += 1) {
                            // a[tile_i + i + j * lda] += temp1 * px[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp1,
                                px[i],
                                a[tile_i + i + j * lda],
                            );

                            // a[tile_i + i + j * lda] += temp2 * py[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp2,
                                py[i],
                                a[tile_i + i + j * lda],
                            );
                        }

                        // a[j + j * lda] += px[j - tile_i] * temp1
                        numeric.fmaInto(
                            &a[j + j * lda],
                            px[j - tile_i],
                            temp1,
                            a[j + j * lda],
                        );

                        // a[j + j * lda] += py[j - tile_i] * temp2
                        numeric.fmaInto(
                            &a[j + j * lda],
                            py[j - tile_i],
                            temp2,
                            a[j + j * lda],
                        );
                    } else {
                        var i: usize = 0;
                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += temp1 * px[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp1,
                                    px[i + u],
                                    a[tile_i + i + u + j * lda],
                                );

                                // a[tile_i + i + u + j * lda] += temp2 * py[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp2,
                                    py[i + u],
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += temp1 * px[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp1,
                                px[i],
                                a[tile_i + i + j * lda],
                            );

                            // a[tile_i + i + j * lda] += temp2 * py[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp2,
                                py[i],
                                a[tile_i + i + j * lda],
                            );
                        }
                    }
                }

                jx += incx;
                jy += incy;
            }
        }
    } else {
        const unroll = 2 * int.min(
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, Y), X, A)) orelse 2,
            std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), Y, A)) orelse 2,
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

                break :blk @as([*]const Y, &local_y);
            };

            var j: usize = 0;
            var jx: isize = kx;
            var jy: isize = ky;
            while (j < tile_i + b_len) : (j += 1) {
                if (numeric.ne(x[numeric.cast(usize, jx)], 0) or numeric.ne(y[numeric.cast(usize, jy)], 0)) {
                    // temp1 = alpha * y[jy]
                    const temp1 = numeric.mul(
                        alpha,
                        y[numeric.cast(usize, jy)],
                    );

                    // temp2 = alpha * x[jx]
                    const temp2 = numeric.mul(
                        alpha,
                        x[numeric.cast(usize, jx)],
                    );

                    if (j >= tile_i) {
                        // a[j + j * lda] += px[j - tile_i]) * temp1
                        numeric.fmaInto(
                            &a[j + j * lda],
                            px[j - tile_i],
                            temp1,
                            a[j + j * lda],
                        );

                        // a[j + j * lda] += py[j - tile_i]) * temp2
                        numeric.fmaInto(
                            &a[j + j * lda],
                            py[j - tile_i],
                            temp2,
                            a[j + j * lda],
                        );

                        var i: usize = j - tile_i + 1;
                        while (i < b_len and i % unroll != 0) : (i += 1) {
                            // a[tile_i + i + j * lda] += temp1 * px[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp1,
                                px[i],
                                a[tile_i + i + j * lda],
                            );

                            // a[tile_i + i + j * lda] += temp2 * py[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp2,
                                py[i],
                                a[tile_i + i + j * lda],
                            );
                        }

                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += temp1 * px[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp1,
                                    px[i + u],
                                    a[tile_i + i + u + j * lda],
                                );

                                // a[tile_i + i + u + j * lda] += temp2 * py[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp2,
                                    py[i + u],
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += temp1 * px[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp1,
                                px[i],
                                a[tile_i + i + j * lda],
                            );

                            // a[tile_i + i + j * lda] += temp2 * py[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp2,
                                py[i],
                                a[tile_i + i + j * lda],
                            );
                        }
                    } else {
                        var i: usize = 0;
                        while (i < (b_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // a[tile_i + i + u + j * lda] += temp1 * px[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp1,
                                    px[i + u],
                                    a[tile_i + i + u + j * lda],
                                );

                                // a[tile_i + i + u + j * lda] += temp2 * py[i + u]
                                numeric.fmaInto(
                                    &a[tile_i + i + u + j * lda],
                                    temp2,
                                    py[i + u],
                                    a[tile_i + i + u + j * lda],
                                );
                            }
                        }

                        while (i < b_len) : (i += 1) {
                            // a[tile_i + i + j * lda] += temp1 * px[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp1,
                                px[i],
                                a[tile_i + i + j * lda],
                            );

                            // a[tile_i + i + j * lda] += temp2 * py[i]
                            numeric.fmaInto(
                                &a[tile_i + i + j * lda],
                                temp2,
                                py[i],
                                a[tile_i + i + j * lda],
                            );
                        }
                    }
                }

                jx += incx;
                jy += incy;
            }
        }
    }
}
