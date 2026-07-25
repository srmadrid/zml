const std = @import("std");

const options = @import("options");

const float = @import("../../../float.zig");
const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Solves a system of linear equations, whose coefficients are in a triangular
/// matrix, defined as:
///
/// ```zig
///     A * x = b,
/// ```
///
/// or
///
/// ```zig
///     conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^T * x = b,
/// ```
///
/// or
///
/// ```zig
///     A^H * x = b,
/// ```
///
/// where `b` and `x` are `n`-element vectors, `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.trsv(layout: matrix.Layout, uplo: matrix.Uplo, transa: linalg.blas.Transpose, diag: matrix.Diag, n: usize, a: [*]const A, lda: usize, x: [*]X, incx: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `uplo` (`matrix.Uplo`): Specifies whether the matrix `A` is an upper or
///   lower triangular matrix.
/// * `transa` (`linalg.blas.Transpose`): Specifies the system of equations to
///   be solved:
///   * `no_transpose`: `A * x = b`
///   * `transpose`: `Aᵀ * x = b`
///   * `conj_no_transpose`: `conj(A) * x = b`
///   * `conj_transpose`: `Aᴴ * x = b`
/// * `diag` (`matrix.Diag`): Specifies whether the matrix `A` is unit
///   triangular.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * n`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, n)`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. On entry, must contain the right hand side
///   vector `b`, on return, contains the solution vector `x`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, n)`, or
///   if `incx` is 0.
pub fn trsv(
    layout: matrix.Layout,
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    diag: matrix.Diag,
    n: usize,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
) !void {
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.trsv: a must be a many-item pointer to numerics, and x must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);

    if (lda < int.max(1, n) or incx == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const eff_transa = if (layout == .col_major) transa else transa.invert();
    const noconj = transa == .no_trans or transa == .trans;

    // Quick return if possible.
    if (n == 0)
        return;

    if (noconj)
        k_trsv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, true)
    else
        k_trsv(eff_uplo, eff_transa, diag, n, a, lda, x, incx, false);
}

pub fn k_trsv(
    uplo: matrix.Uplo,
    transa: linalg.blas.Transpose,
    diag: matrix.Diag,
    n: usize,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    comptime noconj: bool,
) void {
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));

    // Set up the start point in x.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  x = A⁻¹ * x  or  x = conj(A)⁻¹ * x.
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(X, A, X)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        if (uplo == .lower) {
            var tile_j: usize = 0;
            while (tile_j < n) : (tile_j += tile_size) {
                const c_end = int.min(tile_j + tile_size, n);

                var j: usize = tile_j;
                while (j < c_end) : (j += 1) {
                    const jx = kx + numeric.cast(isize, j) * incx;

                    if (diag == .non_unit) {
                        numeric.divInto(
                            &x[numeric.cast(usize, jx)],
                            x[numeric.cast(usize, jx)],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                        );
                    }

                    const temp = numeric.neg(x[numeric.cast(usize, jx)]);

                    if (incx == 1) {
                        var i: usize = j + 1;
                        while (i < (c_end / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                numeric.fmaInto(
                                    &x[i + u],
                                    temp,
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[i + u],
                                );
                            }
                        }

                        while (i < c_end) : (i += 1) {
                            numeric.fmaInto(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                            );
                        }
                    } else {
                        var ix = kx + numeric.cast(isize, j + 1) * incx;
                        var i: usize = j + 1;
                        while (i < c_end) : (i += 1) {
                            numeric.fmaInto(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix += incx;
                        }
                    }
                }

                if (c_end < n) {
                    var tile_i: usize = c_end;
                    while (tile_i < n) : (tile_i += tile_size) {
                        const r_end = int.min(tile_i + tile_size, n);
                        const r_len = r_end - tile_i;
                        var local_x: [tile_size]X = undefined;

                        const px = if (incx == 1)
                            x + tile_i
                        else blk: {
                            linalg.blas.copy(
                                r_len,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                                @as([*]X, &local_x),
                                1,
                            ) catch unreachable;

                            break :blk @as([*]X, &local_x);
                        };

                        j = tile_j;
                        while (j < c_end) : (j += 1) {
                            const temp = numeric.neg(x[numeric.cast(usize, kx + numeric.cast(isize, j) * incx)]);

                            var i: usize = 0;
                            while (i < (r_len / unroll) * unroll) : (i += unroll) {
                                inline for (0..unroll) |u| {
                                    numeric.fmaInto(
                                        &px[i + u],
                                        temp,
                                        if (comptime noconj)
                                            a[tile_i + i + u + j * lda]
                                        else
                                            numeric.conj(a[tile_i + i + u + j * lda]),
                                        px[i + u],
                                    );
                                }
                            }

                            while (i < r_len) : (i += 1) {
                                numeric.fmaInto(
                                    &px[i],
                                    temp,
                                    if (comptime noconj)
                                        a[tile_i + i + j * lda]
                                    else
                                        numeric.conj(a[tile_i + i + j * lda]),
                                    px[i],
                                );
                            }
                        }

                        if (incx != 1) {
                            linalg.blas.copy(
                                r_len,
                                @as([*]const X, &local_x),
                                1,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                            ) catch unreachable;
                        }
                    }
                }
            }
        } else {
            var ct: usize = (n + tile_size - 1) / tile_size;
            while (ct > 0) {
                ct -= 1;
                const tile_j = ct * tile_size;
                const c_end = int.min(tile_j + tile_size, n);

                var j: usize = c_end;
                while (j > tile_j) {
                    j -= 1;
                    const jx = kx + numeric.cast(isize, j) * incx;

                    if (diag == .non_unit) {
                        numeric.divInto(
                            &x[numeric.cast(usize, jx)],
                            x[numeric.cast(usize, jx)],
                            if (comptime noconj)
                                a[j + j * lda]
                            else
                                numeric.conj(a[j + j * lda]),
                        );
                    }

                    const temp = numeric.neg(x[numeric.cast(usize, jx)]);

                    if (incx == 1) {
                        var i: usize = tile_j;
                        while (i < (j / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                numeric.fmaInto(
                                    &x[i + u],
                                    temp,
                                    if (comptime noconj)
                                        a[i + u + j * lda]
                                    else
                                        numeric.conj(a[i + u + j * lda]),
                                    x[i + u],
                                );
                            }
                        }

                        while (i < j) : (i += 1) {
                            numeric.fmaInto(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[i],
                            );
                        }
                    } else {
                        var ix = kx + numeric.cast(isize, tile_j) * incx;
                        var i: usize = tile_j;
                        while (i < j) : (i += 1) {
                            numeric.fmaInto(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    a[i + j * lda]
                                else
                                    numeric.conj(a[i + j * lda]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix += incx;
                        }
                    }
                }

                if (tile_j > 0) {
                    var tile_i: usize = 0;
                    while (tile_i < tile_j) : (tile_i += tile_size) {
                        const r_end = int.min(tile_i + tile_size, tile_j);
                        const r_len = r_end - tile_i;
                        var local_x: [tile_size]X = undefined;

                        const px = if (incx == 1)
                            x + tile_i
                        else blk: {
                            linalg.blas.copy(
                                r_len,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                                @as([*]X, &local_x),
                                1,
                            ) catch unreachable;

                            break :blk @as([*]X, &local_x);
                        };

                        j = tile_j;
                        while (j < c_end) : (j += 1) {
                            const temp = numeric.neg(x[numeric.cast(usize, kx + numeric.cast(isize, j) * incx)]);

                            var i: usize = 0;
                            while (i < (r_len / unroll) * unroll) : (i += unroll) {
                                inline for (0..unroll) |u| {
                                    numeric.fmaInto(
                                        &px[i + u],
                                        temp,
                                        if (comptime noconj)
                                            a[tile_i + i + u + j * lda]
                                        else
                                            numeric.conj(a[tile_i + i + u + j * lda]),
                                        px[i + u],
                                    );
                                }
                            }

                            while (i < r_len) : (i += 1) {
                                numeric.fmaInto(
                                    &px[i],
                                    temp,
                                    if (comptime noconj)
                                        a[tile_i + i + j * lda]
                                    else
                                        numeric.conj(a[tile_i + i + j * lda]),
                                    px[i],
                                );
                            }
                        }

                        if (incx != 1) {
                            linalg.blas.copy(
                                r_len,
                                @as([*]const X, &local_x),
                                1,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_i) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_i + r_len)) * incx),
                                incx,
                            ) catch unreachable;
                        }
                    }
                }
            }
        }
    } else {
        // Form  x = A⁻ᵀ * x  or  x = A⁻ᴴ * x.
        const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, X))) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(meta.Accumulator(numeric.Mul(A, X))) + @sizeOf(X) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        if (uplo == .upper) {
            var tile_i: usize = 0;
            while (tile_i < n) : (tile_i += tile_size) {
                const i_end = int.min(tile_i + tile_size, n);
                const i_len = i_end - tile_i;

                if (tile_i > 0) {
                    var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                    @memset(local_sums[0..i_len], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                    var tile_j: usize = 0;
                    while (tile_j < tile_i) : (tile_j += tile_size) {
                        const j_end = int.min(tile_j + tile_size, tile_i);
                        const j_len = j_end - tile_j;
                        var local_x: [tile_size]X = undefined;

                        const px = if (incx == 1)
                            x + tile_j
                        else blk: {
                            linalg.blas.copy(
                                j_len,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_j) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_j + j_len)) * incx),
                                incx,
                                @as([*]X, &local_x),
                                1,
                            ) catch unreachable;

                            break :blk @as([*]const X, &local_x);
                        };

                        var i: usize = tile_i;
                        while (i < i_end) : (i += 1) {
                            var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                            var j: usize = 0;
                            while (j < (j_len / unroll) * unroll) : (j += unroll) {
                                inline for (0..unroll) |u| {
                                    numeric.fmaInto(
                                        &sums[u],
                                        if (comptime noconj)
                                            a[tile_j + j + u + i * lda]
                                        else
                                            numeric.conj(a[tile_j + j + u + i * lda]),
                                        px[j + u],
                                        sums[u],
                                    );
                                }
                            }

                            var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));
                            inline for (0..unroll) |u| {
                                numeric.addInto(&temp, temp, sums[u]);
                            }

                            while (j < j_len) : (j += 1) {
                                numeric.fmaInto(
                                    &temp,
                                    if (comptime noconj)
                                        a[tile_j + j + i * lda]
                                    else
                                        numeric.conj(a[tile_j + j + i * lda]),
                                    px[j],
                                    temp,
                                );
                            }

                            numeric.addInto(&local_sums[i - tile_i], local_sums[i - tile_i], temp);
                        }
                    }

                    var i: usize = tile_i;
                    while (i < i_end) : (i += 1) {
                        const ix = kx + numeric.cast(isize, i) * incx;
                        var val = x[numeric.cast(usize, ix)];

                        numeric.subInto(
                            &val,
                            val,
                            numeric.cast(X, local_sums[i - tile_i]),
                        );

                        numeric.set(&x[numeric.cast(usize, ix)], val);
                    }
                }

                var i: usize = tile_i;
                while (i < i_end) : (i += 1) {
                    const ix = kx + numeric.cast(isize, i) * incx;
                    var sum = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                    var j: usize = tile_i;
                    while (j < i) : (j += 1) {
                        const jx = kx + numeric.cast(isize, j) * incx;
                        numeric.fmaInto(
                            &sum,
                            if (comptime noconj)
                                a[j + i * lda]
                            else
                                numeric.conj(a[j + i * lda]),
                            x[numeric.cast(usize, jx)],
                            sum,
                        );
                    }

                    var temp = x[numeric.cast(usize, ix)];
                    numeric.subInto(
                        &temp,
                        temp,
                        numeric.cast(X, sum),
                    );

                    if (diag == .non_unit) {
                        numeric.divInto(
                            &temp,
                            temp,
                            if (comptime noconj)
                                a[i + i * lda]
                            else
                                numeric.conj(a[i + i * lda]),
                        );
                    }

                    numeric.set(&x[numeric.cast(usize, ix)], temp);
                }
            }
        } else {
            var rt: usize = (n + tile_size - 1) / tile_size;
            while (rt > 0) {
                rt -= 1;
                const tile_i = rt * tile_size;
                const i_end = int.min(tile_i + tile_size, n);
                const i_len = i_end - tile_i;

                if (i_end < n) {
                    var local_sums: [tile_size]meta.Accumulator(numeric.Mul(A, X)) = undefined;
                    @memset(local_sums[0..i_len], numeric.zero(meta.Accumulator(numeric.Mul(A, X))));

                    var tile_j: usize = i_end;
                    while (tile_j < n) : (tile_j += tile_size) {
                        const j_end = int.min(tile_j + tile_size, n);
                        const j_len = j_end - tile_j;
                        var local_x: [tile_size]X = undefined;

                        const px = if (incx == 1)
                            x + tile_j
                        else blk: {
                            linalg.blas.copy(
                                j_len,
                                x + numeric.cast(usize, if (incx > 0)
                                    numeric.cast(isize, tile_j) * incx
                                else
                                    (-numeric.cast(isize, n) + numeric.cast(isize, tile_j + j_len)) * incx),
                                incx,
                                @as([*]X, &local_x),
                                1,
                            ) catch unreachable;

                            break :blk @as([*]const X, &local_x);
                        };

                        var i: usize = tile_i;
                        while (i < i_end) : (i += 1) {
                            var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                            var j: usize = 0;
                            while (j < (j_len / unroll) * unroll) : (j += unroll) {
                                inline for (0..unroll) |u| {
                                    numeric.fmaInto(
                                        &sums[u],
                                        if (comptime noconj)
                                            a[tile_j + j + u + i * lda]
                                        else
                                            numeric.conj(a[tile_j + j + u + i * lda]),
                                        px[j + u],
                                        sums[u],
                                    );
                                }
                            }

                            var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));
                            inline for (0..unroll) |u| {
                                numeric.addInto(
                                    &temp,
                                    temp,
                                    sums[u],
                                );
                            }

                            while (j < j_len) : (j += 1) {
                                numeric.fmaInto(
                                    &temp,
                                    if (comptime noconj)
                                        a[tile_j + j + i * lda]
                                    else
                                        numeric.conj(a[tile_j + j + i * lda]),
                                    px[j],
                                    temp,
                                );
                            }

                            numeric.addInto(
                                &local_sums[i - tile_i],
                                local_sums[i - tile_i],
                                temp,
                            );
                        }
                    }

                    var i: usize = tile_i;
                    while (i < i_end) : (i += 1) {
                        const ix = kx + numeric.cast(isize, i) * incx;
                        var val = x[numeric.cast(usize, ix)];

                        numeric.subInto(
                            &val,
                            val,
                            numeric.cast(X, local_sums[i - tile_i]),
                        );

                        numeric.set(&x[numeric.cast(usize, ix)], val);
                    }
                }

                var i: usize = i_end;
                while (i > tile_i) {
                    i -= 1;
                    const ix = kx + numeric.cast(isize, i) * incx;
                    var sum = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                    var j: usize = i + 1;
                    while (j < i_end) : (j += 1) {
                        const jx = kx + numeric.cast(isize, j) * incx;
                        numeric.fmaInto(
                            &sum,
                            if (comptime noconj)
                                a[j + i * lda]
                            else
                                numeric.conj(a[j + i * lda]),
                            x[numeric.cast(usize, jx)],
                            sum,
                        );
                    }

                    var temp = x[numeric.cast(usize, ix)];
                    numeric.subInto(&temp, temp, numeric.cast(X, sum));

                    if (diag == .non_unit) {
                        numeric.divInto(
                            &temp,
                            temp,
                            if (comptime noconj)
                                a[i + i * lda]
                            else
                                numeric.conj(a[i + i * lda]),
                        );
                    }

                    numeric.set(&x[numeric.cast(usize, ix)], temp);
                }
            }
        }
    }
}
