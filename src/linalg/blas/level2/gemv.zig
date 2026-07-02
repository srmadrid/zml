const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");
const matrix = @import("../../../matrix.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

/// Computes a matrix-vector product with a general matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᵀ * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᴴ * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are vectors, and `A` is an
/// `m`-by-`n` matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.gemv(layout: matrix.Layout, transa: linalg.blas.Transpose, m: usize, n: usize, alpha: Al, a: [*]const A, lda: usize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `transa` (`linalg.blas.Transpose`): Specifies the operation to be
///   performed on `A`:
///   * `no_transpose`: `y = alpha * A * x + beta * y`
///   * `transpose`: `y = alpha * Aᵀ * x + beta * y`
///   * `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
///   * `conj_transpose`: `y = alpha * Aᴴ * x + beta * y`
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * k`, where `k` is
///   `n` when `layout` is `col_major`, or `m` when `layout` is `row_major`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (m - 1) * abs(incx)` otherwise.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `y` need not be set on input.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (m - 1) * abs(incy)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (n - 1) * abs(incy)` otherwise. On return,
///   contains the result of the operation.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `opts`: Optional parameters:
///   * `num_threads` (`usize = 0`): Number of threads to spawn:
///     * `0`: automatic. The thread count is derived from `m * n` and
///       `parallel_threshold`:
///       ```zig
///       threads = max(1, min(std.Thread.getCpuCount(), options.max_threads, (m * n) / parallel_threshold))
///       ```
///     * 1: force serial execution. parallel_threshold is ignored.
///     * N >= 2: use exactly N threads, clamped by
///       std.Thread.getCpuCount() and options.max_threads as a hard safety
///       ceiling. parallel_threshold is ignored.
///   * parallel_threshold (usize = 2_097_152 / @sizeOf(meta.Child(Y))):
///     Minimum number of matrix elements (`m * n`) required to trigger
///     multithreaded execution.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, m)` or
///   `max(1, n)`, or if `incx` or `incy` is 0.
pub fn gemv(
    layout: matrix.Layout,
    transa: linalg.blas.Transpose,
    m: usize,
    n: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 2_097_152 / @sizeOf(meta.Child(@TypeOf(y))),
    },
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
        @compileError("zsl.linalg.blas.gemv: alpha and beta must be numerics, a and x must be many-item pointers to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);
    Y = meta.Child(Y);

    if (lda < int.max(1, if (layout == .col_major) m else n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_transa = if (layout == .col_major) transa else transa.invert();
    const eff_m = if (layout == .col_major) m else n;
    const eff_n = if (layout == .col_major) n else m;
    const noconj = transa == .no_trans or transa == .trans;
    const no_trans = eff_transa == .no_trans or eff_transa == .conj_no_trans;
    const leny = if (no_trans) eff_m else eff_n;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return;

    if (numeric.eq(alpha, 0)) {
        if (numeric.ne(beta, 1))
            @import("../level1/scal.zig").k_scal(leny, beta, y, incy);

        return;
    }

    if (opts.num_threads == 1)
        return if (noconj)
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, false);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, (eff_m * eff_n) / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);
    num_threads = int.min(num_threads, leny);

    if (num_threads <= 1)
        return if (noconj)
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, false);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);
    num_threads = int.min(num_threads, leny);

    if (num_threads <= 1)
        return if (noconj)
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, false);

    const Worker = struct {
        fn execute(
            worker_transa: linalg.blas.Transpose,
            worker_m: usize,
            worker_n: usize,
            worker_alpha: Al,
            worker_a: [*]const A,
            worker_lda: usize,
            worker_x: [*]const X,
            worker_incx: isize,
            worker_beta: Be,
            worker_y: [*]Y,
            worker_incy: isize,
            worker_noconj: bool,
        ) void {
            return if (worker_noconj)
                k_gemv(worker_transa, worker_m, worker_n, worker_alpha, worker_a, worker_lda, worker_x, worker_incx, worker_beta, worker_y, worker_incy, true)
            else
                k_gemv(worker_transa, worker_m, worker_n, worker_alpha, worker_a, worker_lda, worker_x, worker_incx, worker_beta, worker_y, worker_incy, false);
        }
    };

    var threads: [options.max_threads]std.Thread = undefined;

    const chunk_size = int.div(leny, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) leny else chunk_start + chunk_size;
        const chunk_len = chunk_end - chunk_start;

        if (std.Thread.spawn(.{}, Worker.execute, .{
            eff_transa,
            if (no_trans) chunk_len else eff_m,
            if (no_trans) eff_n else chunk_len,
            alpha,
            a + if (no_trans) chunk_start else (chunk_start * lda),
            lda,
            x,
            incx,
            beta,
            y + numeric.cast(usize, if (incy > 0)
                numeric.cast(isize, chunk_start) * incy
            else
                (-numeric.cast(isize, leny) + numeric.cast(isize, chunk_end)) * incy),
            incy,
            noconj,
        })) |th| {
            threads[i] = th;
            spawned_count += 1;
        } else |err| {
            spawn_err = err;
            break;
        }
    }

    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
    }

    if (spawn_err) |err|
        return err;
}

pub fn k_gemv(transa: linalg.blas.Transpose, m: usize, n: usize, alpha: anytype, a: anytype, lda: usize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Quick return if possible.
    if (m == 0 or n == 0)
        return;

    // Set lenx and leny, the lengths of the vectors x and y.
    const lenx: usize = if (transa == .no_trans or transa == .conj_no_trans) n else m;
    const leny: usize = if (transa == .no_trans or transa == .conj_no_trans) m else n;

    // First form  y = beta * y.
    if (numeric.ne(beta, 1))
        @import("../level1/scal.zig").k_scal(leny, beta, y, incy);

    if (numeric.eq(alpha, 0))
        return;

    // Set up the start points in x and y.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, lenx) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, leny) + 1) * incy else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  y = alpha * A * x + y  or  y = alpha * conj(A) * x + y.
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(Y) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < m) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, m - tile_i);
            var local_y: [tile_size]Y = undefined;

            const py = if (incy == 1)
                y + tile_i
            else blk: {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, m) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                    @as([*]Y, &local_y),
                    1,
                );

                break :blk @as([*]Y, &local_y);
            };

            var jx: isize = kx;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp = alpha * x[jx]
                const temp = numeric.mul(
                    alpha,
                    x[numeric.cast(usize, jx)],
                );

                var i: usize = 0;
                while (i < (b_len / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // py[i + u] += temp * a[tile_i + i + u + j * lda]
                        numeric.fmaInto(
                            &py[i + u],
                            temp,
                            if (comptime noconj)
                                a[tile_i + i + u + j * lda]
                            else
                                numeric.conj(a[tile_i + i + u + j * lda]),
                            py[i + u],
                        );
                    }
                }

                while (i < b_len) : (i += 1) {
                    // py[i] += temp * a[tile_i + i + j * lda]
                    numeric.fmaInto(
                        &py[i],
                        temp,
                        if (comptime noconj)
                            a[tile_i + i + j * lda]
                        else
                            numeric.conj(a[tile_i + i + j * lda]),
                        py[i],
                    );
                }

                jx += incx;
            }

            if (incy != 1) {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    @as([*]Y, &local_y),
                    1,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, m) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                );
            }
        }
    } else {
        // Form  y = alpha * Aᵀ * x + y  or  y = alpha * Aᴴ * x + y.
        const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, X))) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size -|= tile_size % unroll;

        var tile_i: usize = 0;
        while (tile_i < m) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, m - tile_i);
            var local_x: [tile_size]X = undefined;

            const px = if (incx == 1)
                x + tile_i
            else blk: {
                @import("../level1/copy.zig").k_copy(
                    b_len,
                    x + numeric.cast(usize, if (incx > 0)
                        numeric.cast(isize, tile_i) * incx
                    else
                        (-numeric.cast(isize, m) + numeric.cast(isize, tile_i + b_len)) * incx),
                    incx,
                    @as([*]X, &local_x),
                    1,
                );

                break :blk @as([*]const X, &local_x);
            };

            var jy: isize = ky;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                var i: usize = 0;
                while (i < (b_len / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // sums[u] += a[tile_i + i + u + j * lda] * px[i + u]
                        numeric.fmaInto(
                            &sums[u],
                            if (comptime noconj)
                                a[tile_i + i + u + j * lda]
                            else
                                numeric.conj(a[tile_i + i + u + j * lda]),
                            px[i + u],
                            sums[u],
                        );
                    }
                }

                inline for (0..unroll) |u| {
                    numeric.addInto(&temp, temp, sums[u]);
                }

                while (i < b_len) : (i += 1) {
                    // temp += a[tile_i + i + j * lda] * x[i]
                    numeric.fmaInto(
                        &temp,
                        if (comptime noconj)
                            a[tile_i + i + j * lda]
                        else
                            numeric.conj(a[tile_i + i + j * lda]),
                        px[i],
                        temp,
                    );
                }

                // y[jy] += alpha * temp
                numeric.fmaInto(
                    &y[numeric.cast(usize, jy)],
                    alpha,
                    temp,
                    y[numeric.cast(usize, jy)],
                );

                jy += incy;
            }
        }
    }

    return;
}
