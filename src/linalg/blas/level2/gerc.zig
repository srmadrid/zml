const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");
const matrix = @import("../../../matrix.zig");

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

/// Performs a rank-1 update (conjugated) of a general matrix defined as:
///
/// ```zig
/// A = alpha * x * yᴴ + A,
/// ```
///
/// where `alpha` is a numeric, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.gerc(layout: matrix.Layout, m: usize, n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, a: [*]A, lda: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (m - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * k`, where
///   `k` is `n` when `layout` is `col_major`, or `m` when `layout` is
///   `row_major`. On return, contains the result of the operation.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
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
///   * parallel_threshold (usize = 2_097_152 / @sizeOf(meta.Child(A))):
///     Minimum number of matrix elements (`m * n`) required to trigger
///     multithreaded execution.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0, or if `lda`
///   is less than `max(1, m)` or `max(1, n)`.
pub fn gerc(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: usize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 2_097_152 / @sizeOf(meta.Child(@TypeOf(a))),
    },
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or !meta.isNumeric(meta.Child(Y)) or
        !meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.blas.gerc: alpha must be a numeric, x and y must be many-item pointers to numerics, and a must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\ta: " ++ @typeName(A) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);
    A = meta.Child(A);

    if (lda < int.max(1, if (layout == .col_major) m else n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_m = if (layout == .col_major) m else n;
    const eff_n = if (layout == .col_major) n else m;

    // Quick return if possible.
    if (m == 0 or n == 0 or numeric.eq(alpha, 0))
        return;

    if (opts.num_threads == 1)
        return if (layout == .col_major)
            k_gerc(m, n, alpha, x, incx, y, incy, a, lda, true)
        else
            k_gerc(n, m, alpha, y, incy, x, incx, a, lda, false);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, (eff_m * eff_n) / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);
    num_threads = int.min(num_threads, eff_n);

    if (num_threads <= 1)
        return if (layout == .col_major)
            k_gerc(m, n, alpha, x, incx, y, incy, a, lda, true)
        else
            k_gerc(n, m, alpha, y, incy, x, incx, a, lda, false);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);
    num_threads = int.min(num_threads, eff_n);

    if (num_threads <= 1)
        return if (layout == .col_major)
            k_gerc(m, n, alpha, x, incx, y, incy, a, lda, true)
        else
            k_gerc(n, m, alpha, y, incy, x, incx, a, lda, false);

    const Worker = struct {
        fn execute(worker_layout: matrix.Layout, worker_m: usize, worker_n: usize, worker_alpha: Al, worker_x: [*]const X, worker_incx: isize, worker_y: [*]const Y, worker_incy: isize, worker_a: [*]A, worker_lda: usize) void {
            return if (worker_layout == .col_major)
                k_gerc(worker_m, worker_n, worker_alpha, worker_x, worker_incx, worker_y, worker_incy, worker_a, worker_lda, true)
            else
                k_gerc(worker_n, worker_m, worker_alpha, worker_y, worker_incy, worker_x, worker_incx, worker_a, worker_lda, false);
        }
    };

    var threads: [options.max_threads]std.Thread = undefined;

    const chunk_size = int.div(eff_n, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) eff_n else chunk_start + chunk_size;
        const chunk_len = chunk_end - chunk_start;

        if (std.Thread.spawn(.{}, Worker.execute, .{
            layout,
            if (layout == .col_major) m else chunk_len,
            if (layout == .col_major) chunk_len else n,
            alpha,
            if (layout == .col_major) x else x + numeric.cast(usize, if (incx > 0)
                numeric.cast(isize, chunk_start) * incx
            else
                (-numeric.cast(isize, m) + numeric.cast(isize, chunk_end)) * incx),
            incx,
            if (layout == .col_major) y + numeric.cast(usize, if (incy > 0)
                numeric.cast(isize, chunk_start) * incy
            else
                (-numeric.cast(isize, n) + numeric.cast(isize, chunk_end)) * incy) else y,
            incy,
            a + chunk_start * lda,
            lda,
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

fn k_gerc(m: usize, n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize, a: anytype, lda: usize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Quick return if possible.
    if (m == 0 or n == 0 or numeric.eq(alpha, 0))
        return;

    // Form  A = alpha * x * yᴴ + A
    const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(if (comptime noconj) X else numeric.Conj(X), numeric.Mul(Al, if (comptime noconj) numeric.Conj(Y) else Y), A)) orelse 2);
    comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
    tile_size = comptime int.max(1, tile_size -| tile_size % unroll);

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

        var jy: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;
        var j: usize = 0;
        while (j < n) : (j += 1) {
            if (numeric.ne(y[numeric.cast(usize, jy)], 0)) {
                // temp = alpha * conj(y[jy])
                const temp = numeric.mul(
                    alpha,
                    if (comptime noconj)
                        numeric.conj(y[numeric.cast(usize, jy)])
                    else
                        y[numeric.cast(usize, jy)],
                );

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

            jy += incy;
        }
    }

    return;
}
