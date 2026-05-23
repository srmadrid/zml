const std = @import("std");
const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");

/// Size of tiles to use when multithreading.
const tile_size = 128;

/// Performs a rank-2 update of a Hermitian matrix defined as:
///
/// ```zig
/// A = alpha * x * yᴴ + conjg(alpha) * y * xᴴ + A,
/// ```
///
/// where `alpha` is a numeric, `x` and `y` are `n`-element vectors, and `A` is
/// an `n`-by-`n` Hermitian matrix.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.her2(layout: Layout, uplo: Uplo, n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, a: [*]A, lda: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the upper or lower triangular part of
///   the Hermitian packed matrix `A` is used.
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
/// * `opts`: Optional parameters:
///   * `num_threads` (`usize = 0`): Number of threads to spawn:
///     * `0`: automatic. The thread count is derived from `m * n` and
///       `parallel_threshold`:
///       ```zig
///       threads = max(1, min(std.Thread.getCpuCount(), options.max_threads, (n * n) / parallel_threshold))
///       ```
///     * 1: force serial execution. parallel_threshold is ignored.
///     * N >= 2: use exactly N threads, clamped by
///       std.Thread.getCpuCount() and options.max_threads as a hard safety
///       ceiling. parallel_threshold is ignored.
///   * parallel_threshold (usize = 4_194_304 / @sizeOf(meta.Child(A))):
///     Minimum number of matrix elements (`n * n`) required to trigger
///     multithreaded execution.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0, or if `lda`
///   is less than `max(1, n)`.
pub fn her2(
    layout: Layout,
    uplo: Uplo,
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
        parallel_threshold: usize = 4_194_304 / @sizeOf(meta.Child(@TypeOf(a))),
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
        @compileError("zsl.linalg.blas.her2: alpha must be a real numeric, x and y must be many-item pointers to numerics, and a must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\ta: " ++ @typeName(A) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);
    A = meta.Child(A);

    if (lda < int.max(1, n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == X and Al == Y and Al == A) {
        switch (comptime meta.numericType(Al)) {
            .complex => {
                if (comptime meta.Scalar(Al) == f32)
                    return linalg.cblas.cher(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), &alpha, x, incx, y, incy, a, numeric.cast(isize, lda))
                else if (comptime meta.Scalar(Al) == f64)
                    return linalg.cblas.zher(layout.toInt(c_int), uplo.toInt(c_int), numeric.cast(isize, n), &alpha, x, incx, y, incy, a, numeric.cast(isize, lda));
            },
            else => {},
        }
    }

    const eff_uplo = if (layout == .col_major) uplo else uplo.invert();
    const noconj = layout == .col_major;

    // Quick return if possible.
    if (n == 0 or numeric.eq(alpha, 0))
        return;

    if (opts.num_threads == 1)
        return if (noconj)
            k_her2(eff_uplo, n, alpha, x, incx, y, incy, a, lda, true)
        else
            k_her2(eff_uplo, n, alpha, x, incx, y, incy, a, lda, false);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, (n * n) / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return if (noconj)
            k_her2(eff_uplo, n, alpha, x, incx, y, incy, a, lda, true)
        else
            k_her2(eff_uplo, n, alpha, x, incx, y, incy, a, lda, false);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);
    num_threads = int.min(num_threads, n);

    if (num_threads <= 1)
        return if (noconj)
            k_her2(eff_uplo, n, alpha, x, incx, y, incy, a, lda, true)
        else
            k_her2(eff_uplo, n, alpha, x, incx, y, incy, a, lda, false);

    const k = (n + tile_size - 1) / tile_size;
    const num_tiles = k * (k + 1) / 2;

    var atomic_counter = std.atomic.Value(usize).init(0);
    var threads: [options.max_threads]std.Thread = undefined;

    const Worker = struct {
        fn execute(
            worker_uplo: Uplo,
            worker_n: usize,
            worker_alpha: Al,
            worker_a: [*]A,
            worker_lda: usize,
            worker_x: [*]const X,
            worker_incx: isize,
            worker_y: [*]const Y,
            worker_incy: isize,
            comptime worker_noconj: bool,
            counter: *std.atomic.Value(usize),
            comptime worker_tile_size: comptime_int,
            worker_num_tiles: usize,
        ) void {
            while (true) {
                const idx = counter.fetchAdd(1, .monotonic);

                if (idx >= worker_num_tiles) // When all tiles have been assigned, break.
                    break;

                // Map 1D atomic index to 2D upper triangular coordinates (tile_i, tile_j) using triangular numbers.
                var tile_j = numeric.cast(usize, (float.sqrt(1.0 + 8.0 * numeric.cast(f64, idx)) - 1.0) / 2.0);

                while (tile_j * (tile_j + 1) / 2 > idx)
                    tile_j -= 1;

                while ((tile_j + 1) * (tile_j + 2) / 2 <= idx)
                    tile_j += 1;

                const tile_i = idx - tile_j * (tile_j + 1) / 2;

                const phys_r = if (worker_uplo == .upper) tile_i else tile_j;
                const phys_c = if (worker_uplo == .upper) tile_j else tile_i;

                const r_start = phys_r * worker_tile_size;
                const c_start = phys_c * worker_tile_size;
                const r_len = int.min(worker_tile_size, worker_n - r_start);
                const c_len = int.min(worker_tile_size, worker_n - c_start);

                if (tile_i == tile_j) {
                    // Diagonal tile
                    var local_x: [worker_tile_size]X = undefined;
                    var local_y: [worker_tile_size]Y = undefined;

                    @import("copy.zig").k_copy(
                        r_len,
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, r_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx),
                        worker_incx,
                        @as([*]X, &local_x),
                        1,
                    );

                    @import("copy.zig").k_copy(
                        c_len,
                        worker_y + numeric.cast(usize, if (worker_incy > 0)
                            numeric.cast(isize, r_start) * worker_incy
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incy),
                        worker_incy,
                        @as([*]Y, &local_y),
                        1,
                    );

                    k_her2(
                        worker_uplo,
                        r_len,
                        worker_alpha,
                        @as([*]const X, &local_x),
                        1,
                        @as([*]const Y, &local_y),
                        1,
                        worker_a + r_start + c_start * worker_lda,
                        worker_lda,
                        worker_noconj,
                    );
                } else {
                    const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(if (comptime worker_noconj) numeric.Conj(X) else X, numeric.Mul(Al, if (comptime worker_noconj) X else numeric.Conj(X)), A)) orelse 2);

                    // Off-diagonal tile
                    var local_x_r: [worker_tile_size]X = undefined;
                    var local_x_c: [worker_tile_size]X = undefined;

                    @import("copy.zig").k_copy(
                        r_len,
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, r_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, r_start + r_len)) * worker_incx),
                        worker_incx,
                        @as([*]X, &local_x_r),
                        1,
                    );

                    @import("copy.zig").k_copy(
                        c_len,
                        worker_x + numeric.cast(usize, if (worker_incx > 0)
                            numeric.cast(isize, c_start) * worker_incx
                        else
                            (-numeric.cast(isize, worker_n) + numeric.cast(isize, c_start + c_len)) * worker_incx),
                        worker_incx,
                        @as([*]X, &local_x_c),
                        1,
                    );

                    var j: usize = 0;
                    while (j < c_len) : (j += 1) {
                        // temp = worker_alpha * conj(local_x_c[j])
                        const temp = numeric.mul(
                            worker_alpha,
                            if (worker_noconj)
                                numeric.conj(local_x_c[j])
                            else
                                local_x_c[j],
                        );

                        var i: usize = 0;
                        while (i < (r_len / unroll) * unroll) : (i += unroll) {
                            inline for (0..unroll) |u| {
                                // worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda] += temp * local_x_r[i + u]
                                numeric.fma_(
                                    &worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda],
                                    temp,
                                    if (comptime worker_noconj)
                                        local_x_r[i + u]
                                    else
                                        numeric.conj(local_x_r[i + u]),
                                    worker_a[r_start + c_start * worker_lda + i + u + j * worker_lda],
                                );
                            }
                        }

                        while (i < r_len) : (i += 1) {
                            // worker_a[r_start + c_start * worker_lda + i + j * worker_lda] += temp * local_x_r[i]
                            numeric.fma_(
                                &worker_a[r_start + c_start * worker_lda + i + j * worker_lda],
                                temp,
                                if (comptime worker_noconj)
                                    local_x_r[i]
                                else
                                    numeric.conj(local_x_r[i]),
                                worker_a[r_start + c_start * worker_lda + i + j * worker_lda],
                            );
                        }
                    }
                }
            }
        }
    };

    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;

    while (i < num_threads) : (i += 1) {
        if (if (noconj)
            std.Thread.spawn(.{}, Worker.execute, .{
                eff_uplo,
                n,
                alpha,
                a,
                lda,
                x,
                incx,
                true,
                &atomic_counter,
                tile_size,
                num_tiles,
            })
        else
            std.Thread.spawn(.{}, Worker.execute, .{
                eff_uplo,
                n,
                alpha,
                a,
                lda,
                x,
                incx,
                false,
                &atomic_counter,
                tile_size,
                num_tiles,
            })) |th|
        {
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

fn k_her2(
    uplo: Uplo,
    n: isize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: isize,
    noconj: bool,
    ctx: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = types.Child(@TypeOf(x));
    const C2: type = types.Coerce(Al, X);
    const Y: type = types.Child(@TypeOf(y));
    const C1: type = types.Coerce(Al, Y);
    const A: type = types.Child(@TypeOf(a));
    const CC: type = types.Coerce(Al, types.Coerce(X, types.Coerce(Y, A)));

    if (n < 0 or lda < int.max(1, n) or incx == 0 or incy == 0)
        return blas.Error.InvalidArgument;

    // Quick return if possible.
    if (n == 0 or ops.eq(alpha, 0, ctx) catch unreachable)
        return;

    const kx: isize = if (incx < 0) (-n + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-n + 1) * incy else 0;

    if (comptime !types.isArbitraryPrecision(CC)) {
        if (uplo == .upper) {
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[j])
                                ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[j])
                                x[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += x[i] * temp1 + y[i] * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, i)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, i)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[j] * temp1 + y[j] * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, j)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, j)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[jy])
                                ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[jx])
                                x[scast(u32, jx)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var ix: isize = kx;
                            var iy: isize = ky;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += x[ix] * temp1 + y[iy] * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, ix)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, iy)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                                iy += incy;
                            }

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[jx] * temp1 + y[jy] * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, jx)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, jy)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        jy += incy;
                    }
                }
            } else {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[j]
                                y[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[j]))
                                ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += conj(x[i]) * temp1 + conj(y[i]) * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, i)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, i)], ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[j]) * temp1 + conj(y[j]) * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[jy]
                                y[scast(u32, jy)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[jx]))
                                ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            var ix: isize = kx;
                            var iy: isize = ky;
                            var i: isize = 0;
                            while (i < j) : (i += 1) {
                                ops.add_( // a[i + j * lda] += conj(x[ix]) * temp1 + conj(y[iy]) * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, iy)], ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;

                                ix += incx;
                                iy += incy;
                            }

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[jx]) * temp1 + conj(y[jy]) * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        jy += incy;
                    }
                }
            }
        } else {
            if (noconj) {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[j])
                                ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[j])
                                x[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[j] * temp1 + y[j] * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, j)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, j)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // a[i + j * lda] += x[i] * temp1 + y[i] * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, i)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, i)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * conj(y[jx])
                                ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * x[jx])
                                x[scast(u32, jx)],
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(x[jx] * temp1 + y[jy] * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        x[scast(u32, jx)],
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        y[scast(u32, jy)],
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = jx;
                            var iy: isize = jy;
                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ix += incx;
                                iy += incy;

                                ops.add_( // a[i + j * lda] += x[ix] * temp1 + y[iy] * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            x[scast(u32, ix)],
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            y[scast(u32, iy)],
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        jy += incy;
                    }
                }
            } else {
                if (incx == 1 and incy == 1) {
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, j)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, j)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[j]
                                y[scast(u32, j)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[j]))
                                ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[j]) * temp1 + conj(y[j]) * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, j)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, j)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ops.add_( // a[i + j * lda] += conj(x[i]) * temp1 + conj(y[i]) * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, i)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, i)], ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }
                    }
                } else {
                    var jx: isize = kx;
                    var jy: isize = ky;
                    var j: isize = 0;
                    while (j < n) : (j += 1) {
                        if (ops.ne(x[scast(u32, jx)], 0, ctx) catch unreachable or
                            ops.ne(y[scast(u32, jy)], 0, ctx) catch unreachable)
                        {
                            const temp1: C1 = ops.mul( // temp1 = alpha * y[jy]
                                y[scast(u32, jy)],
                                alpha,
                                ctx,
                            ) catch unreachable;
                            const temp2: C2 = ops.conj(ops.mul( // temp2 = conj(alpha * conj(x[jx]))
                                ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                alpha,
                                ctx,
                            ) catch unreachable, ctx) catch unreachable;

                            ops.add_( // a[j + j * lda] = re(a[j + j * lda]) + re(conj(x[j]) * temp1 + conj(y[j]) * temp2)
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ops.re(ops.add(
                                    ops.mul(
                                        ops.conj(x[scast(u32, jx)], ctx) catch unreachable,
                                        temp1,
                                        ctx,
                                    ) catch unreachable,
                                    ops.mul(
                                        ops.conj(y[scast(u32, jy)], ctx) catch unreachable,
                                        temp2,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable, ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;

                            var ix: isize = jx;
                            var iy: isize = jy;
                            var i: isize = j + 1;
                            while (i < n) : (i += 1) {
                                ix += incx;
                                iy += incy;

                                ops.add_( // a[i + j * lda] += conj(x[ix]) * temp1 + conj(y[iy]) * temp2
                                    &a[scast(u32, i + j * lda)],
                                    a[scast(u32, i + j * lda)],
                                    ops.add(
                                        ops.mul(
                                            ops.conj(x[scast(u32, ix)], ctx) catch unreachable,
                                            temp1,
                                            ctx,
                                        ) catch unreachable,
                                        ops.mul(
                                            ops.conj(y[scast(u32, iy)], ctx) catch unreachable,
                                            temp2,
                                            ctx,
                                        ) catch unreachable,
                                        ctx,
                                    ) catch unreachable,
                                    ctx,
                                ) catch unreachable;
                            }
                        } else {
                            ops.set( // a[i + j * lda] = re(a[i + j * lda])
                                &a[scast(u32, j + j * lda)],
                                ops.re(a[scast(u32, j + j * lda)], ctx) catch unreachable,
                                ctx,
                            ) catch unreachable;
                        }

                        jx += incx;
                        jy += incy;
                    }
                }
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.blas.her2 not implemented for arbitrary precision types yet");
    }

    return;
}
