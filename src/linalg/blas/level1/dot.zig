const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

pub fn Dot(X: type, Y: type) type {
    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.dot: x and y must be many-item pointers to numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.Mul(meta.Child(X), meta.Child(Y));
}

/// Computes a vector dot product:
///
/// ```zig
/// x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function.
///
/// ## Signature
/// ```zig
/// linalg.blas.dot(n: isize, x: [*]const X, incx: isize, y: [*]const Y, incy: isize) !linalg.blas.Dot([*]const X, [*]const Y)
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
///   Must be greater than 0.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `opts`: Optional parameters:
///   * `num_threads` (`usize = 0`): Number of threads to spawn:
///     * `0`: automatic. The thread count is derived from `n` and
///       `parallel_threshold`:
///       ```zig
///       threads = max(1, min(std.Thread.getCpuCount(), options.max_threads, n / parallel_threshold))
///       ```
///     * `1`: force serial execution. `parallel_threshold` is ignored.
///     * `N >= 2`: use exactly `N` threads, clamped by
///       `std.Thread.getCpuCount()` and`options.max_threads` as a hard safety
///       ceiling. `parallel_threshold` is ignored.
///   * `parallel_threshold` (`usize = 8_388_608 / @sizeOf(meta.Child(X))`):
///     Minimum number of elements required to trigger multithreaded execution.
///
/// ## Returns
/// `Dot(@TypeOf(x), @TypeOf(y))`: The dot product of `x` and `y`.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `n`, `incx` or `incy` is equal to
///   0.
pub fn dot(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 8_388_608 / @sizeOf(meta.Child(@TypeOf(x))),
    },
) !linalg.blas.Dot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if ((comptime options.link_cblas != null) and !@import("builtin").is_test and meta.Child(X) == meta.Child(Y)) {
        switch (comptime meta.numericType(meta.Child(X))) {
            .float => {
                if (comptime meta.Child(X) == f32)
                    return linalg.cblas.sdot(numeric.cast(isize, n), x, incx, y, incy)
                else if (comptime meta.Child(X) == f64)
                    return linalg.cblas.ddot(numeric.cast(isize, n), x, incx, y, incy);
            },
            .complex => {
                var result: linalg.blas.Dot(X, Y) = undefined;
                if (comptime meta.Scalar(meta.Child(X)) == f32)
                    linalg.cblas.cdotu_sub(numeric.cast(isize, n), x, incx, y, incy, &result)
                else if (comptime meta.Scalar(meta.Child(X)) == f64)
                    linalg.cblas.zdotu_sub(numeric.cast(isize, n), x, incx, y, incy, &result);

                return result;
            },
            else => {},
        }
    }

    if (n == 0)
        return numeric.zero(linalg.blas.Dot(X, Y));

    if (opts.num_threads == 1)
        return numeric.cast(linalg.blas.Dot(X, Y), k_dot(n, x, incx, y, incy));

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, n / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return numeric.cast(linalg.blas.Dot(X, Y), k_dot(n, x, incx, y, incy));

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return numeric.cast(linalg.blas.Dot(X, Y), k_dot(n, x, incx, y, incy));

    var threads: [options.max_threads]std.Thread = undefined;
    var sums: [options.max_threads]meta.Accumulator(linalg.blas.Dot(X, Y)) = .{numeric.zero(meta.Accumulator(linalg.blas.Dot(X, Y)))} ** options.max_threads;

    const Worker = struct {
        fn execute(out: *meta.Accumulator(linalg.blas.Dot(X, Y)), worker_n: usize, worker_x: X, worker_incx: isize, worker_y: Y, worker_incy: isize) void {
            out.* = k_dot(worker_n, worker_x, worker_incx, worker_y, worker_incy);
        }
    };

    const chunk_size = int.div(n, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) n else chunk_start + chunk_size;

        if (std.Thread.spawn(.{}, Worker.execute, .{
            &sums[i],
            chunk_end - chunk_start,
            x + numeric.cast(usize, if (incx > 0)
                numeric.cast(isize, chunk_start) * incx
            else
                (-numeric.cast(isize, n) + numeric.cast(isize, chunk_end)) * incx),
            incx,
            y + numeric.cast(usize, if (incy > 0)
                numeric.cast(isize, chunk_start) * incy
            else
                (-numeric.cast(isize, n) + numeric.cast(isize, chunk_end)) * incy),
            incy,
        })) |th| {
            threads[i] = th;
            spawned_count += 1;
        } else |err| {
            spawn_err = err;
            break;
        }
    }

    var sum = numeric.zero(meta.Accumulator(linalg.blas.Dot(X, Y)));
    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
        numeric.add_(&sum, sum, sums[t]);
    }

    if (spawn_err) |err|
        return err;

    return numeric.cast(linalg.blas.Dot(X, Y), sum);
}

fn k_dot(n: usize, x: anytype, incx: isize, y: anytype, incy: isize) meta.Accumulator(linalg.blas.Dot(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (n == 0)
        return numeric.zero(linalg.blas.Dot(X, Y));

    const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(linalg.blas.Dot(X, Y))) orelse 2);

    var sums: [unroll]meta.Accumulator(linalg.blas.Dot(X, Y)) = .{numeric.zero(meta.Accumulator(linalg.blas.Dot(X, Y)))} ** unroll;

    if (incx == 1 and incy == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // sums[u] += x[i + u] * y[i + u]
                numeric.fma_(
                    &sums[u],
                    x[i + u],
                    y[i + u],
                    sums[u],
                );
            }
        }

        while (i < n) : (i += 1) {
            // sums[0] += x[i] * y[i]
            numeric.fma_(
                &sums[0],
                x[i],
                y[i],
                sums[0],
            );
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        var iy: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // sums[u] += x[ix + u * incx] * y[iy + u * incy]
                numeric.fma_(
                    &sums[u],
                    x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                    y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                    sums[u],
                );
            }

            ix += numeric.cast(isize, unroll) * incx;
            iy += numeric.cast(isize, unroll) * incy;
        }

        while (i < n) : (i += 1) {
            // sums[0] += x[ix] * y[ix]
            numeric.fma_(
                &sums[0],
                x[numeric.cast(usize, ix)],
                y[numeric.cast(usize, iy)],
                sums[0],
            );

            ix += incx;
            iy += incy;
        }
    }

    var sum = numeric.zero(meta.Accumulator(linalg.blas.Dot(X, Y)));
    inline for (0..unroll) |u| {
        // sum += sums[u]
        numeric.add_(&sum, sum, sums[u]);
    }

    return sum;
}
