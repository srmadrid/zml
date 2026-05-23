const std = @import("std");
const options = @import("options");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

pub fn Asum(X: type) type {
    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.asum: x must be a many-item pointer to numerics, got \n\tx: " ++ @typeName(X) ++ "\n");

    return numeric.Abs1(meta.Child(X));
}

/// Computes the sum of magnitudes of the elements of a real vector, or the sum
/// of magnitudes of the real and imaginary parts of elements of a complex
/// vector:
///
/// ```zig
/// abs(x[0].re) + abs(x[0].im) + abs(x[1].re) + abs(x[1].im) + ... + abs(x[n - 1].re) + abs(x[n - 1].im),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function.
///
/// ## Signature
/// ```zig
/// linalg.blas.asum(n: usize, x: [*]const X, incx: isize, opts: Opts) !linalg.blas.Asum([*]const X)
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vector `x`. Must be
///   greater than 0.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
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
/// `Asum(@TypeOf(x))`: The sum of magnitudes of real and imaginary parts of all
/// elements of the vector.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is equal to 0.
pub fn asum(
    n: usize,
    x: anytype,
    incx: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 8_388_608 / @sizeOf(meta.Child(@TypeOf(x))),
    },
) !linalg.blas.Asum(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and !@import("builtin").is_test) {
        switch (comptime meta.numericType(meta.Child(X))) {
            .float => {
                if (comptime meta.Child(X) == f32)
                    return linalg.cblas.sasum(numeric.cast(isize, n), x, incx)
                else if (comptime meta.Child(X) == f64)
                    return linalg.cblas.dasum(numeric.cast(isize, n), x, incx);
            },
            .complex => {
                if (comptime meta.Scalar(meta.Child(X)) == f32)
                    return linalg.cblas.scasum(numeric.cast(isize, n), x, incx)
                else if (comptime meta.Scalar(meta.Child(X)) == f64)
                    return linalg.cblas.dzasum(numeric.cast(isize, n), x, incx);
            },
            else => {},
        }
    }

    if (n == 0)
        return numeric.zero(linalg.blas.Asum(X));

    if (opts.num_threads == 1)
        return numeric.cast(linalg.blas.Asum(X), k_asum(n, x, incx));

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, n / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return numeric.cast(linalg.blas.Asum(X), k_asum(n, x, incx));

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return numeric.cast(linalg.blas.Asum(X), k_asum(n, x, incx));

    var threads: [options.max_threads]std.Thread = undefined;
    var sums: [options.max_threads]meta.Accumulator(linalg.blas.Asum(X)) = .{numeric.zero(meta.Accumulator(linalg.blas.Asum(X)))} ** options.max_threads;

    const Worker = struct {
        fn execute(out: *meta.Accumulator(linalg.blas.Asum(X)), worker_n: usize, worker_x: X, worker_incx: isize) void {
            out.* = k_asum(worker_n, worker_x, worker_incx);
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
        })) |th| {
            threads[i] = th;
            spawned_count += 1;
        } else |err| {
            spawn_err = err;
            break;
        }
    }

    var sum = numeric.zero(meta.Accumulator(linalg.blas.Asum(X)));
    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
        numeric.add_(&sum, sum, sums[t]);
    }

    if (spawn_err) |err|
        return err;

    return numeric.cast(linalg.blas.Asum(X), sum);
}

fn k_asum(n: usize, x: anytype, incx: isize) meta.Accumulator(linalg.blas.Asum(@TypeOf(x))) {
    const X: type = @TypeOf(x);

    if (n == 0)
        return numeric.zero(linalg.blas.Asum(X));

    const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(linalg.blas.Asum(X))) orelse 2);

    var sums: [unroll]meta.Accumulator(linalg.blas.Asum(X)) = .{numeric.zero(meta.Accumulator(linalg.blas.Asum(X)))} ** unroll;

    if (incx == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // sums[u] += abs1(x[i + u])
                numeric.add_(
                    &sums[u],
                    sums[u],
                    numeric.abs1(x[i + u]),
                );
            }
        }

        while (i < n) : (i += 1) {
            // sums[0] += abs1(x[i])
            numeric.add_(
                &sums[0],
                sums[0],
                numeric.abs1(x[i]),
            );
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // sums[u] += abs1(x[ix + u * incx])
                numeric.add_(
                    &sums[u],
                    sums[u],
                    numeric.abs1(x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)]),
                );
            }

            ix += numeric.cast(isize, unroll) * incx;
        }

        while (i < n) : (i += 1) {
            // sums[0] += abs1(x[ix])
            numeric.add_(
                &sums[0],
                sums[0],
                numeric.abs1(x[numeric.cast(usize, ix)]),
            );

            ix += incx;
        }
    }

    var sum = numeric.zero(meta.Accumulator(linalg.blas.Asum(X)));
    inline for (0..unroll) |u| {
        // sum += sums[u]
        numeric.add_(&sum, sum, sums[u]);
    }

    return sum;
}
