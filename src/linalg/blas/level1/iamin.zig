const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

/// Finds the (0-based) index of the element with the smallest magnitude:
///
/// ```zig
/// argmin_i abs1(x[i]),   for i in 0..n
/// ```
///
/// If multiple elements share the minimum value, the smallest index is
/// returned.
///
/// ## Signature
/// ```zig
/// linalg.blas.iamin(n: usize, x: [*]const X, incx: isize) !usize
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Number of elements in `x`. Must be greater than 0.
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
///       `std.Thread.getCpuCount()` and `options.max_threads` as a hard
///       safety ceiling. `parallel_threshold` is ignored.
///   * `parallel_threshold` (`usize = 2_097_152 / @sizeOf(meta.Child(X))`):
///     Minimum number of elements required to trigger multithreaded
///     execution.
///
/// ## Returns
/// `usize`: The 0-based index of the first element with the smallest magnitude.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is equal to 0.
pub fn iamin(
    n: usize,
    x: anytype,
    incx: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 2_097_152 / @sizeOf(meta.Child(@TypeOf(x))),
    },
) !usize {
    const X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.iamin: x must be a many-item pointer to numerics, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return 0;

    if (opts.num_threads == 1)
        return k_iamin(n, x, incx).index;

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, n / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return k_iamin(n, x, incx).index;

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return k_iamin(n, x, incx).index;

    var threads: [options.max_threads]std.Thread = undefined;
    var results: [options.max_threads]IaminResult(numeric.Abs1(meta.Child(X))) = .{IaminResult(numeric.Abs1(meta.Child(X))){ .value = numeric.zero(numeric.Abs1(meta.Child(X))), .index = 0 }} ** options.max_threads;
    var chunk_bases: [options.max_threads]usize = .{0} ** options.max_threads;

    const Worker = struct {
        fn execute(out: *IaminResult(numeric.Abs1(meta.Child(X))), worker_n: usize, worker_x: X, worker_incx: isize) void {
            out.* = k_iamin(worker_n, worker_x, worker_incx);
        }
    };

    const chunk_size = int.div(n, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) n else chunk_start + chunk_size;

        chunk_bases[i] = chunk_start;

        if (std.Thread.spawn(.{}, Worker.execute, .{
            &results[i],
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

    var best: IaminResult(numeric.Abs1(meta.Child(X))) = .{ .value = numeric.zero(numeric.Abs1(meta.Child(X))), .index = 0 };
    var best_found = false;
    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
        if (!best_found or numeric.lt(results[t].value, best.value)) {
            best = .{ .value = results[t].value, .index = chunk_bases[t] + results[t].index };
            best_found = true;
        }
    }

    if (spawn_err) |err|
        return err;

    return best.index;
}

pub fn IaminResult(N: type) type {
    return struct {
        value: N,
        index: usize,
    };
}

fn k_iamin(n: usize, x: anytype, incx: isize) IaminResult(numeric.Abs1(meta.Child(@TypeOf(x)))) {
    if (n == 0)
        return .{ .value = numeric.zero(numeric.Abs1(meta.Child(@TypeOf(x)))), .index = 0 };

    var best_value = if (incx == 1)
        numeric.abs1(x[0])
    else
        numeric.abs1(x[numeric.cast(usize, if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0)]);
    var best_index: usize = 0;

    if (incx == 1) {
        var i: usize = 1;
        while (i < n) : (i += 1) {
            const temp = numeric.abs1(x[i]);
            if (numeric.lt(temp, best_value)) {
                best_value = temp;
                best_index = i;
            }
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

        ix += incx;

        var i: usize = 1;
        while (i < n) : (i += 1) {
            const temp = numeric.abs1(x[numeric.cast(usize, ix)]);
            if (numeric.lt(temp, best_value)) {
                best_value = temp;
                best_index = i;
            }
            ix += incx;
        }
    }

    return .{ .value = best_value, .index = best_index };
}
