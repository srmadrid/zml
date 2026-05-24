const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

/// Copies a vector to another vector:
///
/// ```zig
/// y = x,
/// ```
///
/// where `x` and `y` are vectors.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function.
///
/// ## Signature
/// ```zig
/// linalg.blas.copy(n: usize, x: [*]const X, incx: isize, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
///   Must be greater than 0.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`. On return contains the updated vector `y`.
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
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `n`,`incx` or `incy` is equal to
///   0.
pub fn copy(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 8_388_608 / @sizeOf(meta.Child(@TypeOf(x))),
    },
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.copy: x must be a many-item pointer to numerics, and y must be a mutable many-item pointer to numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and !@import("builtin").is_test and X == Y) {
        switch (comptime meta.numericType(X)) {
            .float => {
                if (comptime X == f32)
                    return linalg.cblas.scopy(numeric.cast(isize, n), x, incx, y, incy)
                else if (comptime X == f64)
                    return linalg.cblas.dcopy(numeric.cast(isize, n), x, incx, y, incy);
            },
            .complex => {
                if (comptime meta.Scalar(X) == f32)
                    return linalg.cblas.ccopy(numeric.cast(isize, n), x, incx, y, incy)
                else if (comptime meta.Scalar(X) == f64)
                    return linalg.cblas.zcopy(numeric.cast(isize, n), x, incx, y, incy);
            },
            else => {},
        }
    }

    if (n == 0)
        return;

    if (opts.num_threads == 1)
        return k_copy(n, x, incx, y, incy);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, n / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return k_copy(n, x, incx, y, incy);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return k_copy(n, x, incx, y, incy);

    var threads: [options.max_threads]std.Thread = undefined;

    const chunk_size = int.div(n, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) n else chunk_start + chunk_size;

        if (std.Thread.spawn(.{}, k_copy, .{
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

    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
    }

    if (spawn_err) |err|
        return err;
}

pub fn k_copy(n: usize, x: anytype, incx: isize, y: anytype, incy: isize) void {
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    if (n == 0)
        return;

    const unroll = 2 * (int.min(
        std.simd.suggestVectorLength(X) orelse 2,
        std.simd.suggestVectorLength(Y) orelse 2,
    ));

    if (incx == 1 and incy == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // y[i + u] = x[i + u]
                numeric.set(
                    &y[i + u],
                    x[i + u],
                );
            }
        }

        while (i < n) : (i += 1) {
            // y[i] = x[i]
            numeric.set(
                &y[i],
                x[i],
            );
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        var iy: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // y[iy + u * incy] = x[ix + u * incx]
                numeric.set(
                    &y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                    x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                );
            }

            ix += numeric.cast(isize, unroll) * incx;
            iy += numeric.cast(isize, unroll) * incy;
        }

        while (i < n) : (i += 1) {
            // y[iy] = x[ix]
            numeric.set(
                &y[numeric.cast(usize, iy)],
                x[numeric.cast(usize, ix)],
            );

            ix += incx;
            iy += incy;
        }
    }
}
