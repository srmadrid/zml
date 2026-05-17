const std = @import("std");
const options = @import("options");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Computes the product of a vector by a numeric:
///
/// ```zig
/// x = alpha * x,
/// ```
///
/// where `alpha` is a numeric, and `x` is a vector with `n` elements.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function.
///
/// ## Signature
/// ```zig
/// linalg.blas.scal(n: usize, alpha: Al, x: [*]X, incx: isize) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
///   Must be greater than 0.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
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
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is equal to 0.
pub fn scal(
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 8_388_608 / @sizeOf(meta.Child(@TypeOf(x))),
    },
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.scal: alpha must be a numeric, and x must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    X = meta.Child(X);

    if (n == 0 or incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == X) {
        switch (comptime meta.numericType(Al)) {
            .float => {
                if (comptime Al == f32)
                    return linalg.cblas.sscal(numeric.cast(isize, n), alpha, x, incx)
                else if (comptime Al == f64)
                    return linalg.cblas.dscal(numeric.cast(isize, n), alpha, x, incx);
            },
            .complex => {
                if (comptime meta.Scalar(Al) == f32)
                    return linalg.cblas.cscal(numeric.cast(isize, n), &alpha, x, incx)
                else if (comptime meta.Scalar(Al) == f64)
                    return linalg.cblas.zscal(numeric.cast(isize, n), &alpha, x, incx);
            },
            else => {},
        }
    }

    if (opts.num_threads == 1)
        return k_scal(n, alpha, x, incx);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, n / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return k_scal(n, alpha, x, incx);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return k_scal(n, alpha, x, incx);

    var threads: [options.max_threads]std.Thread = undefined;

    const chunk_size = int.div(n, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) n else chunk_start + chunk_size;

        if (std.Thread.spawn(.{}, k_scal, .{
            chunk_end - chunk_start,
            alpha,
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

    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
    }

    if (spawn_err) |err|
        return err;
}

pub fn k_scal(n: usize, alpha: anytype, x: anytype, incx: isize) void {
    const Al: type = @TypeOf(alpha);
    const X: type = meta.Child(@TypeOf(x));

    const unroll = 2 * (std.simd.suggestVectorLength(numeric.Mul(Al, X)) orelse 2);

    if (incx == 1) {
        if (numeric.eq(alpha, 0)) {
            var i: usize = 0;
            while (i < (n / unroll) * unroll) : (i += unroll) {
                inline for (0..unroll) |u| {
                    // x[i + u] = 0
                    x[i + u] = numeric.zero(X);
                }
            }

            while (i < n) : (i += 1) {
                // x[i] = 0
                x[i] = numeric.zero(X);
            }
        } else {
            var i: usize = 0;
            while (i < (n / unroll) * unroll) : (i += unroll) {
                inline for (0..unroll) |u| {
                    // x[i + u] *= alpha
                    numeric.mul_(
                        &x[i + u],
                        alpha,
                        x[i + u],
                    );
                }
            }

            while (i < n) : (i += 1) {
                // x[i] *= alpha
                numeric.mul_(
                    &x[i],
                    alpha,
                    x[i],
                );
            }
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        if (numeric.eq(alpha, 0)) {
            var i: usize = 0;
            while (i < (n / unroll) * unroll) : (i += unroll) {
                inline for (0..unroll) |u| {
                    // x[ix + u * incx] = 0
                    x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)] = numeric.zero(X);
                }

                ix += numeric.cast(isize, unroll) * incx;
            }

            while (i < n) : (i += 1) {
                // x[ix] *= alpha
                x[numeric.cast(usize, ix)] = numeric.zero(X);

                ix += incx;
            }
        } else {
            var i: usize = 0;
            while (i < (n / unroll) * unroll) : (i += unroll) {
                inline for (0..unroll) |u| {
                    // x[ix + u * incx] *= alpha
                    numeric.mul_(
                        &x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                        alpha,
                        x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                    );
                }

                ix += numeric.cast(isize, unroll) * incx;
            }

            while (i < n) : (i += 1) {
                // x[ix] *= alpha
                numeric.mul_(
                    &x[numeric.cast(usize, ix)],
                    alpha,
                    x[numeric.cast(usize, ix)],
                );

                ix += incx;
            }
        }
    }
}
