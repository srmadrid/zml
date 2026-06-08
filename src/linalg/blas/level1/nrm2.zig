const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

pub fn Nrm2(X: type) type {
    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.Nrm2: X must be a many-item pointer type to numerics, got \n\tX = " ++ @typeName(X) ++ "\n");

    return numeric.Sqrt(numeric.Abs2(meta.Child(X)));
}

/// Computes the Euclidean norm of a vector:
///
/// ```zig
/// sqrt(|x[0]|^2 + |x[1]|^2 + ... + |x[n - 1]|^2),
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// ## Signature
/// ```zig
/// linalg.blas.nrm2(n: usize, x: [*]const X, incx: isize) !linalg.blas.Nrm2([*]const X)
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
///       `std.Thread.getCpuCount()` and `options.max_threads` as a hard
///       safety ceiling. `parallel_threshold` is ignored.
///   * `parallel_threshold` (`usize = 8_388_608 / @sizeOf(meta.Child(X))`):
///     Minimum number of elements required to trigger multithreaded
///     execution.
///
/// ## Returns
/// `Nrm2(@TypeOf(x))`: The Euclidean norm of `x`.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `n` or `incx` is equal to 0.
pub fn nrm2(
    n: usize,
    x: anytype,
    incx: isize,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 8_388_608 / @sizeOf(meta.Child(@TypeOf(x))),
    },
) !linalg.blas.Nrm2(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return numeric.zero(linalg.blas.Nrm2(X));

    if (opts.num_threads == 1)
        return numeric.cast(linalg.blas.Nrm2(X), numeric.sqrt(k_nrm2_ssq(n, x, incx)));

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, n / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return numeric.cast(linalg.blas.Nrm2(X), numeric.sqrt(k_nrm2_ssq(n, x, incx)));

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return numeric.cast(linalg.blas.Nrm2(X), numeric.sqrt(k_nrm2_ssq(n, x, incx)));

    const Worker = struct {
        fn execute(out: *meta.Accumulator(linalg.blas.Nrm2(X)), worker_n: usize, worker_x: X, worker_incx: isize) void {
            out.* = k_nrm2_ssq(worker_n, worker_x, worker_incx);
        }
    };

    var threads: [options.max_threads]std.Thread = undefined;
    var sums: [options.max_threads]meta.Accumulator(linalg.blas.Nrm2(X)) = .{numeric.zero(meta.Accumulator(linalg.blas.Nrm2(X)))} ** options.max_threads;

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

    var ssq = numeric.zero(meta.Accumulator(linalg.blas.Nrm2(X)));
    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
        numeric.addInto(&ssq, ssq, sums[t]);
    }

    if (spawn_err) |err|
        return err;

    return numeric.cast(linalg.blas.Nrm2(X), numeric.sqrt(ssq));
}

fn k_nrm2_ssq(n: usize, x: anytype, incx: isize) meta.Accumulator(linalg.blas.Nrm2(@TypeOf(x))) {
    const X: type = @TypeOf(x);

    if (n == 0)
        return numeric.zero(linalg.blas.Nrm2(X));

    const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(linalg.blas.Nrm2(X))) orelse 2);

    var sums: [unroll]meta.Accumulator(linalg.blas.Nrm2(X)) = .{numeric.zero(meta.Accumulator(linalg.blas.Nrm2(X)))} ** unroll;

    if (incx == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                if (comptime meta.isComplex(meta.Child(X))) {
                    // sums[u] += re(x[i + u])^2
                    numeric.fmaInto(
                        &sums[u],
                        numeric.re(x[i + u]),
                        numeric.re(x[i + u]),
                        sums[u],
                    );

                    // sums[u] += im(x[i + u])^2
                    numeric.fmaInto(
                        &sums[u],
                        numeric.im(x[i + u]),
                        numeric.im(x[i + u]),
                        sums[u],
                    );
                } else {
                    // sums[u] += x[i + u]^2
                    numeric.fmaInto(
                        &sums[u],
                        x[i + u],
                        x[i + u],
                        sums[u],
                    );
                }
            }
        }

        while (i < n) : (i += 1) {
            if (comptime meta.isComplex(meta.Child(X))) {
                // sums[0] += re(x[i])^2
                numeric.fmaInto(
                    &sums[0],
                    numeric.re(x[i]),
                    numeric.re(x[i]),
                    sums[0],
                );

                // sums[0] += im(x[i])^2
                numeric.fmaInto(
                    &sums[0],
                    numeric.im(x[i]),
                    numeric.im(x[i]),
                    sums[0],
                );
            } else {
                // sums[0] += x[i]^2
                numeric.fmaInto(
                    &sums[0],
                    x[i],
                    x[i],
                    sums[0],
                );
            }
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                const idx = numeric.cast(usize, ix + numeric.cast(isize, u) * incx);
                if (comptime meta.isComplex(meta.Child(X))) {
                    // sums[0] += re(x[idx])^2
                    numeric.fmaInto(
                        &sums[0],
                        numeric.re(x[idx]),
                        numeric.re(x[idx]),
                        sums[0],
                    );

                    // sums[0] += im(x[idx])^2
                    numeric.fmaInto(
                        &sums[0],
                        numeric.im(x[idx]),
                        numeric.im(x[idx]),
                        sums[0],
                    );
                } else {
                    // sums[0] += x[idx]^2
                    numeric.fmaInto(
                        &sums[0],
                        x[idx],
                        x[idx],
                        sums[0],
                    );
                }
            }

            ix += numeric.cast(isize, unroll) * incx;
        }

        while (i < n) : (i += 1) {
            const idx = numeric.cast(usize, ix);
            if (comptime meta.isComplex(meta.Child(X))) {
                // sums[0] += re(x[idx])^2
                numeric.fmaInto(
                    &sums[0],
                    numeric.re(x[idx]),
                    numeric.re(x[idx]),
                    sums[0],
                );

                // sums[0] += im(x[idx])^2
                numeric.fmaInto(
                    &sums[0],
                    numeric.im(x[idx]),
                    numeric.im(x[idx]),
                    sums[0],
                );
            } else {
                // sums[0] += x[idx]^2
                numeric.fmaInto(
                    &sums[0],
                    x[idx],
                    x[idx],
                    sums[0],
                );
            }

            ix += incx;
        }
    }

    var ssq = numeric.zero(meta.Accumulator(linalg.blas.Nrm2(X)));
    inline for (0..unroll) |u| {
        numeric.addInto(&ssq, ssq, sums[u]);
    }

    return ssq;
}
