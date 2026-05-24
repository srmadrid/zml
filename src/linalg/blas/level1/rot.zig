const std = @import("std");
const options = @import("options");

const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

const int = @import("../../../int.zig");

const linalg = @import("../../../linalg.zig");

/// Applies a Givens plane rotation to the vectors `x` and `y`:
///
/// ```zig
/// x[i] = c * x[i] + s * y[i]
/// y[i] = c * y[i] - s * x[i]
/// ```
///
/// for `i in 0..n`. `c` is a real numeric, but `s` may be either a real or a
/// complex numeric.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function.
///
/// ## Signature
/// ```zig
/// linalg.blas.rot(n: usize, x: [*]X, incx: isize, y: [*]Y, incy: isize, c: C, s: S) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Number of elements in `x` and `y`. Must be greater than 0.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `c` (`anytype`): Cosine of the rotation angle. Real numeric.
/// * `s` (`anytype`): Sine of the rotation angle. Real or complex numeric.
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
///   * `parallel_threshold` (`usize = 4_194_304 / @sizeOf(meta.Child(X))`):
///     Minimum number of elements required to trigger multithreaded
///     execution.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `n`, `incx` or `incy` is equal to
///   zero.
pub fn rot(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    c: anytype,
    s: anytype,
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 4_194_304 / @sizeOf(meta.Child(@TypeOf(x))),
    },
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    const C: type = @TypeOf(c);
    const S: type = @TypeOf(s);

    comptime if (!meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)) or
        !meta.isNumeric(C) or !meta.isReal(C) or !meta.isNumeric(S))
        @compileError("zsl.linalg.blas.rot: x and y must be mutable many-item pointers to numerics, c must be a real numeric, and s must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tc: " ++ @typeName(C) ++ "\n\ts: " ++ @typeName(S) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);

    if (n == 0 or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and X == Y and meta.Real(X) == C and meta.Real(X) == S) {
        switch (comptime meta.numericType(X)) {
            .float => {
                if (comptime X == f32)
                    return linalg.cblas.srot(numeric.cast(isize, n), x, incx, y, incy, c, s)
                else if (comptime X == f64)
                    return linalg.cblas.drot(numeric.cast(isize, n), x, incx, y, incy, c, s);
            },
            .complex => {
                if (comptime !meta.isComplex(S)) {
                    if (comptime meta.Scalar(X) == f32)
                        return linalg.cblas.csrot(numeric.cast(isize, n), x, incx, y, incy, c, s)
                    else if (comptime meta.Scalar(X) == f64)
                        return linalg.cblas.zdrot(numeric.cast(isize, n), x, incx, y, incy, c, s);
                } else {
                    if (comptime meta.Scalar(X) == f32)
                        return linalg.cblas.crot(numeric.cast(isize, n), x, incx, y, incy, c, &s)
                    else if (comptime meta.Scalar(X) == f64)
                        return linalg.cblas.zrot(numeric.cast(isize, n), x, incx, y, incy, c, &s);
                }
            },
            else => {},
        }
    }

    if (opts.num_threads == 1)
        return k_rot(n, x, incx, y, incy, c, s);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, n / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return k_rot(n, x, incx, y, incy, c, s);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return k_rot(n, x, incx, y, incy, c, s);

    var threads: [options.max_threads]std.Thread = undefined;

    const chunk_size = int.div(n, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) n else chunk_start + chunk_size;

        if (std.Thread.spawn(.{}, k_rot, .{
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
            c,
            s,
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

fn k_rot(n: usize, x: anytype, incx: isize, y: anytype, incy: isize, c: anytype, s: anytype) void {
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));
    const C: type = @TypeOf(c);
    const S: type = @TypeOf(s);

    const unroll = 2 * int.min(
        std.simd.suggestVectorLength(numeric.Fma(C, X, numeric.Mul(S, Y))) orelse 2,
        std.simd.suggestVectorLength(numeric.Fma(C, Y, numeric.Neg(numeric.Mul(numeric.Conj(S), X)))) orelse 2,
    );

    if (incx == 1 and incy == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                const xi = x[i + u];

                // x[i + u] = c * xi + s * y[i + u]
                numeric.fma_(
                    &x[i + u],
                    c,
                    xi,
                    numeric.mul(s, y[i + u]),
                );

                // y[i + u] = c * y[i + u] - conj(s) * xi
                numeric.fma_(
                    &y[i + u],
                    c,
                    y[i + u],
                    numeric.neg(numeric.mul(numeric.conj(s), xi)),
                );
            }
        }

        while (i < n) : (i += 1) {
            const xi = x[i];

            // x[i] = c * xi + s * y[i]
            numeric.fma_(
                &x[i],
                c,
                xi,
                numeric.mul(s, y[i]),
            );

            // y[i] = c * y[i] - conj(s) * xi
            numeric.fma_(
                &y[i],
                c,
                y[i],
                numeric.neg(numeric.mul(numeric.conj(s), xi)),
            );
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        var iy: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                const x_idx = numeric.cast(usize, ix + numeric.cast(isize, u) * incx);
                const y_idx = numeric.cast(usize, iy + numeric.cast(isize, u) * incy);
                const xi = x[x_idx];

                // x[x_idx] = c * xi + s * y[y_idx]
                numeric.fma_(
                    &x[x_idx],
                    c,
                    xi,
                    numeric.mul(s, y[y_idx]),
                );

                // y[y_idx] = c * y[y_idx] - conj(s) * xi
                numeric.fma_(
                    &y[y_idx],
                    c,
                    y[y_idx],
                    numeric.neg(numeric.mul(numeric.conj(s), xi)),
                );
            }

            ix += numeric.cast(isize, unroll) * incx;
            iy += numeric.cast(isize, unroll) * incy;
        }

        while (i < n) : (i += 1) {
            const x_idx = numeric.cast(usize, ix);
            const y_idx = numeric.cast(usize, iy);
            const xi = x[x_idx];

            // x[x_idx] = c * xi + s * y[y_idx]
            numeric.fma_(
                &x[x_idx],
                c,
                xi,
                numeric.mul(s, y[y_idx]),
            );

            // y[y_idx] = c * y[y_idx] - conj(s) * xi
            numeric.fma_(
                &y[y_idx],
                c,
                y[y_idx],
                numeric.neg(numeric.mul(numeric.conj(s), xi)),
            );

            ix += incx;
            iy += incy;
        }
    }
}
