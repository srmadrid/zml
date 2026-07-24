const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

pub fn Dot(X: type, Y: type) type {
    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.Dot: X and Y must be many-item pointer types to numerics, got \n\tX = " ++ @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

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
/// ## Signature
/// ```zig
/// linalg.blas.dot(n: isize, x: [*]const X, incx: isize, y: [*]const Y, incy: isize) !linalg.blas.Dot([*]const X, [*]const Y)
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
///
/// ## Returns
/// `Dot(@TypeOf(x), @TypeOf(y))`: The dot product of `x` and `y`.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn dot(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
) !linalg.blas.Dot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return numeric.zero(linalg.blas.Dot(X, Y));

    return numeric.cast(linalg.blas.Dot(X, Y), k_dot(n, x, incx, y, incy));
}

/// Computes a vector dot product:
///
/// ```zig
/// x[0] * y[0] + x[1] * y[1] + ... + x[n - 1] * y[n - 1],
/// ```
///
/// where `x` and `y` are vectors, splitting the work across the worker threads
/// of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.dotParallel(n: isize, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, pool: *thread.Pool) !linalg.blas.Dot([*]const X, [*]const Y)
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `Dot(@TypeOf(x), @TypeOf(y))`: The dot product of `x` and `y`.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn dotParallel(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    pool: *thread.Pool,
) !linalg.blas.Dot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return numeric.zero(linalg.blas.Dot(X, Y));

    const Ctx = struct {
        n: usize,
        x: X,
        incx: isize,
        y: Y,
        incy: isize,
        sums: *[thread.max_workers]meta.Accumulator(linalg.blas.Dot(X, Y)),

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            numeric.addInto(
                &ctx.sums[worker_id],
                ctx.sums[worker_id],
                k_dot(
                    end - start,
                    ctx.x + numeric.cast(usize, if (ctx.incx > 0)
                        numeric.cast(isize, start) * ctx.incx
                    else
                        (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incx),
                    ctx.incx,
                    ctx.y + numeric.cast(usize, if (ctx.incy > 0)
                        numeric.cast(isize, start) * ctx.incy
                    else
                        (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incy),
                    ctx.incy,
                ),
            );
        }
    };

    var sums: [thread.max_workers]meta.Accumulator(linalg.blas.Dot(X, Y)) = @splat(numeric.zero(meta.Accumulator(linalg.blas.Dot(X, Y))));

    pool.parallelFor(
        n,
        Ctx{
            .n = n,
            .x = x,
            .incx = incx,
            .y = y,
            .incy = incy,
            .sums = &sums,
        },
        Ctx.kernel,
    );

    var sum = numeric.zero(meta.Accumulator(linalg.blas.Dot(X, Y)));
    for (0..int.min(pool.idCount(), thread.max_workers)) |i| {
        numeric.addInto(&sum, sum, sums[i]);
    }

    return numeric.cast(linalg.blas.Dot(X, Y), sum);
}

fn k_dot(n: usize, x: anytype, incx: isize, y: anytype, incy: isize) meta.Accumulator(linalg.blas.Dot(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(linalg.blas.Dot(X, Y))) orelse 2);

    var sums: [unroll]meta.Accumulator(linalg.blas.Dot(X, Y)) = .{numeric.zero(meta.Accumulator(linalg.blas.Dot(X, Y)))} ** unroll;

    if (incx == 1 and incy == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // sums[u] += x[i + u] * y[i + u]
                numeric.fmaInto(
                    &sums[u],
                    x[i + u],
                    y[i + u],
                    sums[u],
                );
            }
        }

        while (i < n) : (i += 1) {
            // sums[0] += x[i] * y[i]
            numeric.fmaInto(
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
                numeric.fmaInto(
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
            numeric.fmaInto(
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
        numeric.addInto(&sum, sum, sums[u]);
    }

    return sum;
}
