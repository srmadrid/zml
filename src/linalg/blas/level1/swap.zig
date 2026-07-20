const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Swaps a vector with another vector.
///
/// ## Signature
/// ```zig
/// linalg.blas.swap(n: usize, x: [*]X, incx: isize, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. On return contains the updated vector `x`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`. On return
///   contains the updated vector `y`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn swap(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.swap: x and y must be mutable many-item pointers to numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return;

    return k_swap(n, x, incx, y, incy);
}

/// Swaps a vector with another vector, splitting the work across the worker
/// threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.swapParallel(n: usize, x: [*]X, incx: isize, y: [*]Y, incy: isize, pool: *thread.Pool) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. On return contains the updated vector `x`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`. On return
///   contains the updated vector `y`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn swapParallel(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    pool: *thread.Pool,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.swap: x and y must be mutable many-item pointers to numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return;

    const Ctx = struct {
        n: usize,
        x: [*]X,
        incx: isize,
        y: [*]Y,
        incy: isize,

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            _ = worker_id;

            k_swap(
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
            );
        }
    };

    pool.parallelFor(
        n,
        Ctx{
            .n = n,
            .x = x,
            .incx = incx,
            .y = y,
            .incy = incy,
        },
        Ctx.kernel,
    );
}

fn k_swap(n: usize, x: anytype, incx: isize, y: anytype, incy: isize) void {
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    const unroll = 2 * (int.min(
        std.simd.suggestVectorLength(X) orelse 2,
        std.simd.suggestVectorLength(Y) orelse 2,
    ));

    if (incx == 1 and incy == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                const temp = x[i + u];

                // x[i + u] = y[i + u]
                numeric.set(
                    &x[i + u],
                    y[i + u],
                );

                // y[i + u] = temp
                numeric.set(
                    &y[i + u],
                    temp,
                );
            }
        }

        while (i < n) : (i += 1) {
            const temp = x[i];

            // x[i] = y[i]
            numeric.set(
                &x[i],
                y[i],
            );

            // y[i] = temp
            numeric.set(
                &y[i],
                temp,
            );
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        var iy: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                const temp = x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)];

                // x[ix + u * incx] = y[iy + u * incy]
                numeric.set(
                    &x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                    y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                );

                // y[iy + u * incy] = temp
                numeric.set(
                    &y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                    temp,
                );
            }

            ix += numeric.cast(isize, unroll) * incx;
            iy += numeric.cast(isize, unroll) * incy;
        }

        while (i < n) : (i += 1) {
            const temp = x[numeric.cast(usize, ix)];

            // x[ix] = y[iy]
            numeric.set(
                &x[numeric.cast(usize, ix)],
                y[numeric.cast(usize, iy)],
            );

            // y[iy] = temp
            numeric.set(
                &y[numeric.cast(usize, iy)],
                temp,
            );

            ix += incx;
            iy += incy;
        }
    }
}
