const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Computes a vector-scalar product and adds the result to a vector:
///
/// ```zig
/// y = alpha * x + y,
/// ```
///
/// where `alpha` is a numeric, and `x` and `y` are vectors each with `n`
/// elements.
///
/// ## Signature
/// ```zig
/// linalg.blas.axpy(n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`. On return contains the updated vector `y`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn axpy(
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.axpy: alpha must be a numeric, x must be a many-item pointer to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0 or numeric.eq(alpha, 0))
        return;

    return k_axpy(n, alpha, x, incx, y, incy);
}

/// Computes a vector-scalar product and adds the result to a vector:
///
/// ```zig
/// y = alpha * x + y,
/// ```
///
/// where `alpha` is a numeric, and `x` and `y` are vectors each with `n`
/// elements, splitting the work across the worker threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.axpy(n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]Y, incy: isize, pool: *thread.Pool) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`. On return contains the updated vector `y`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn axpyParallel(
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    pool: *thread.Pool,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.axpyParallel: alpha must be a numeric, x must be a many-item pointer to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0 or numeric.eq(alpha, 0))
        return;

    const Ctx = struct {
        n: usize,
        alpha: Al,
        x: [*]const X,
        incx: isize,
        y: [*]Y,
        incy: isize,

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            _ = worker_id;

            k_axpy(
                end - start,
                ctx.alpha,
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
            .alpha = alpha,
            .x = x,
            .incx = incx,
            .y = y,
            .incy = incy,
        },
        Ctx.kernel,
    );
}

fn k_axpy(n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize) void {
    const Al: type = @TypeOf(alpha);
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(Al, X, Y)) orelse 2);

    if (incx == 1 and incy == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // y[i + u] += alpha * x[i + u]
                numeric.fmaInto(
                    &y[i + u],
                    alpha,
                    x[i + u],
                    y[i + u],
                );
            }
        }

        while (i < n) : (i += 1) {
            // y[i] += alpha * x[i]
            numeric.fmaInto(
                &y[i],
                alpha,
                x[i],
                y[i],
            );
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;
        var iy: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // y[iy + u * incy] += alpha * x[ix + u * incx]
                numeric.fmaInto(
                    &y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                    alpha,
                    x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                    y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                );
            }

            ix += numeric.cast(isize, unroll) * incx;
            iy += numeric.cast(isize, unroll) * incy;
        }

        while (i < n) : (i += 1) {
            // y[iy] += alpha * x[ix]
            numeric.fmaInto(
                &y[numeric.cast(usize, iy)],
                alpha,
                x[numeric.cast(usize, ix)],
                y[numeric.cast(usize, iy)],
            );

            ix += incx;
            iy += incy;
        }
    }
}
