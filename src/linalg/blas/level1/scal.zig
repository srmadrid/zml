const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Computes the product of a vector by a numeric:
///
/// ```zig
/// x = alpha * x,
/// ```
///
/// where `alpha` is a numeric, and `x` is a vector with `n` elements.
///
/// ## Signature
/// ```zig
/// linalg.blas.scal(n: usize, alpha: Al, x: [*]X, incx: isize) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn scal(n: usize, alpha: anytype, x: anytype, incx: isize) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.scal: alpha must be a numeric, and x must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    X = meta.Child(X);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return;

    return k_scal(n, alpha, x, incx);
}

/// Computes the product of a vector by a numeric:
///
/// ```zig
/// x = alpha * x,
/// ```
///
/// where `alpha` is a numeric, and `x` is a vector with `n` elements, splitting
/// the work across the worker threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.scalParallel(n: usize, alpha: Al, x: [*]X, incx: isize, pool: *thread.Pool) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vectors `x` and `y`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn scalParallel(
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    pool: *thread.Pool,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.scalParallel: alpha must be a numeric, and x must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    X = meta.Child(X);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return;

    const Ctx = struct {
        n: usize,
        alpha: Al,
        x: [*]X,
        incx: isize,

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            _ = worker_id;

            k_scal(
                end - start,
                ctx.alpha,
                ctx.x + numeric.cast(usize, if (ctx.incx > 0)
                    numeric.cast(isize, start) * ctx.incx
                else
                    (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incx),
                ctx.incx,
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
        },
        Ctx.kernel,
    );
}

fn k_scal(n: usize, alpha: anytype, x: anytype, incx: isize) void {
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
                    numeric.mulInto(
                        &x[i + u],
                        alpha,
                        x[i + u],
                    );
                }
            }

            while (i < n) : (i += 1) {
                // x[i] *= alpha
                numeric.mulInto(
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
                    numeric.mulInto(
                        &x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                        alpha,
                        x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                    );
                }

                ix += numeric.cast(isize, unroll) * incx;
            }

            while (i < n) : (i += 1) {
                // x[ix] *= alpha
                numeric.mulInto(
                    &x[numeric.cast(usize, ix)],
                    alpha,
                    x[numeric.cast(usize, ix)],
                );

                ix += incx;
            }
        }
    }
}
