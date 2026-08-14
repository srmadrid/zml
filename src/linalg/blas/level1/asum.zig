const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

pub fn Asum(X: type) type {
    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.Asum: X must be a many-item pointer type to numerics, got \n\tX = " ++ @typeName(X) ++ "\n");

    return numeric.Abs1(meta.Child(X));
}

/// Computes the sum of magnitudes of the elements of a real vector, or the sum
/// of magnitudes of the real and imaginary parts of elements of a complex
/// vector:
///
/// ```zig
/// ∑ᵢ ‖xᵢ‖₁,
/// ```
///
/// where `x` is a vector with `n` elements.
///
/// ## Signature
/// ```zig
/// linalg.blas.asum(n: usize, x: [*]const X, incx: isize) !linalg.blas.Asum([*]const X)
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vector `x`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `Asum(@TypeOf(x))`: The sum of magnitudes of real and imaginary parts of all
/// elements of the vector.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn asum(n: usize, x: anytype, incx: isize) !linalg.blas.Asum(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return numeric.cast(linalg.blas.Asum(X), 0);

    return numeric.cast(linalg.blas.Asum(X), k_asum(n, x, incx));
}

/// Computes the sum of magnitudes of the elements of a real vector, or the sum
/// of magnitudes of the real and imaginary parts of elements of a complex
/// vector:
///
/// ```zig
/// ∑ᵢ ‖xᵢ‖₁,
/// ```
///
/// where `x` is a vector with `n` elements, splitting the work across the
/// worker threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.asumParallel(n: usize, x: [*]const X, incx: isize, pool: *thread.Pool) !linalg.blas.Asum([*]const X)
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the number of elements in vector `x`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `Asum(@TypeOf(x))`: The sum of magnitudes of real and imaginary parts of all
/// elements of the vector.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn asumParallel(n: usize, x: anytype, incx: isize, pool: *thread.Pool) !linalg.blas.Asum(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return numeric.cast(linalg.blas.Asum(X), 0);

    const Ctx = struct {
        n: usize,
        x: X,
        incx: isize,
        sums: *[thread.max_workers]meta.Accumulator(linalg.blas.Asum(X)),

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            numeric.addInto(
                &ctx.sums[worker_id],
                ctx.sums[worker_id],
                k_asum(
                    end - start,
                    ctx.x + numeric.cast(usize, if (ctx.incx > 0)
                        numeric.cast(isize, start) * ctx.incx
                    else
                        (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incx),
                    ctx.incx,
                ),
            );
        }
    };

    var sums: [thread.max_workers]meta.Accumulator(linalg.blas.Asum(X)) = @splat(numeric.cast(meta.Accumulator(linalg.blas.Asum(X)), 0));

    pool.parallelFor(
        n,
        Ctx{
            .n = n,
            .x = x,
            .incx = incx,
            .sums = &sums,
        },
        Ctx.kernel,
    );

    var sum = numeric.cast(meta.Accumulator(linalg.blas.Asum(X)), 0);
    for (0..int.min(pool.idCount(), thread.max_workers)) |i| {
        numeric.addInto(&sum, sum, sums[i]);
    }

    return numeric.cast(linalg.blas.Asum(X), sum);
}

fn k_asum(n: usize, x: anytype, incx: isize) meta.Accumulator(linalg.blas.Asum(@TypeOf(x))) {
    const X: type = @TypeOf(x);

    const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(linalg.blas.Asum(X))) orelse 2);

    var sums: [unroll]meta.Accumulator(linalg.blas.Asum(X)) = @splat(numeric.cast(meta.Accumulator(linalg.blas.Asum(X)), 0));

    if (incx == 1) {
        var i: usize = 0;
        while (i < (n / unroll) * unroll) : (i += unroll) {
            inline for (0..unroll) |u| {
                // sums[u] += abs1(x[i + u])
                numeric.addInto(
                    &sums[u],
                    sums[u],
                    numeric.abs1(x[i + u]),
                );
            }
        }

        while (i < n) : (i += 1) {
            // sums[0] += abs1(x[i])
            numeric.addInto(
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
                numeric.addInto(
                    &sums[u],
                    sums[u],
                    numeric.abs1(x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)]),
                );
            }

            ix += numeric.cast(isize, unroll) * incx;
        }

        while (i < n) : (i += 1) {
            // sums[0] += abs1(x[ix])
            numeric.addInto(
                &sums[0],
                sums[0],
                numeric.abs1(x[numeric.cast(usize, ix)]),
            );

            ix += incx;
        }
    }

    var sum = numeric.cast(meta.Accumulator(linalg.blas.Asum(X)), 0);
    inline for (0..unroll) |u| {
        // sum += sums[u]
        numeric.addInto(&sum, sum, sums[u]);
    }

    return sum;
}
