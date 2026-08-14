const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

pub fn Nrm2(X: type) type {
    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.Nrm2: X must be a many-item pointer type to numerics, got \n\tX = " ++ @typeName(X) ++ "\n");

    return numeric.Sqrt(numeric.Abs2(meta.Child(X)));
}

/// Computes the Euclidean norm of a vector:
///
/// ```zig
/// ‖x‖₂ = √(∑ᵢ |xᵢ|²),
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
/// * `n` (`usize`): Specifies the number of elements in vector `x`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `Nrm2(@TypeOf(x))`: The Euclidean norm of `x`.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn nrm2(n: usize, x: anytype, incx: isize) !linalg.blas.Nrm2(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return numeric.cast(linalg.blas.Nrm2(X), 0);

    return numeric.cast(linalg.blas.Nrm2(X), numeric.sqrt(k_nrm2_ssq(n, x, incx)));
}

/// Computes the Euclidean norm of a vector:
///
/// ```zig
/// ‖x‖₂ = √(∑ᵢ |xᵢ|²),
/// ```
///
/// where `x` is a vector with `n` elements, splitting the work across the
/// worker threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.nrm2Parallel(n: usize, x: [*]const X, incx: isize, pool: *thread.Pool) !linalg.blas.Nrm2([*]const X)
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
/// `Nrm2(@TypeOf(x))`: The Euclidean norm of `x`.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn nrm2Parallel(n: usize, x: anytype, incx: isize, pool: *thread.Pool) !linalg.blas.Nrm2(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return numeric.cast(linalg.blas.Nrm2(X), 0);

    const Ctx = struct {
        n: usize,
        x: X,
        incx: isize,
        sums: *[thread.max_workers]meta.Accumulator(linalg.blas.Asum(X)),

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            numeric.addInto(
                &ctx.sums[worker_id],
                ctx.sums[worker_id],
                k_nrm2_ssq(
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

    return numeric.cast(linalg.blas.Asum(X), numeric.sqrt(sum));
}

fn k_nrm2_ssq(n: usize, x: anytype, incx: isize) meta.Accumulator(linalg.blas.Nrm2(@TypeOf(x))) {
    const X: type = @TypeOf(x);

    const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(linalg.blas.Nrm2(X))) orelse 2);

    var sums: [unroll]meta.Accumulator(linalg.blas.Nrm2(X)) = @splat(numeric.cast(meta.Accumulator(linalg.blas.Nrm2(X)), 0));

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

    var ssq = numeric.cast(meta.Accumulator(linalg.blas.Nrm2(X)), 0);
    inline for (0..unroll) |u| {
        numeric.addInto(&ssq, ssq, sums[u]);
    }

    return ssq;
}
