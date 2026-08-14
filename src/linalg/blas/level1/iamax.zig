const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Finds the (0-based) index of the element with the largest magnitude:
///
/// ```zig
/// argmaxᵢ ‖xᵢ‖₁.
/// ```
///
/// If multiple elements share the maximum value, the smallest index is
/// returned.
///
/// ## Signature
/// ```zig
/// linalg.blas.iamax(n: usize, x: [*]const X, incx: isize) !usize
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Number of elements in `x`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `usize`: The 0-based index of the first element with the largest magnitude.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn iamax(n: usize, x: anytype, incx: isize) !usize {
    const X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.iamax: x must be a many-item pointer to numerics, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return 0;

    return k_iamax(n, x, incx).index;
}

/// Finds the (0-based) index of the element with the largest magnitude:
///
/// ```zig
/// argmaxᵢ ‖xᵢ‖₁,
/// ```
///
/// splitting the work across the worker threads of `pool`. If multiple elements
/// share the maximum value, the smallest index is returned.
///
/// ## Signature
/// ```zig
/// linalg.blas.iamaxParallel(n: usize, x: [*]const X, incx: isize, pool: *thread.Pool) !usize
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Number of elements in `x`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `usize`: The 0-based index of the first element with the largest magnitude.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn iamaxParallel(n: usize, x: anytype, incx: isize, pool: *thread.Pool) !usize {
    const X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.iamaxParallel: x must be a many-item pointer to numerics, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return 0;

    const Ctx = struct {
        n: usize,
        x: X,
        incx: isize,
        maxes: *[thread.max_workers]?IamaxResult(numeric.Abs1(meta.Child(X))),

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            var max = k_iamax(
                end - start,
                ctx.x + numeric.cast(usize, if (ctx.incx > 0)
                    numeric.cast(isize, start) * ctx.incx
                else
                    (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incx),
                ctx.incx,
            );

            max.index += start;

            if (ctx.maxes[worker_id]) |current_max| {
                if (numeric.gt(max.value, current_max.value) or
                    (numeric.eq(max.value, current_max.value) and max.index < current_max.index))
                {
                    ctx.maxes[worker_id] = .{
                        .value = max.value,
                        .index = max.index,
                    };
                }
            } else {
                ctx.maxes[worker_id] = .{
                    .value = max.value,
                    .index = max.index,
                };
            }
        }
    };

    var maxes: [thread.max_workers]?IamaxResult(numeric.Abs1(meta.Child(X))) = @splat(null);

    pool.parallelFor(
        n,
        Ctx{
            .n = n,
            .x = x,
            .incx = incx,
            .maxes = &maxes,
        },
        Ctx.kernel,
    );

    var max: ?IamaxResult(numeric.Abs1(meta.Child(X))) = null;
    for (0..int.min(pool.idCount(), thread.max_workers)) |i| {
        if (maxes[i]) |max_i| {
            if (max != null) {
                if (numeric.gt(max_i.value, max.?.value) or
                    (max_i.value == max.?.value and max_i.index < max.?.index))
                {
                    max = .{
                        .value = max_i.value,
                        .index = max_i.index,
                    };
                }
            } else {
                max = .{
                    .value = max_i.value,
                    .index = max_i.index,
                };
            }
        }
    }

    return max.?.index;
}

fn IamaxResult(N: type) type {
    return struct {
        value: N,
        index: usize,
    };
}

fn k_iamax(n: usize, x: anytype, incx: isize) IamaxResult(numeric.Abs1(meta.Child(@TypeOf(x)))) {
    if (n == 0)
        return .{ .value = numeric.cast(numeric.Abs1(meta.Child(@TypeOf(x))), 0), .index = 0 };

    var best_value = if (incx == 1)
        numeric.abs1(x[0])
    else
        numeric.abs1(x[numeric.cast(usize, if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0)]);
    var best_index: usize = 0;

    if (incx == 1) {
        var i: usize = 1;
        while (i < n) : (i += 1) {
            const temp = numeric.abs1(x[i]);
            if (numeric.gt(temp, best_value)) {
                best_value = temp;
                best_index = i;
            }
        }
    } else {
        var ix: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

        ix += incx;

        var i: usize = 1;
        while (i < n) : (i += 1) {
            const temp = numeric.abs1(x[numeric.cast(usize, ix)]);
            if (numeric.gt(temp, best_value)) {
                best_value = temp;
                best_index = i;
            }
            ix += incx;
        }
    }

    return .{ .value = best_value, .index = best_index };
}
