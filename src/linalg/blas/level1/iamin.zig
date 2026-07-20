const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Finds the (0-based) index of the element with the smallest magnitude:
///
/// ```zig
/// argmin_i abs1(x[i]),   for i in 0..n
/// ```
///
/// If multiple elements share the minimum value, the smallest index is
/// returned.
///
/// ## Signature
/// ```zig
/// linalg.blas.iamin(n: usize, x: [*]const X, incx: isize) !usize
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Number of elements in `x`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `usize`: The 0-based index of the first element with the smallest magnitude.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn iamin(n: usize, x: anytype, incx: isize) !usize {
    const X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.iamin: x must be a many-item pointer to numerics, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return 0;

    return k_iamin(n, x, incx).index;
}

/// Finds the (0-based) index of the element with the smallest magnitude:
///
/// ```zig
/// argmin_i abs1(x[i]),   for i in 0..n
/// ```
///
/// splitting the work across the worker threads of `pool`. If multiple elements
/// share the minimum value, the smallest index is returned.
///
/// ## Signature
/// ```zig
/// linalg.blas.iaminParallel(n: usize, x: [*]const X, incx: isize, pool: *thread.Pool) !usize
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
/// `usize`: The 0-based index of the first element with the smallest magnitude.
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is equal to 0.
pub fn iaminParallel(n: usize, x: anytype, incx: isize, pool: *thread.Pool) !usize {
    const X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.iaminParallel: x must be a many-item pointer to numerics, got \n\tx: " ++ @typeName(X) ++ "\n");

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return 0;

    const Ctx = struct {
        n: usize,
        x: X,
        incx: isize,
        mins: *[thread.max_workers]?IaminResult(numeric.Abs1(meta.Child(X))),

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            var min = k_iamin(
                end - start,
                ctx.x + numeric.cast(usize, if (ctx.incx > 0)
                    numeric.cast(isize, start) * ctx.incx
                else
                    (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incx),
                ctx.incx,
            );

            min.index += start;

            if (ctx.mins[worker_id]) |current_min| {
                if (numeric.lt(min.value, current_min.value) or
                    (numeric.eq(min.value, current_min.value) and min.index < current_min.index))
                {
                    ctx.mins[worker_id] = .{
                        .value = min.value,
                        .index = min.index,
                    };
                }
            } else {
                ctx.mins[worker_id] = .{
                    .value = min.value,
                    .index = min.index,
                };
            }
        }
    };

    var mins: [thread.max_workers]?IaminResult(numeric.Abs1(meta.Child(X))) = @splat(null);

    pool.parallelFor(
        n,
        Ctx{
            .n = n,
            .x = x,
            .incx = incx,
            .mins = &mins,
        },
        Ctx.kernel,
    );

    var min: ?IaminResult(numeric.Abs1(meta.Child(X))) = null;
    for (0..int.min(pool.idCount(), thread.max_workers)) |i| {
        if (mins[i]) |min_i| {
            if (min != null) {
                if (numeric.lt(min_i.value, min.?.value) or
                    (min_i.value == min.?.value and min_i.index < min.?.index))
                {
                    min = .{
                        .value = min_i.value,
                        .index = min_i.index,
                    };
                }
            } else {
                min = .{
                    .value = min_i.value,
                    .index = min_i.index,
                };
            }
        }
    }

    return min.?.index;
}

fn IaminResult(N: type) type {
    return struct {
        value: N,
        index: usize,
    };
}

fn k_iamin(n: usize, x: anytype, incx: isize) IaminResult(numeric.Abs1(meta.Child(@TypeOf(x)))) {
    if (n == 0)
        return .{ .value = numeric.zero(numeric.Abs1(meta.Child(@TypeOf(x)))), .index = 0 };

    var best_value = if (incx == 1)
        numeric.abs1(x[0])
    else
        numeric.abs1(x[numeric.cast(usize, if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0)]);
    var best_index: usize = 0;

    if (incx == 1) {
        var i: usize = 1;
        while (i < n) : (i += 1) {
            const temp = numeric.abs1(x[i]);
            if (numeric.lt(temp, best_value)) {
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
            if (numeric.lt(temp, best_value)) {
                best_value = temp;
                best_index = i;
            }
            ix += incx;
        }
    }

    return .{ .value = best_value, .index = best_index };
}
