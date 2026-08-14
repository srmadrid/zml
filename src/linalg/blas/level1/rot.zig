const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Applies a Givens plane rotation to the vectors `x` and `y`:
///
/// ```zig
/// xᵢ = c * xᵢ + s * yᵢ,
/// yᵢ = c * yᵢ - s * xᵢ,
/// ```
///
/// where `c` is a real numeric, but `s` may be either a real or a complex
/// numeric.
///
/// ## Signature
/// ```zig
/// linalg.blas.rot(n: usize, x: [*]X, incx: isize, y: [*]Y, incy: isize, c: C, s: S) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Number of elements in `x` and `y`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `c` (`anytype`): Cosine of the rotation angle. Real numeric.
/// * `s` (`anytype`): Sine of the rotation angle. Real or complex numeric.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn rot(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    c: anytype,
    s: anytype,
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

    if (incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (n == 0)
        return;

    return k_rot(n, x, incx, y, incy, c, s);
}

/// Applies a Givens plane rotation to the vectors `x` and `y`:
///
/// ```zig
/// xᵢ = c * xᵢ + s * yᵢ,
/// yᵢ = c * yᵢ - s * xᵢ,
/// ```
///
/// where `c` is a real numeric, but `s` may be either a real or a complex
/// numeric, splitting the work across the worker threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.rotParallel(n: usize, x: [*]X, incx: isize, y: [*]Y, incy: isize, c: C, s: S, pool: *thread.Pool) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Number of elements in `x` and `y`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `c` (`anytype`): Cosine of the rotation angle. Real numeric.
/// * `s` (`anytype`): Sine of the rotation angle. Real or complex numeric.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is equal to 0.
pub fn rotParallel(
    n: usize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    c: anytype,
    s: anytype,
    pool: *thread.Pool,
) !void {
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    const C: type = @TypeOf(c);
    const S: type = @TypeOf(s);

    comptime if (!meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)) or
        !meta.isNumeric(C) or !meta.isReal(C) or !meta.isNumeric(S))
        @compileError("zsl.linalg.blas.rotParallel: x and y must be mutable many-item pointers to numerics, c must be a real numeric, and s must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tc: " ++ @typeName(C) ++ "\n\ts: " ++ @typeName(S) ++ "\n");

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
        c: C,
        s: S,

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            _ = worker_id;

            k_rot(
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
                ctx.c,
                ctx.s,
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
            .c = c,
            .s = s,
        },
        Ctx.kernel,
    );
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
                numeric.fmaInto(
                    &x[i + u],
                    c,
                    xi,
                    numeric.mul(s, y[i + u]),
                );

                // y[i + u] = c * y[i + u] - conj(s) * xi
                numeric.fmaInto(
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
            numeric.fmaInto(
                &x[i],
                c,
                xi,
                numeric.mul(s, y[i]),
            );

            // y[i] = c * y[i] - conj(s) * xi
            numeric.fmaInto(
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
                numeric.fmaInto(
                    &x[x_idx],
                    c,
                    xi,
                    numeric.mul(s, y[y_idx]),
                );

                // y[y_idx] = c * y[y_idx] - conj(s) * xi
                numeric.fmaInto(
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
            numeric.fmaInto(
                &x[x_idx],
                c,
                xi,
                numeric.mul(s, y[y_idx]),
            );

            // y[y_idx] = c * y[y_idx] - conj(s) * xi
            numeric.fmaInto(
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
