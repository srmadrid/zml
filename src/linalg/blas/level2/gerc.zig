const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Performs a rank-1 update (conjugated) of a general matrix defined as:
///
/// ```zig
/// A = alpha * x * yᴴ + A,
/// ```
///
/// where `alpha` is a numeric, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.gerc(layout: matrix.Layout, m: usize, n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, a: [*]A, lda: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (m - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * k`, where
///   `k` is `n` when `layout` is `col_major`, or `m` when `layout` is
///   `row_major`. On return, contains the result of the operation.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0, or if `lda`
///   is less than `max(1, m)` or `max(1, n)`.
pub fn gerc(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: usize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or !meta.isNumeric(meta.Child(Y)) or
        !meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.blas.gerc: alpha must be a numeric, x and y must be many-item pointers to numerics, and a must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\ta: " ++ @typeName(A) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);
    A = meta.Child(A);

    if (lda < int.max(1, if (layout == .col_major) m else n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0 or numeric.eq(alpha, 0))
        return;

    if (layout == .col_major)
        k_gerc(m, n, alpha, x, incx, y, incy, a, lda, true)
    else
        k_gerc(n, m, alpha, y, incy, x, incx, a, lda, false);
}

/// Performs a rank-1 update (conjugated) of a general matrix defined as:
///
/// ```zig
/// A = alpha * x * yᴴ + A,
/// ```
///
/// where `alpha` is a numeric, `x` is an `m`-element vector, `y` is an
/// `n`-element vector, and `A` is an `m`-by-`n` general matrix, splitting the
/// work across the worker threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.gercParallel(layout: matrix.Layout, m: usize, n: usize, alpha: Al, x: [*]const X, incx: isize, y: [*]const Y, incy: isize, a: [*]A, lda: usize, pool: *thread.Pool) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (m - 1) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `y` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incy)`.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * k`, where
///   `k` is `n` when `layout` is `col_major`, or `m` when `layout` is
///   `row_major`. On return, contains the result of the operation.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` or `incy` is 0, or if `lda`
///   is less than `max(1, m)` or `max(1, n)`.
pub fn gercParallel(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    a: anytype,
    lda: usize,
    pool: *thread.Pool,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var X: type = @TypeOf(x);
    comptime var Y: type = @TypeOf(y);
    comptime var A: type = @TypeOf(a);

    comptime if (!meta.isNumeric(Al) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or !meta.isNumeric(meta.Child(Y)) or
        !meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.blas.gercParallel: alpha must be a numeric, x and y must be many-item pointers to numerics, and a must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\ta: " ++ @typeName(A) ++ "\n");

    X = meta.Child(X);
    Y = meta.Child(Y);
    A = meta.Child(A);

    if (lda < int.max(1, if (layout == .col_major) m else n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    const eff_m = if (layout == .col_major) m else n;
    const eff_n = if (layout == .col_major) n else m;

    // Quick return if possible.
    if (m == 0 or n == 0 or numeric.eq(alpha, 0))
        return;

    const Ctx = struct {
        m: usize,
        n: usize,
        alpha: Al,
        x: [*]const X,
        incx: isize,
        y: [*]const Y,
        incy: isize,
        a: [*]A,
        lda: usize,
        iscolm: bool,

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            _ = worker_id;

            const chunk_len = end - start;
            const worker_m = ctx.m;
            const worker_n = chunk_len;
            const worker_x = if (ctx.iscolm)
                ctx.x
            else
                ctx.x + numeric.cast(usize, if (ctx.incx > 0)
                    numeric.cast(isize, start) * ctx.incx
                else
                    (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incx);
            const worker_y = if (ctx.iscolm)
                ctx.y + numeric.cast(usize, if (ctx.incy > 0)
                    numeric.cast(isize, start) * ctx.incy
                else
                    (-numeric.cast(isize, ctx.n) + numeric.cast(isize, end)) * ctx.incy)
            else
                ctx.y;
            const worker_a = ctx.a + start * ctx.lda;

            if (ctx.iscolm)
                k_gerc(worker_m, worker_n, ctx.alpha, worker_x, ctx.incx, worker_y, ctx.incy, worker_a, ctx.lda, true)
            else
                k_gerc(worker_m, worker_n, ctx.alpha, worker_y, ctx.incy, worker_x, ctx.incx, worker_a, ctx.lda, false);
        }
    };

    pool.parallelFor(
        eff_n,
        Ctx{
            .m = eff_m,
            .n = eff_n,
            .alpha = alpha,
            .x = x,
            .incx = incx,
            .y = y,
            .incy = incy,
            .a = a,
            .lda = lda,
            .iscolm = layout == .col_major,
        },
        Ctx.kernel,
    );
}

fn k_gerc(m: usize, n: usize, alpha: anytype, x: anytype, incx: isize, y: anytype, incy: isize, a: anytype, lda: usize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Form  A = alpha * x * yᴴ + A
    const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(if (comptime noconj) X else numeric.Conj(X), numeric.Mul(Al, if (comptime noconj) numeric.Conj(Y) else Y), A)) orelse 2);
    comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
    tile_size = comptime int.max(1, tile_size -| tile_size % unroll);

    var tile_i: usize = 0;
    while (tile_i < m) : (tile_i += tile_size) {
        const b_len = int.min(tile_size, m - tile_i);
        var local_x: [tile_size]X = undefined;

        const px = if (incx == 1)
            x + tile_i
        else blk: {
            linalg.blas.copy(
                b_len,
                x + numeric.cast(usize, if (incx > 0)
                    numeric.cast(isize, tile_i) * incx
                else
                    (-numeric.cast(isize, m) + numeric.cast(isize, tile_i + b_len)) * incx),
                incx,
                @as([*]X, &local_x),
                1,
            ) catch unreachable;

            break :blk @as([*]const X, &local_x);
        };

        var jy: isize = if (incy < 0) (-numeric.cast(isize, n) + 1) * incy else 0;
        var j: usize = 0;
        while (j < n) : (j += 1) {
            if (numeric.ne(y[numeric.cast(usize, jy)], 0)) {
                // temp = alpha * conj(y[jy])
                const temp = numeric.mul(
                    alpha,
                    if (comptime noconj)
                        numeric.conj(y[numeric.cast(usize, jy)])
                    else
                        y[numeric.cast(usize, jy)],
                );

                var i: usize = 0;
                while (i < (b_len / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // a[tile_i + i + u + j * lda] += px[i + u] * temp
                        numeric.fmaInto(
                            &a[tile_i + i + u + j * lda],
                            if (comptime noconj)
                                px[i + u]
                            else
                                numeric.conj(px[i + u]),
                            temp,
                            a[tile_i + i + u + j * lda],
                        );
                    }
                }

                while (i < b_len) : (i += 1) {
                    // a[tile_i + i + j * lda] += px[i] * temp
                    numeric.fmaInto(
                        &a[tile_i + i + j * lda],
                        if (comptime noconj)
                            px[i]
                        else
                            numeric.conj(px[i]),
                        temp,
                        a[tile_i + i + j * lda],
                    );
                }
            }

            jy += incy;
        }
    }

    return;
}
