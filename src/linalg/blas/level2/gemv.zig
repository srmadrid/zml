const std = @import("std");

const options = @import("options");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const thread = @import("../../../thread.zig");

/// Computes a matrix-vector product with a general matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᵀ * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᴴ * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are vectors, and `A` is an
/// `m`-by-`n` matrix.
///
/// ## Signature
/// ```zig
/// linalg.blas.gemv(layout: matrix.Layout, transa: linalg.blas.Transpose, m: usize, n: usize, alpha: Al, a: [*]const A, lda: usize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `transa` (`linalg.blas.Transpose`): Specifies the operation to be
///   performed on `A`:
///   * `no_transpose`: `y = alpha * A * x + beta * y`
///   * `transpose`: `y = alpha * Aᵀ * x + beta * y`
///   * `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
///   * `conj_transpose`: `y = alpha * Aᴴ * x + beta * y`
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * k`, where `k` is
///   `n` when `layout` is `col_major`, or `m` when `layout` is `row_major`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (m - 1) * abs(incx)` otherwise.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `y` need not be set on input.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (m - 1) * abs(incy)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (n - 1) * abs(incy)` otherwise. On return,
///   contains the result of the operation.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, m)` or
///   `max(1, n)`, or if `incx` or `incy` is 0.
pub fn gemv(
    layout: matrix.Layout,
    transa: linalg.blas.Transpose,
    m: usize,
    n: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.gemv: alpha and beta must be numerics, a and x must be many-item pointers to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);
    Y = meta.Child(Y);

    const eff_transa = if (layout == .col_major) transa else transa.invert();
    const eff_m = if (layout == .col_major) m else n;
    const eff_n = if (layout == .col_major) n else m;
    const notrans = eff_transa == .no_trans or eff_transa == .conj_no_trans;
    const noconj = transa == .no_trans or transa == .trans;
    const leny = if (notrans) eff_m else eff_n;

    if (lda < int.max(1, eff_m) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return;

    if (numeric.eq(alpha, 0)) {
        if (numeric.ne(beta, 1))
            linalg.blas.scal(leny, beta, y, incy) catch unreachable;

        return;
    }

    if (noconj)
        k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, true)
    else
        k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, false);
}

/// Computes a matrix-vector product with a general matrix defined as:
///
/// ```zig
/// y = alpha * A * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᵀ * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * conj(A) * x + beta * y,
/// ```
///
/// or
///
/// ```zig
/// y = alpha * Aᴴ * x + beta * y,
/// ```
///
/// where `alpha` and `beta` are numerics, `x` and `y` are vectors, and `A` is an
/// `m`-by-`n` matrix, splitting the work across the worker threads of `pool`.
///
/// ## Signature
/// ```zig
/// linalg.blas.gemvParallel(layout: matrix.Layout, transa: linalg.blas.Transpose, m: usize, n: usize, alpha: Al, a: [*]const A, lda: usize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize, pool: *thread.Pool) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `transa` (`linalg.blas.Transpose`): Specifies the operation to be
///   performed on `A`:
///   * `no_transpose`: `y = alpha * A * x + beta * y`
///   * `transpose`: `y = alpha * Aᵀ * x + beta * y`
///   * `conj_no_transpose`: `y = alpha * conj(A) * x + beta * y`
///   * `conj_transpose`: `y = alpha * Aᴴ * x + beta * y`
/// * `m` (`usize`): Specifies the number of rows of the matrix `A`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `alpha` (`anytype`): Specifies the numeric `alpha`.
/// * `a` (`anytype`): Many-item pointer, size at least `lda * k`, where `k` is
///   `n` when `layout` is `col_major`, or `m` when `layout` is `row_major`.
/// * `lda` (`usize`): Specifies the leading dimension of `a` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
/// * `x` (`anytype`): Many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (m - 1) * abs(incx)` otherwise.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `beta` (`anytype`): Specifies the numeric `beta`. When `beta` is 0, then
///   `y` need not be set on input.
/// * `y` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (m - 1) * abs(incy)` when `transa` is `no_transpose` or
///   `conj_no_transpose`, or `1 + (n - 1) * abs(incy)` otherwise. On return,
///   contains the result of the operation.
/// * `incy` (`isize`): Indexing increment for `y`. Must be different from 0.
/// * `pool` (`*thread.Pool`): Thread pool used to parallelize the computation.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, m)` or
///   `max(1, n)`, or if `incx` or `incy` is 0.
pub fn gemvParallel(
    layout: matrix.Layout,
    transa: linalg.blas.Transpose,
    m: usize,
    n: usize,
    alpha: anytype,
    a: anytype,
    lda: usize,
    x: anytype,
    incx: isize,
    beta: anytype,
    y: anytype,
    incy: isize,
    pool: *thread.Pool,
) !void {
    const Al: type = @TypeOf(alpha);
    comptime var A: type = @TypeOf(a);
    comptime var X: type = @TypeOf(x);
    const Be: type = @TypeOf(beta);
    comptime var Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(Al) or !meta.isNumeric(Be) or
        !meta.isManyItemPointer(A) or !meta.isNumeric(meta.Child(A)) or
        !meta.isManyItemPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isManyItemPointer(Y) or meta.isConstPointer(Y) or !meta.isNumeric(meta.Child(Y)))
        @compileError("zsl.linalg.blas.gemvParallel: alpha and beta must be numerics, a and x must be many-item pointers to numerics, and y must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\ta: " ++ @typeName(A) ++ "\n\tx: " ++ @typeName(X) ++ "\n\tbeta: " ++ @typeName(Be) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    A = meta.Child(A);
    X = meta.Child(X);
    Y = meta.Child(Y);

    const eff_transa = if (layout == .col_major) transa else transa.invert();
    const eff_m = if (layout == .col_major) m else n;
    const eff_n = if (layout == .col_major) n else m;
    const notrans = eff_transa == .no_trans or eff_transa == .conj_no_trans;
    const noconj = transa == .no_trans or transa == .trans;
    const leny = if (notrans) eff_m else eff_n;

    if (lda < int.max(1, eff_m) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    // Quick return if possible.
    if (m == 0 or n == 0)
        return;

    if (numeric.eq(alpha, 0)) {
        if (numeric.ne(beta, 1))
            linalg.blas.scal(leny, beta, y, incy) catch unreachable;

        return;
    }

    const Ctx = struct {
        transa: linalg.blas.Transpose,
        m: usize,
        n: usize,
        alpha: Al,
        a: [*]const A,
        lda: usize,
        x: [*]const X,
        incx: isize,
        beta: Be,
        y: [*]Y,
        incy: isize,
        notrans: bool,
        noconj: bool,
        leny: usize,

        fn kernel(ctx: @This(), start: usize, end: usize, worker_id: usize) void {
            _ = worker_id;

            const chunk_len = end - start;
            const worker_m = if (ctx.notrans) chunk_len else ctx.m;
            const worker_n = if (ctx.notrans) ctx.n else chunk_len;
            const worker_a = ctx.a + if (ctx.notrans) start else (start * ctx.lda);
            const worker_y = ctx.y + numeric.cast(usize, if (ctx.incy > 0)
                numeric.cast(isize, start) * ctx.incy
            else
                (-numeric.cast(isize, ctx.leny) + numeric.cast(isize, end)) * ctx.incy);

            if (ctx.noconj)
                k_gemv(ctx.transa, worker_m, worker_n, ctx.alpha, worker_a, ctx.lda, ctx.x, ctx.incx, ctx.beta, worker_y, ctx.incy, true)
            else
                k_gemv(ctx.transa, worker_m, worker_n, ctx.alpha, worker_a, ctx.lda, ctx.x, ctx.incx, ctx.beta, worker_y, ctx.incy, false);
        }
    };

    pool.parallelFor(
        leny,
        Ctx{
            .transa = eff_transa,
            .m = eff_m,
            .n = eff_n,
            .alpha = alpha,
            .a = a,
            .lda = lda,
            .x = x,
            .incx = incx,
            .beta = beta,
            .y = y,
            .incy = incy,
            .notrans = notrans,
            .noconj = noconj,
            .leny = leny,
        },
        Ctx.kernel,
    );
}

fn k_gemv(transa: linalg.blas.Transpose, m: usize, n: usize, alpha: anytype, a: anytype, lda: usize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Set lenx and leny, the lengths of the vectors x and y.
    const lenx: usize = if (transa == .no_trans or transa == .conj_no_trans) n else m;
    const leny: usize = if (transa == .no_trans or transa == .conj_no_trans) m else n;

    // First form  y = beta * y.
    if (numeric.ne(beta, 1))
        linalg.blas.scal(leny, beta, y, incy) catch unreachable;

    // Set up the start points in x and y.
    const kx: isize = if (incx < 0) (-numeric.cast(isize, lenx) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, leny) + 1) * incy else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        // Form  y = alpha * A * x + y  or  y = alpha * conj(A) * x + y.
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(Y) + @sizeOf(A)));
        tile_size = comptime int.max(1, tile_size -| (tile_size % unroll));

        var tile_i: usize = 0;
        while (tile_i < m) : (tile_i += tile_size) {
            const b_len = int.min(tile_size, m - tile_i);
            var local_y: [tile_size]Y = undefined;

            const py = if (incy == 1)
                y + tile_i
            else blk: {
                linalg.blas.copy(
                    b_len,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, m) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                    @as([*]Y, &local_y),
                    1,
                ) catch unreachable;

                break :blk @as([*]Y, &local_y);
            };

            var jx: isize = kx;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                // temp = alpha * x[jx]
                const temp = numeric.mul(
                    alpha,
                    x[numeric.cast(usize, jx)],
                );

                var i: usize = 0;
                while (i < (b_len / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // py[i + u] += temp * a[tile_i + i + u + j * lda]
                        numeric.fmaInto(
                            &py[i + u],
                            temp,
                            if (comptime noconj)
                                a[tile_i + i + u + j * lda]
                            else
                                numeric.conj(a[tile_i + i + u + j * lda]),
                            py[i + u],
                        );
                    }
                }

                while (i < b_len) : (i += 1) {
                    // py[i] += temp * a[tile_i + i + j * lda]
                    numeric.fmaInto(
                        &py[i],
                        temp,
                        if (comptime noconj)
                            a[tile_i + i + j * lda]
                        else
                            numeric.conj(a[tile_i + i + j * lda]),
                        py[i],
                    );
                }

                jx += incx;
            }

            if (incy != 1) {
                linalg.blas.copy(
                    b_len,
                    @as([*]Y, &local_y),
                    1,
                    y + numeric.cast(usize, if (incy > 0)
                        numeric.cast(isize, tile_i) * incy
                    else
                        (-numeric.cast(isize, m) + numeric.cast(isize, tile_i + b_len)) * incy),
                    incy,
                ) catch unreachable;
            }
        }
    } else {
        // Form  y = alpha * Aᵀ * x + y  or  y = alpha * Aᴴ * x + y.
        const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, X))) orelse 2);
        comptime var tile_size = int.max(1, ((3 * options.l1_size) / 4) / (@sizeOf(X) + @sizeOf(A)));
        tile_size -|= tile_size % unroll;

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

            var jy: isize = ky;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

                var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                var i: usize = 0;
                while (i < (b_len / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // sums[u] += a[tile_i + i + u + j * lda] * px[i + u]
                        numeric.fmaInto(
                            &sums[u],
                            if (comptime noconj)
                                a[tile_i + i + u + j * lda]
                            else
                                numeric.conj(a[tile_i + i + u + j * lda]),
                            px[i + u],
                            sums[u],
                        );
                    }
                }

                inline for (0..unroll) |u| {
                    numeric.addInto(&temp, temp, sums[u]);
                }

                while (i < b_len) : (i += 1) {
                    // temp += a[tile_i + i + j * lda] * x[i]
                    numeric.fmaInto(
                        &temp,
                        if (comptime noconj)
                            a[tile_i + i + j * lda]
                        else
                            numeric.conj(a[tile_i + i + j * lda]),
                        px[i],
                        temp,
                    );
                }

                // y[jy] += alpha * temp
                numeric.fmaInto(
                    &y[numeric.cast(usize, jy)],
                    alpha,
                    temp,
                    y[numeric.cast(usize, jy)],
                );

                jy += incy;
            }
        }
    }

    return;
}
