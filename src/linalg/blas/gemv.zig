const std = @import("std");
const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

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
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.gemv(layout: Layout, transa: linalg.Transpose, m: usize, n: usize, alpha: Al, a: [*]const A, lda: usize, x: [*]const X, incx: isize, beta: Be, y: [*]Y, incy: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `transa` (`linalg.Transpose`): Specifies the operation to be performed on
///   `A`:
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
/// * `opts`: Optional parameters:
///   * `num_threads` (`usize = 0`): Number of threads to spawn:
///     * `0`: automatic. The thread count is derived from `m * n` and
///       `parallel_threshold`:
///       ```zig
///       threads = max(1, min(std.Thread.getCpuCount(), options.max_threads, (m * n) / parallel_threshold))
///       ```
///     * 1: force serial execution. parallel_threshold is ignored.
///     * N >= 2: use exactly N threads, clamped by
///       std.Thread.getCpuCount() and options.max_threads as a hard safety
///       ceiling. parallel_threshold is ignored.
///   * parallel_threshold (usize = 8_388_608 / @sizeOf(meta.Child(Y))):
///     Minimum number of matrix elements (`m * n`) required to trigger
///     multithreaded execution.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `lda` is less than `max(1, m)` or
///   `max(1, n)`, or if `incx` or `incy` is 0.
pub fn gemv(
    layout: Layout,
    transa: linalg.Transpose,
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
    opts: struct {
        num_threads: usize = 0,
        parallel_threshold: usize = 8_388_608 / @sizeOf(meta.Child(@TypeOf(y))),
    },
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

    if (lda < int.max(1, if (layout == .col_major) m else n) or incx == 0 or incy == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Al == A and Al == X and Al == Y and Al == Be) {
        switch (comptime meta.numericType(Al)) {
            .float => {
                if (comptime Al == f32)
                    return linalg.cblas.sgemv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), alpha, a, numeric.cast(isize, lda), x, incx, beta, y, incy)
                else if (comptime Al == f64)
                    return linalg.cblas.dgemv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), alpha, a, numeric.cast(isize, lda), x, incx, beta, y, incy);
            },
            .complex => {
                if (comptime meta.Scalar(Al) == f32)
                    return linalg.cblas.cgemv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), &alpha, a, numeric.cast(isize, lda), x, incx, &beta, y, incy)
                else if (comptime meta.Scalar(Al) == f64)
                    return linalg.cblas.zgemv(layout.toInt(c_int), transa.toInt(c_int), numeric.cast(isize, m), numeric.cast(isize, n), &alpha, a, numeric.cast(isize, lda), x, incx, &beta, y, incy);
            },
            else => {},
        }
    }

    const eff_transa = if (layout == .col_major) transa else transa.invert();
    const eff_m = if (layout == .col_major) m else n;
    const eff_n = if (layout == .col_major) n else m;
    const noconj = transa == .no_trans or transa == .trans;

    if (opts.num_threads == 1)
        return if (noconj)
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, false);

    var num_threads: usize = if (opts.num_threads == 0) blk: {
        if (opts.parallel_threshold == 0)
            break :blk options.max_threads;

        break :blk int.max(1, (eff_m * eff_n) / opts.parallel_threshold);
    } else opts.num_threads;

    num_threads = int.min(num_threads, options.max_threads);

    if (num_threads <= 1)
        return if (noconj)
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, false);

    num_threads = int.min(num_threads, std.Thread.getCpuCount() catch 1);

    if (num_threads <= 1)
        return if (noconj)
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, true)
        else
            k_gemv(eff_transa, eff_m, eff_n, alpha, a, lda, x, incx, beta, y, incy, false);

    var threads: [options.max_threads]std.Thread = undefined;

    const Worker = struct {
        fn execute(worker_transa: linalg.Transpose, worker_m: usize, worker_n: usize, worker_alpha: Al, worker_a: [*]const A, worker_lda: usize, worker_x: [*]const X, worker_incx: isize, worker_beta: Be, worker_y: [*]Y, worker_incy: isize, worker_noconj: bool) void {
            return if (worker_noconj)
                k_gemv(worker_transa, worker_m, worker_n, worker_alpha, worker_a, worker_lda, worker_x, worker_incx, worker_beta, worker_y, worker_incy, true)
            else
                k_gemv(worker_transa, worker_m, worker_n, worker_alpha, worker_a, worker_lda, worker_x, worker_incx, worker_beta, worker_y, worker_incy, false);
        }
    };

    const no_trans = eff_transa == .no_trans or eff_transa == .conj_no_trans;
    const leny = if (no_trans) eff_m else eff_n;
    const chunk_size = int.div(leny, num_threads);
    var spawn_err: ?anyerror = null;
    var spawned_count: usize = 0;
    var i: usize = 0;
    while (i < num_threads) : (i += 1) {
        const chunk_start = i * chunk_size;
        const chunk_end = if (i == num_threads - 1) leny else chunk_start + chunk_size;
        const chunk_len = chunk_end - chunk_start;

        if (std.Thread.spawn(.{}, Worker.execute, .{
            eff_transa,
            if (no_trans) chunk_len else eff_m,
            if (no_trans) eff_n else chunk_len,
            alpha,
            a + if (no_trans) chunk_start else (chunk_start * lda),
            lda,
            x,
            incx,
            beta,
            y + numeric.cast(usize, if (incy > 0)
                numeric.cast(isize, chunk_start) * incy
            else
                (-numeric.cast(isize, leny) + numeric.cast(isize, chunk_end)) * incy),
            incy,
            noconj,
        })) |th| {
            threads[i] = th;
            spawned_count += 1;
        } else |err| {
            spawn_err = err;
            break;
        }
    }

    var t: usize = 0;
    while (t < spawned_count) : (t += 1) {
        threads[t].join();
    }

    if (spawn_err) |err|
        return err;
}

fn k_gemv(transa: linalg.Transpose, m: usize, n: usize, alpha: anytype, a: anytype, lda: usize, x: anytype, incx: isize, beta: anytype, y: anytype, incy: isize, comptime noconj: bool) void {
    const Al: type = @TypeOf(alpha);
    const A: type = meta.Child(@TypeOf(a));
    const X: type = meta.Child(@TypeOf(x));
    const Y: type = meta.Child(@TypeOf(y));

    // Quick return if possible.
    if (m == 0 or n == 0 or (numeric.eq(alpha, 0) and numeric.eq(beta, 1)))
        return;

    // Set lenx and leny, the lengths of the vectors x and y, and set up the
    // start points in x and y.
    const lenx: usize = if (transa == .no_trans or transa == .conj_no_trans) n else m;
    const leny: usize = if (transa == .no_trans or transa == .conj_no_trans) m else n;

    const kx: isize = if (incx < 0) (-numeric.cast(isize, lenx) + 1) * incx else 0;
    const ky: isize = if (incy < 0) (-numeric.cast(isize, leny) + 1) * incy else 0;

    // First form  y = beta * y.
    if (numeric.ne(beta, 1))
        @import("scal.zig").k_scal(leny, beta, y, incy);

    if (numeric.eq(alpha, 0))
        return;

    if (transa == .no_trans or transa == .conj_no_trans) {
        const unroll = 2 * (std.simd.suggestVectorLength(numeric.Fma(numeric.Mul(Al, X), A, Y)) orelse 2);

        // Form  y = alpha * A * x + y  or  y = alpha * conj(A) * x + y.
        var jx: isize = kx;
        var j: usize = 0;
        while (j < n) : (j += 1) {
            // temp = alpha * x[jx]
            const temp = numeric.mul(
                alpha,
                x[numeric.cast(usize, jx)],
            );

            var i: usize = 0;
            if (incy == 1) {
                while (i < (m / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // y[i + u] += temp * a[i + u + j * lda]
                        numeric.fma_(
                            &y[i + u],
                            temp,
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            y[i + u],
                        );
                    }
                }

                while (i < m) : (i += 1) {
                    // y[i] += temp * a[i + j * lda]
                    numeric.fma_(
                        &y[i],
                        temp,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[i],
                    );
                }
            } else {
                var iy: isize = ky;
                while (i < (m / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // y[iy + u * incy] += temp * a[i + u + j * lda]
                        numeric.fma_(
                            &y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                            temp,
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            y[numeric.cast(usize, iy + numeric.cast(isize, u) * incy)],
                        );
                    }

                    iy += numeric.cast(isize, unroll) * incy;
                }

                while (i < m) : (i += 1) {
                    // y[iy] += temp * a[i + j * lda]
                    numeric.fma_(
                        &y[numeric.cast(usize, iy)],
                        temp,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        y[numeric.cast(usize, iy)],
                    );

                    iy += incy;
                }
            }

            jx += incx;
        }
    } else {
        const unroll = 2 * (std.simd.suggestVectorLength(meta.Accumulator(numeric.Mul(A, X))) orelse 2);

        // Form  y = alpha * Aᵀ * x + y  or  y = alpha * Aᴴ * x + y.
        var jy: isize = ky;
        var j: usize = 0;
        while (j < n) : (j += 1) {
            var temp = numeric.zero(meta.Accumulator(numeric.Mul(A, X)));

            var i: usize = 0;
            if (incx == 1) {
                var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                while (i < (m / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // sums[u] += a[i + u + j * lda] * x[i + u]
                        numeric.fma_(
                            &sums[u],
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            x[i + u],
                            sums[u],
                        );
                    }
                }

                inline for (0..unroll) |u| {
                    numeric.add_(&temp, temp, sums[u]);
                }

                while (i < m) : (i += 1) {
                    // temp += a[i + j * lda] * x[i]
                    numeric.fma_(
                        &temp,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[i],
                        temp,
                    );
                }
            } else {
                var sums: [unroll]meta.Accumulator(numeric.Mul(A, X)) = .{numeric.zero(meta.Accumulator(numeric.Mul(A, X)))} ** unroll;

                var ix: isize = kx;
                while (i < (m / unroll) * unroll) : (i += unroll) {
                    inline for (0..unroll) |u| {
                        // sums[u] += a[i + u + j * lda] * x[ix + u * incx]
                        numeric.fma_(
                            &sums[u],
                            if (comptime noconj)
                                a[i + u + j * lda]
                            else
                                numeric.conj(a[i + u + j * lda]),
                            x[numeric.cast(usize, ix + numeric.cast(isize, u) * incx)],
                            sums[u],
                        );
                    }

                    ix += numeric.cast(isize, unroll) * incx;
                }

                inline for (0..unroll) |u| {
                    numeric.add_(&temp, temp, sums[u]);
                }

                while (i < m) : (i += 1) {
                    // temp += a[i + j * lda] * x[ix]
                    numeric.fma_(
                        &temp,
                        if (comptime noconj)
                            a[i + j * lda]
                        else
                            numeric.conj(a[i + j * lda]),
                        x[numeric.cast(usize, ix)],
                        temp,
                    );

                    ix += incx;
                }
            }

            // y[jy] += alpha * temp
            numeric.fma_(
                &y[numeric.cast(usize, jy)],
                alpha,
                temp,
                y[numeric.cast(usize, jy)],
            );

            jy += incy;
        }
    }

    return;
}
