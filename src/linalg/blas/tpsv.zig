const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;
const Diag = meta.Diag;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Solves a system of linear equations whose coefficients are in a triangular
/// packed matrix defined as:
///
/// ```zig
/// A * x = b,
/// ```
///
/// or
///
/// ```zig
/// conj(A) * x = b,
/// ```
///
/// or
///
/// ```zig
/// Aᵀ * x = b,
/// ```
///
/// or
///
/// ```zig
/// Aᴴ * x = b,
/// ```
///
/// where `b` and `x` are `n`-element vectors, and `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix, supplied in packed form.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.tpsv(layout: Layout, uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, ap: [*]const Ap, x: [*]X, incx: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
///   triangular matrix.
/// * `transa` (`linalg.Transpose`): Specifies the system of equations to be
///   solved:
///   * `no_transpose`: `A * x = b`
///   * `transpose`: `Aᵀ * x = b`
///   * `conj_no_transpose`: `conj(A) * x = b`
///   * `conj_transpose`: `Aᴴ * x = b`
/// * `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `ap` (`anytype`): Many-item pointer, size at least `(n * (n + 1)) / 2`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. Before entry, the incremented array `x` must
///   contain the n-element right-hand side vector `b`. On return, it contains
///   the solution vector `x`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is 0.
pub fn tpsv(layout: Layout, uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, ap: anytype, x: anytype, incx: isize) !void {
    comptime var Ap: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(Ap) or !meta.isNumeric(meta.Child(Ap)) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.tpsv: ap must be a many-item pointer to numerics, and x must be a mutable many-item pointer to numerics, got \n\tap: " ++ @typeName(Ap) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    Ap = meta.Child(Ap);
    X = meta.Child(X);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Ap == X) {
        switch (comptime meta.numericType(Ap)) {
            .float => {
                if (comptime Ap == f32)
                    return linalg.cblas.stpsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx)
                else if (comptime Ap == f64)
                    return linalg.cblas.dtpsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx);
            },
            .complex => {
                if (comptime meta.Scalar(Ap) == f32)
                    return linalg.cblas.ctpsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx)
                else if (comptime meta.Scalar(Ap) == f64)
                    return linalg.cblas.ztpsv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx);
            },
            else => {},
        }
    }

    if (layout == .col_major) {
        return if (transa == .no_trans or transa == .trans)
            k_tpsv(uplo, transa, diag, n, ap, x, incx, true)
        else
            k_tpsv(uplo, transa, diag, n, ap, x, incx, false);
    } else {
        return if (transa == .no_trans or transa == .trans)
            k_tpsv(uplo.invert(), transa.invert(), diag, n, ap, x, incx, true)
        else
            k_tpsv(uplo.invert(), transa.invert(), diag, n, ap, x, incx, false);
    }
}

fn k_tpsv(uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, ap: anytype, x: anytype, incx: isize, comptime noconj: bool) void {
    const Ap: type = meta.Child(@TypeOf(ap));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (n == 0)
        return;

    const nounit: bool = diag == .non_unit;

    var kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        if (uplo == .upper) {
            var kk: usize = int.div(n * (n + 1), 2) - 1;
            if (incx == 1) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    if (numeric.ne(x[j], 0)) {
                        if (nounit) {
                            // x[j] /= ap[kk]
                            numeric.div_(
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    ap[kk]
                                else
                                    numeric.conj(ap[kk]),
                            );
                        }

                        const temp = x[j];

                        var k: usize = kk - 1;
                        var i: usize = j;
                        while (i > 0) {
                            i -= 1;

                            // x[i] -= temp * ap[k]
                            numeric.fma_(
                                &x[i],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[i],
                            );

                            k -= 1;
                        }
                    }

                    kk -= j + 1;
                }
            } else {
                var jx: usize = kx + (n - 1) * incx;
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        if (nounit) {
                            // x[jx] /= ap[kk]
                            numeric.div_(
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    ap[kk]
                                else
                                    numeric.conj(ap[kk]),
                            );
                        }

                        const temp = x[numeric.cast(usize, jx)];

                        var ix: isize = jx;
                        var k: usize = kk - 1;
                        while (k >= kk - j) : (k -= 1) {
                            ix -= incx;

                            // x[ix] -= temp * ap[k]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[numeric.cast(usize, ix)],
                            );
                        }
                    }

                    jx -= incx;
                    kk -= j + 1;
                }
            }
        } else {
            var kk: usize = 0;
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    if (numeric.ne(x[j], 0)) {
                        if (nounit) {
                            // x[j] /= ap[kk]
                            numeric.div_(
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    ap[kk]
                                else
                                    numeric.conj(ap[kk]),
                            );
                        }

                        const temp = x[j];

                        var k: usize = kk + 1;
                        var i: usize = j + 1;
                        while (i < n) : (i += 1) {
                            // x[i] -= temp * ap[k]
                            numeric.fma_(
                                &x[i],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[i],
                            );

                            k += 1;
                        }
                    }

                    kk += n - j;
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        if (nounit) {
                            // x[jx] /= ap[kk]
                            numeric.div_(
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    ap[kk]
                                else
                                    numeric.conj(ap[kk]),
                            );
                        }

                        const temp = x[numeric.cast(usize, jx)];

                        var ix: usize = jx;
                        var k: usize = kk + 1;
                        while (k < kk + n - j) : (k += 1) {
                            ix += incx;

                            // x[ix] -= temp * ap[k]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                numeric.neg(temp),
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[numeric.cast(usize, ix)],
                            );
                        }
                    }

                    jx += incx;
                    kk += n - j;
                }
            }
        }
    } else {
        if (uplo == .upper) {
            var kk: usize = 0;
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));
                    numeric.set(&temp, x[j]);

                    var k: usize = kk;
                    var i: usize = 0;
                    while (i < j) : (i += 1) {
                        // temp -= ap[k] * x[i]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            numeric.neg(x[i]),
                            temp,
                        );

                        k += 1;
                    }

                    if (nounit) {
                        // temp /= ap[kk + j]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                ap[kk + j]
                            else
                                numeric.conj(ap[kk + j]),
                        );
                    }

                    // x[j] = temp
                    numeric.set(
                        &x[j],
                        temp,
                    );

                    kk += j + 1;
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));
                    numeric.set(&temp, x[numeric.cast(usize, jx)]);

                    var ix: isize = kx;
                    var k: usize = kk;
                    while (k < kk + j) : (k += 1) {
                        // temp -= ap[k] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            numeric.neg(x[numeric.cast(usize, ix)]),
                            temp,
                        );

                        ix += incx;
                    }

                    if (nounit) {
                        // temp /= ap[kk + j]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                ap[kk + j]
                            else
                                numeric.conj(ap[kk + j]),
                        );
                    }

                    // x[jx] = temp
                    numeric.set(
                        &x[numeric.cast(usize, jx)],
                        temp,
                    );

                    jx += incx;
                    kk += j + 1;
                }
            }
        } else {
            var kk: usize = int.div(n * (n + 1), 2) - 1;
            if (incx == 1) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));
                    numeric.set(&temp, x[j]);

                    var k: usize = kk;
                    var i: usize = n;
                    while (i > j) {
                        i -= 1;

                        // temp -= ap[k] * x[i + 1]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            numeric.neg(x[i + 1]),
                            temp,
                        );

                        k -= 1;
                    }

                    if (nounit) {
                        // temp /= ap[kk - (n - 1) + j]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                ap[kk - (n - 1) + j]
                            else
                                numeric.conj(ap[kk - (n - 1) + j]),
                        );
                    }

                    // x[j] = temp
                    numeric.set(
                        &x[j],
                        temp,
                    );

                    kk -= n - j;
                }
            } else {
                kx += (n - 1) * incx;

                var jx: isize = kx;
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));
                    numeric.set(&temp, x[numeric.cast(usize, jx)]);

                    var ix: isize = kx;
                    var k: usize = kk;
                    while (k > kk - (n - (j + 1))) : (k -= 1) {
                        // temp -= ap[k] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            numeric.neg(x[numeric.cast(usize, ix)]),
                            temp,
                        );

                        ix -= incx;
                    }

                    if (nounit) {
                        // temp /= ap[kk - (n - 1) + j]
                        numeric.div_(
                            &temp,
                            temp,
                            if (comptime noconj)
                                ap[kk - (n - 1) + j]
                            else
                                numeric.conj(ap[kk - (n - 1) + j]),
                        );
                    }

                    // x[jx] = temp
                    numeric.set(
                        &x[numeric.cast(usize, jx)],
                        temp,
                    );

                    jx -= incx;
                    kk -= n - j;
                }
            }
        }
    }

    return;
}
