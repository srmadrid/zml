const options = @import("options");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;
const Diag = meta.Diag;

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");

/// Computes a matrix-vector product using a triangular packed matrix defined as:
///
/// ```zig
/// x = A * x,
/// ```
///
/// or
///
/// ```zig
/// x = conj(A) * x,
/// ```
///
/// or
///
/// ```zig
/// x = Aᵀ * x,
/// ```
///
/// or
///
/// ```zig
/// x = Aᴴ * x,
/// ```
///
/// where `x` is an `n`-element vector, and `A` is an `n`-by-`n` unit, or
/// non-unit, upper or lower triangular matrix, supplied in packed form.
///
/// If the `link_cblas` option is not `null`, the function will try to call the
/// corresponding CBLAS function, if available.
///
/// ## Signature
/// ```zig
/// linalg.blas.tpmv(layout: Layout, uplo: Uplo, transa: Transpose, diag: Diag, n: usize, ap: [*]const Ap, x: [*]X, incx: isize, ctx: anytype) !void
/// ```
///
/// ## Arguments
/// * `layout` (`Layout`): Specifies whether two-dimensional array storage is
///   col-major or row-major.
/// * `uplo` (`Uplo`): Specifies whether the matrix `A` is an upper or lower
///   triangular matrix.
/// * `transa` (`linalg.Transpose`): Specifies the operation to be performed on
///   `A`:
///   * `no_transpose`: `x = A * x`
///   * `transpose`: `x = Aᵀ * x`
///   * `conj_no_transpose`: `x = conj(A) * x`
///   * `conj_transpose`: `x = Aᴴ * x`
/// * `diag` (`Diag`): Specifies whether the matrix `A` is unit triangular.
/// * `n` (`usize`): Specifies the size of the matrix `A`.
/// * `ap` (`anytype`): Many-item pointer, size at least `(n * (n + 1)) / 2`.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 1) * abs(incx)`. On return, contains the result of the
///   operation.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.blas.Error.InvalidArgument`: If `incx` is 0.
pub fn tpmv(layout: Layout, uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, ap: anytype, x: anytype, incx: isize) !void {
    comptime var Ap: type = @TypeOf(ap);
    comptime var X: type = @TypeOf(x);

    comptime if (!meta.isManyItemPointer(Ap) or !meta.isNumeric(meta.Child(Ap)) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)))
        @compileError("zsl.linalg.blas.tpmv: ap must be a many-item pointer to numerics, and x must be a mutable many-item pointer to numerics, got \n\tap: " ++ @typeName(Ap) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    Ap = meta.Child(Ap);
    X = meta.Child(X);

    if (incx == 0)
        return linalg.blas.Error.InvalidArgument;

    if (comptime options.link_cblas != null and Ap == X) {
        switch (comptime meta.numericType(Ap)) {
            .float => {
                if (comptime Ap == f32)
                    return linalg.cblas.stpmv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx)
                else if (comptime Ap == f64)
                    return linalg.cblas.dtpmv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx);
            },
            .complex => {
                if (comptime meta.Scalar(Ap) == f32)
                    return linalg.cblas.ctpmv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx)
                else if (comptime meta.Scalar(Ap) == f64)
                    return linalg.cblas.ztpmv(layout.toInt(c_int), uplo.toInt(c_int), transa.toInt(c_int), diag.toInt(c_int), numeric.cast(isize, n), ap, x, incx);
            },
            else => {},
        }
    }

    if (layout == .col_major) {
        return if (transa == .no_trans or transa == .trans)
            k_tpmv(uplo, transa, diag, n, ap, x, incx, true)
        else
            k_tpmv(uplo, transa, diag, n, ap, x, incx, false);
    } else {
        return if (transa == .no_trans or transa == .trans)
            k_tpmv(uplo.invert(), transa.invert(), diag, n, ap, x, incx, true)
        else
            k_tpmv(uplo.invert(), transa.invert(), diag, n, ap, x, incx, false);
    }
}

fn k_tpmv(uplo: Uplo, transa: linalg.Transpose, diag: Diag, n: usize, ap: anytype, x: anytype, incx: isize, comptime noconj: bool) void {
    const Ap: type = meta.Child(@TypeOf(ap));
    const X: type = meta.Child(@TypeOf(x));

    // Quick return if possible.
    if (n == 0)
        return;

    const nounit: bool = diag == .non_unit;

    var kx: isize = if (incx < 0) (-numeric.cast(isize, n) + 1) * incx else 0;

    if (transa == .no_trans or transa == .conj_no_trans) {
        if (uplo == .upper) {
            var kk: usize = 0;
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    if (numeric.ne(x[j], 0)) {
                        const temp = x[j];

                        var k: usize = kk;
                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            // x[i] += temp * ap[k]
                            numeric.fma_(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[i],
                            );

                            k += 1;
                        }

                        if (nounit) {
                            // x[j] *= ap[kk + j]
                            numeric.mul_(
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    ap[kk + j]
                                else
                                    numeric.conj(ap[kk + j]),
                            );
                        }
                    }

                    kk += j + 1;
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        const temp = x[numeric.cast(usize, jx)];

                        var ix: isize = kx;
                        var k: usize = kk;
                        while (k < kk + j) : (k += 1) {
                            // x[ix] += temp * ap[k]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix += incx;
                        }

                        if (nounit) {
                            // x[jx] *= ap[kk + j]
                            numeric.mul_(
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    ap[kk + j]
                                else
                                    numeric.conj(ap[kk + j]),
                            );
                        }
                    }

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

                    if (numeric.ne(x[j], 0)) {
                        const temp = x[j];

                        var k: usize = kk;
                        var i: usize = n - 1;
                        while (i > j) : (i -= 1) {
                            // x[i] += temp * ap[k]
                            numeric.fma_(
                                &x[i],
                                temp,
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[i],
                            );

                            k -= 1;
                        }
                        if (nounit) {
                            // x[j] *= ap[kk - (n - 1) + j]
                            numeric.mul_(
                                &x[j],
                                x[j],
                                if (comptime noconj)
                                    ap[kk - (n - 1) + j]
                                else
                                    numeric.conj(ap[kk - (n - 1) + j]),
                            );
                        }
                    }

                    kk -= n - j;
                }
            } else {
                kx += (numeric.cast(isize, n) - 1) * incx;

                var jx: isize = kx;
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    if (numeric.ne(x[numeric.cast(usize, jx)], 0)) {
                        const temp = x[numeric.cast(usize, jx)];

                        var ix: isize = kx;
                        var k: usize = kk;
                        while (k > kk - (n - (j + 1))) : (k -= 1) {
                            // x[ix] += temp * ap[k]
                            numeric.fma_(
                                &x[numeric.cast(usize, ix)],
                                temp,
                                if (comptime noconj)
                                    ap[k]
                                else
                                    numeric.conj(ap[k]),
                                x[numeric.cast(usize, ix)],
                            );

                            ix -= incx;
                        }

                        if (nounit) {
                            // x[jx] *= ap[kk - (n - 1) + j]
                            numeric.mul_(
                                &x[numeric.cast(usize, jx)],
                                x[numeric.cast(usize, jx)],
                                if (comptime noconj)
                                    ap[kk - (n - 1) + j]
                                else
                                    numeric.conj(ap[kk - (n - 1) + j]),
                            );
                        }
                    }

                    jx -= incx;
                    kk -= n - j;
                }
            }
        }
    } else {
        if (uplo == .upper) {
            var kk: usize = int.div(n * (n + 1), 2) - 1;
            if (incx == 1) {
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                    // temp += ap[kk] * x[j]
                    numeric.fma_(
                        &temp,
                        if (nounit)
                            if (comptime noconj)
                                ap[kk]
                            else
                                numeric.conj(ap[kk])
                        else
                            numeric.one(Ap),
                        x[j],
                        temp,
                    );

                    var k: usize = kk - 1;
                    var i: usize = j;
                    while (i > 0) {
                        i -= 1;

                        // temp += ap[k] * x[i]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            x[i],
                            temp,
                        );

                        k -= 1;
                    }

                    // x[j] = temp
                    numeric.set(
                        &x[j],
                        temp,
                    );

                    kk -= j + 1;
                }
            } else {
                var jx: isize = kx + (numeric.cast(isize, n) - 1) * incx;
                var j: usize = n;
                while (j > 0) {
                    j -= 1;

                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                    // temp += ap[kk] * x[jx]
                    numeric.fma_(
                        &temp,
                        if (nounit)
                            if (comptime noconj)
                                ap[kk]
                            else
                                numeric.conj(ap[kk])
                        else
                            numeric.one(Ap),
                        x[numeric.cast(usize, jx)],
                        temp,
                    );

                    var ix: isize = jx;
                    var k: usize = kk - 1;
                    while (k >= kk - j) : (k -= 1) {
                        ix -= incx;

                        // temp += ap[k] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            x[numeric.cast(usize, ix)],
                            temp,
                        );
                    }

                    // x[jx] = temp
                    numeric.set(
                        &x[numeric.cast(usize, jx)],
                        temp,
                    );

                    jx -= incx;
                    kk -= j + 1;
                }
            }
        } else {
            var kk: usize = 0;
            if (incx == 1) {
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                    // temp += ap[kk] * x[j]
                    numeric.fma_(
                        &temp,
                        if (nounit)
                            if (comptime noconj)
                                ap[kk]
                            else
                                numeric.conj(ap[kk])
                        else
                            numeric.one(Ap),
                        x[j],
                        temp,
                    );

                    var k: usize = kk + 1;
                    var i: usize = j + 1;
                    while (i < n) : (i += 1) {
                        // temp += ap[k] * x[i]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            x[i],
                            temp,
                        );

                        k += 1;
                    }

                    // x[j] = temp
                    numeric.set(
                        &x[j],
                        temp,
                    );

                    kk += n - j;
                }
            } else {
                var jx: isize = kx;
                var j: usize = 0;
                while (j < n) : (j += 1) {
                    var temp = numeric.zero(meta.Accumulator(numeric.Mul(if (comptime noconj) Ap else numeric.Conj(Ap), X)));

                    // temp += ap[kk] * x[jx]
                    numeric.fma_(
                        &temp,
                        if (nounit)
                            if (comptime noconj)
                                ap[kk]
                            else
                                numeric.conj(ap[kk])
                        else
                            numeric.one(Ap),
                        x[numeric.cast(usize, jx)],
                        temp,
                    );

                    var ix: isize = jx;
                    var k: usize = kk + 1;
                    while (k < kk + n - j) : (k += 1) {
                        ix += incx;

                        // temp += ap[k] * x[ix]
                        numeric.fma_(
                            &temp,
                            if (comptime noconj)
                                ap[k]
                            else
                                numeric.conj(ap[k]),
                            x[numeric.cast(usize, ix)],
                            temp,
                        );
                    }

                    // x[jx] = temp
                    numeric.set(
                        &x[numeric.cast(usize, jx)],
                        temp,
                    );

                    jx += incx;
                    kk += n - j;
                }
            }
        }
    }

    return;
}
