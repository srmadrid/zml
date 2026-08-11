const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Applies an elementary reflector `H` to an `m × n` matrix `C`, from either
/// the left or the right. `H` is represented in the form:
///
/// ```zig
/// H = 𝕀 - τ v vᴴ,
/// ```
///
/// where tau is a `τ` numeric and `v` is a vector.
///
/// If `τ = 0`, then `H` is taken to be the identity matrix.
///
/// To apply `Hᴴ`, use `̅τ` instead of `τ`.
///
/// ## Signature
/// ```zig
/// linalg.lapack.larf1f(layout: matrix.Layout, m: usize, n: usize, v: [*]const V, incv: isize, tau: Ta, c: [*]C, ldc: usize, work: [*]W) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `side` (`linalg.blas.Side`): Specifies whether the elementary reflector is
///   applied from the left or from the right:
///   * `left`: `C = H * C`.
///   * `right`: `C = C * H`.
/// * `m` (`usize`): Specifies the number of rows of the matrix `C`.
/// * `n` (`usize`): Specifies the number of columns of the matrix `C`.
/// * `v` (`anytype`): Many-item pointer, size at least
///   `1 + (m - 1) * abs(incv)` when `side` is `left`, or
///   `1 + (n - 1) * abs(incv)` otherwise.
/// * `incv` (`usize`): Indexing increment for `v`. Must be different from 0.
/// * `tau` (`anytype`): Specifies the numeric `τ`.
/// * `c` (`anytype`): Mutable many-item pointer, size at least `lda * k`, where
///   `k` is `n` when `layout` is `col_major`, or `m` when `layout` is
///   `row_major`.
/// * `ldc` (`usize`): Specifies the leading dimension of `c as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
/// * `work` (`anytype`): Mutable many-item pointer, size at least `max(1, n)`
///   when `side` is `left`, or `max(1, m)` otherwise.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `incv` is 0, or if `lda` is less
///   than `max(1, m)` or `max(1, n)`.
pub fn larf1f(
    layout: matrix.Layout,
    side: linalg.blas.Side,
    m: usize,
    n: usize,
    v: anytype,
    incv: usize,
    tau: anytype,
    c: anytype,
    ldc: usize,
    work: anytype,
) !void {
    const V: type = @TypeOf(v);
    const Ta: type = @TypeOf(tau);
    const C: type = @TypeOf(c);

    comptime if (!meta.isManyItemPointer(V) or !meta.isNumeric(meta.Child(V)) or
        !meta.isNumeric(Ta) or
        !meta.isManyItemPointer(C) or meta.isConstPointer(C) or !meta.isNumeric(meta.Child(C)))
        @compileError("zsl.linalg.lapack.larf1f: v must be a many-item pointer to numerics, tau must be a numeric, and c must be a mutable many-item pointer to numerics, got \n\tv: " ++ @typeName(V) ++ "\n\ttau: " ++ @typeName(Ta) ++ "\n\tc: " ++ @typeName(C) ++ "\n");

    if (incv == 0 or ldc < int.max(1, if (layout == .col_major) m else n))
        return linalg.lapack.Error.InvalidArgument;

    var lastv: usize = 1;
    var lastc: usize = 0;
    if (numeric.ne(tau, 0)) {
        // Set up variables for scanning `v`. lastv begins pointing to the end
        // of `v` up to `v[0]`.
        lastv = if (side == .left)
            m
        else
            n;

        // Look for the last non-zero row in `v`.
        var i: usize = (lastv - 1) * incv;
        while (lastv > 1 and numeric.eq(v[i], 0)) {
            lastv -= 1;
            i -= incv;
        }

        if (side == .left)
            // Scan for the last non-zero column in `C[0..lastv, ..]`.
            lastc = linalg.lapack.ilalc(layout, lastv, n, c, ldc) catch unreachable
        else
            // Scan for the last non-zero row in `C[.., 0..lastv]`.
            lastc = linalg.lapack.ilalr(layout, m, lastv, c, ldc) catch unreachable;
    }

    if (lastc == 0)
        return;

    if (side == .left) {
        // Form `H C`.
        if (lastv == 1) {
            // `C[0, 0..lastc - 1] = (1 - tau) C[0, 0..lastc - 1]`.
            linalg.blas.scal(
                lastc,
                numeric.sub(1, tau),
                c,
                utils.row_ld(layout, ldc),
            ) catch unreachable;
        } else {
            // `w[0..lastc - 1] = C[1..lastv - 1, 0..lastc - 1]ᴴ v[1..lastv - 1]`.
            linalg.blas.gemv(
                layout,
                .conj_trans,
                lastv - 1,
                lastc,
                1,
                c + utils.index(layout, 1, 0, ldc),
                ldc,
                v + incv,
                numeric.cast(isize, incv),
                0,
                work,
                1,
            ) catch unreachable;

            // `w[0..lastc - 1] += v[0] * C[0, 0..lastc - 1]ᴴ`.
            for (0..lastc) |i| {
                numeric.addInto(
                    &work[i],
                    work[i],
                    numeric.conj(c[utils.index(layout, 0, i, ldc)]),
                );
            }

            // `C[0, 0..lastc - 1] -= tau w[0..lastc - 1]ᴴ`.
            for (0..lastc) |i| {
                numeric.subInto(
                    &c[utils.index(layout, 0, i, ldc)],
                    c[utils.index(layout, 0, i, ldc)],
                    numeric.mul(tau, numeric.conj(work[i])),
                );
            }

            // `C[1..lastv - 1, 0..lastc - 1] -= tau v[1..lastv - 1] w[0..lastc - 1]ᴴ`.
            linalg.blas.gerc(
                layout,
                lastv - 1,
                lastc,
                numeric.neg(tau),
                v + incv,
                numeric.cast(isize, incv),
                work,
                1,
                c + utils.index(layout, 1, 0, ldc),
                ldc,
            ) catch unreachable;
        }
    } else {
        // Form `C H`.
        if (lastv == 1) {
            // `C[0..lastc - 1, 0] = (1 - tau) C[0..lastc - 1, 0]`.
            linalg.blas.scal(
                lastc,
                numeric.sub(1, tau),
                c,
                utils.col_ld(layout, ldc),
            ) catch unreachable;
        } else {
            // `w[0..lastc - 1] = C[0..lastc - 1, 1..lastv - 1] v[1..lastv - 1]`.
            linalg.blas.gemv(
                layout,
                .no_trans,
                lastc,
                lastv - 1,
                1,
                c + utils.index(layout, 0, 1, ldc),
                ldc,
                v + incv,
                numeric.cast(isize, incv),
                0,
                work,
                1,
            ) catch unreachable;

            // `w[0..lastc - 1] += v[0] C[0..lastc - 1, 0]`.
            linalg.blas.axpy(
                lastc,
                1,
                c,
                utils.col_ld(layout, ldc),
                work,
                1,
            ) catch unreachable;

            // `C[0..lastc - 1, 0] -= tau v[0] w[0..lastc - 1]`.
            linalg.blas.axpy(
                lastc,
                numeric.neg(tau),
                work,
                1,
                c,
                utils.col_ld(layout, ldc),
            ) catch unreachable;

            // `C[0..lastc - 1, 1..lastv - 1] += -tau w[0..lastc - 1] * v[1..lastv - 1]ᴴ`.
            linalg.blas.gerc(
                layout,
                lastc,
                lastv - 1,
                numeric.neg(tau),
                work,
                1,
                v + incv,
                numeric.cast(isize, incv),
                c + utils.index(layout, 0, 1, ldc),
                ldc,
            ) catch unreachable;
        }
    }
}
