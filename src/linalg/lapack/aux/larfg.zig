const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Generates an elementary reflector `H` of size `n`, such that:
///
/// ```zig
///      ┌      ┐
///      │ α  β │
/// Hᴴ = │ x  0 │,  Hᴴ * H = 𝕀,
///      └      ┘
/// ```
///
/// where `α` and `β` are numerics, with `β` real, and `x` is an `n - 1`-element
/// vector. `H` is represented in the form:
///
/// ```zig
///             ┌   ┐
///             │ 1 │   ┌       ┐
/// H = 𝕀 - τ * │ v │ * │ 1  vᴴ │,
///             └   ┘   └       ┘
/// ```
///
/// where `τ` is a numeric and `v` is an `n - 1`-element vector. Note that `H`
/// is not Hermitian.
///
/// If the elements of `x` are all zero and `α` is real, then `τ = 0` and `H` is
/// taken to be the identity matrix. Otherwise `1 <= Re(tau) <= 2` and
/// `|τ - 1| <= 1`.
///
/// ## Signature
/// ```zig
/// linalg.lapack.larfg(n: usize, alpha: *Al, x: [*]X, incx: isize, tau: *Ta) !void
/// ```
///
/// ## Arguments
/// * `n` (`usize`): Specifies the size of the elementary reflector `H`.
/// * `alpha` (`anytype`): Mutable one-item pointer.
/// * `x` (`anytype`): Mutable many-item pointer, size at least
///   `1 + (n - 2) * abs(incx)`.
/// * `incx` (`isize`): Indexing increment for `x`. Must be different from 0.
/// * `tau` (`anytype`): Mutable one-item pointer.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `incx` is 0.
pub fn larfg(
    n: usize,
    alpha: anytype,
    x: anytype,
    incx: isize,
    tau: anytype,
) !void {
    const Al: type = @TypeOf(alpha);
    const X: type = @TypeOf(x);
    const Ta: type = @TypeOf(tau);

    comptime if (!meta.isPointer(Al) or meta.isConstPointer(Al) or !meta.isNumeric(meta.Child(Al)) or
        !meta.isManyItemPointer(X) or meta.isConstPointer(X) or !meta.isNumeric(meta.Child(X)) or
        !meta.isPointer(Ta) or meta.isConstPointer(Ta) or !meta.isNumeric(meta.Child(Ta)))
        @compileError("zsl.linalg.lapack.larfg: alpha and tau must be mutable one-item pointers to a numeric, and x must be a mutable many-item pointer to numerics, got \n\talpha: " ++ @typeName(Al) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ttau: " ++ @typeName(Ta) ++ "\n");

    if (incx == 0)
        return linalg.lapack.Error.InvalidArgument;

    if (n == 0) {
        numeric.set(tau, 0);

        return;
    }

    var xnorm = linalg.blas.nrm2(n - 1, x, incx) catch unreachable;
    if (numeric.eq(xnorm, 0) and numeric.eq(numeric.im(alpha.*), 0)) {
        // H = 𝕀
        numeric.set(tau, 0);
    } else {
        // General case.
        var beta = numeric.neg(
            numeric.copysign(
                numeric.hypot(
                    numeric.hypot(
                        numeric.re(alpha.*),
                        numeric.im(alpha.*),
                    ),
                    xnorm,
                ),
                numeric.re(alpha.*),
            ),
        );

        const safmin = numeric.div(
            numeric.smallest(meta.Real(meta.Child(X))),
            numeric.eps(meta.Real(meta.Child(X))),
        );
        const rsafmn = numeric.div(1, safmin);

        var knt: usize = 0;
        if (numeric.lt(numeric.abs(beta), safmin)) {
            // xnorm and beta may be inaccurate; scale x and recompute them.
            while (numeric.lt(numeric.abs(beta), safmin) and knt < 20) : (knt += 1) {
                linalg.blas.scal(
                    n - 1,
                    rsafmn,
                    x,
                    incx,
                ) catch unreachable;

                numeric.mulInto(&beta, beta, rsafmn);
                numeric.mulInto(alpha, alpha.*, rsafmn);
            }

            // New beta is at most 1 and at least safmin.
            xnorm = linalg.blas.nrm2(n - 1, x, incx) catch unreachable;

            beta = numeric.neg(
                numeric.copysign(
                    numeric.hypot(
                        numeric.hypot(
                            numeric.re(alpha.*),
                            numeric.im(alpha.*),
                        ),
                        xnorm,
                    ),
                    numeric.re(alpha.*),
                ),
            );
        }

        numeric.set(
            tau,
            numeric.div(
                numeric.sub(beta, alpha.*),
                beta,
            ),
        );

        numeric.divInto(alpha, 1, numeric.sub(alpha.*, beta));

        linalg.blas.scal(
            n - 1,
            alpha.*,
            x,
            incx,
        ) catch unreachable;

        // If alpha is subnormal, it may lose relative accuracy.
        for (0..knt) |_| {
            numeric.mulInto(&beta, beta, safmin);
        }

        numeric.set(alpha, beta);
    }

    return;
}
