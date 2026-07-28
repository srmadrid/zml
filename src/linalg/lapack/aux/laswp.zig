const std = @import("std");

const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Performs a series of row interchanges on a general matrix, defined as:
///
/// ```zig
/// Aᵢⱼ ↔ Aₚ₍ᵢ₎ⱼ
/// ```
///
/// for all `i = k₀, .., k₁ − 1` and `j = 0, ..., n - 1`, where `p(i) = ipivᵢ`
/// and `A` is an `m × n` matrix.
///
/// ## Signature
/// ```zig
/// linalg.lapack.laswp(layout: matrix.Layout, n: usize, a: [*]A, lda: usize, k0: usize, k1: usize, ipiv: [*]const usize, inci: isize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `n` (`usize`): Specifies the number of columns of the matrix `A`.
/// * `a` (`anytype`): Mutable many-item pointer, size at least `lda * k`, where
///   `k` is `n` when `layout` is `col_major`, or `m` when `layout` is
///   `row_major`. `m` is not less than the maximum of
///   `ipiv[k0 + j * abs(inci)]`, for  `j = 0..k1 - k0 - 1`.
/// * `lda` (`usize`): Specifies the leading dimension of `c` as declared in the
///   calling (sub)program. Must be greater than or equal to `max(1, m)` when
///   `layout` is `col_major`, or `max(1, n)` when `layout` is `row_major`.
/// * `k0` (`usize`): The starting index (inclusive) in `ipiv` at which row
///   interchanges begin.
/// * `k1` (`usize`): The ending index (exclusive) in `ipiv`, defining the
///   half-open range `k0..k1` of row interchanges to perform.
/// * `ipiv` (`[*]usize`): Mutable many-item pointer, size at least
///   `k0 + (k1 - k0 - 1) * abs(inci) + 1`. Only elements in the range `k0..k1`
///   of `ipiv` are accessed. `ipiv[i] = j` implies rows `i` and `j` are to be
///   interchanged.
/// * `inci` (`isize`): Indexing increment for `ipiv`. Must be different from 0.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `linalg.lapack.Error.InvalidArgument`: If `lda` is less than `max(1, m)`
///   or `max(1, n)`, or if `inci` is 0.
pub fn laswp(
    layout: matrix.Layout,
    n: usize,
    a: anytype,
    lda: usize,
    k0: usize,
    k1: usize,
    ipiv: [*]const usize,
    inci: isize,
) !void {
    const A: type = @TypeOf(a);

    comptime if (!meta.isManyItemPointer(A) or meta.isConstPointer(A) or !meta.isNumeric(meta.Child(A)))
        @compileError("zsl.linalg.lapack.getrf2: a must be a mutable many-item pointer to numerics, got \n\ta: " ++ @typeName(A) ++ "\n");

    if ((layout == .row_major and lda < int.max(1, n)) or inci == 0)
        return linalg.lapack.Error.InvalidArgument;

    if (k0 >= k1 or n == 0)
        return;

    return k_laswp(layout, n, a, lda, k0, k1, ipiv, inci);
}

fn k_laswp(
    layout: matrix.Layout,
    n: usize,
    a: anytype,
    lda: usize,
    k0: usize,
    k1: usize,
    ipiv: [*]const usize,
    inci: isize,
) !void {
    const A: type = meta.Child(@TypeOf(a));

    const unroll = 2 * (std.simd.suggestVectorLength(A) orelse 2);

    const ii0: isize = if (inci < 0)
        numeric.cast(isize, k0) - (numeric.cast(isize, k1 - k0 - 1) * inci)
    else
        numeric.cast(isize, k0);
    const I0: isize = if (inci < 0) numeric.cast(isize, k1) - 1 else numeric.cast(isize, k0);
    const I1: isize = if (inci < 0) numeric.cast(isize, k0) - 1 else numeric.cast(isize, k1);
    const inc: isize = if (inci < 0) -1 else 1;

    var j: usize = 0;
    while (j < (n / unroll) * unroll) : (j += unroll) {
        var ii: isize = ii0;
        var i: isize = I0;
        while (i != I1) {
            const ip: usize = ipiv[numeric.cast(usize, ii)];
            const r: usize = numeric.cast(usize, i);
            if (ip != r) {
                inline for (0..unroll) |u| {
                    const temp = a[utils.index(layout, r, j + u, lda)];
                    a[utils.index(layout, r, j + u, lda)] = a[utils.index(layout, ip, j + u, lda)];
                    a[utils.index(layout, ip, j + u, lda)] = temp;
                }
            }

            i += inc;
            ii += inci;
        }
    }

    while (j < n) : (j += 1) {
        var ii: isize = ii0;
        var i: isize = I0;
        while (i != I1) {
            const ip: usize = ipiv[numeric.cast(usize, ii)];
            const r: usize = numeric.cast(usize, i);
            if (ip != r) {
                const temp = a[utils.index(layout, r, j, lda)];
                a[utils.index(layout, r, j, lda)] = a[utils.index(layout, ip, j, lda)];
                a[utils.index(layout, ip, j, lda)] = temp;
            }

            i += inc;
            ii += inci;
        }
    }
}
