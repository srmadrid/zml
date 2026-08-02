const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

pub fn ilalc(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
) !usize {
    if (m == 0 or n == 0)
        return 0;

    if (numeric.ne(a[utils.index(layout, 0, n - 1, lda)], 0) or
        numeric.ne(a[utils.index(layout, m - 1, n - 1, lda)], 0))
    {
        return n;
    } else {
        var j: usize = n;
        while (j > 0) {
            j -= 1;

            var i: usize = 0;
            while (i < m) : (i += 1) {
                if (numeric.ne(a[utils.index(layout, i, j, lda)], 0))
                    return j + 1;
            }
        }

        return 0;
    }
}
