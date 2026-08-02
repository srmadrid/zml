const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

pub fn ilalr(
    layout: matrix.Layout,
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
) !usize {
    if (m == 0 or n == 0)
        return 0;

    if (numeric.ne(a[utils.index(layout, m - 1, 0, lda)], 0) or
        numeric.ne(a[utils.index(layout, m - 1, n - 1, lda)], 0))
    {
        return m;
    } else {
        var result: usize = 0;

        var j: usize = 0;
        while (j < n) : (j += 1) {
            var i: usize = m;

            while (i > 0) {
                i -= 1;

                if (numeric.ne(a[utils.index(layout, int.max(0, i), j, lda)], 0)) {
                    result = i + 1;
                    break;
                }
            }
        }

        return result;
    }
}
