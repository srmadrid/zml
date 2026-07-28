const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Uplo = types.Uplo;

const utils = @import("../utils.zig");

pub fn ilalc(
    order: Order,
    m: i32,
    n: i32,
    a: anytype,
    lda: i32,
) !i32 {
    if (m == 0 or n == 0) {
        return 0;
    } else if (try ops.ne(a[utils.index(order, 0, n - 1, lda)], 0, .{}) or
        try ops.ne(a[utils.index(order, m - 1, n - 1, lda)], 0, .{}))
    {
        return n;
    } else {
        var j: i32 = n - 1;
        while (j >= 0) : (j -= 1) {
            var i: i32 = 0;
            while (i < m) : (i += 1) {
                if (try ops.ne(a[utils.index(order, i, j, lda)], 0, .{}))
                    return j + 1;
            }
        }

        return 0;
    }
}
