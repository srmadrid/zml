const std = @import("std");

const types = @import("../../types.zig");
const Scalar = types.Scalar;
const Child = types.Child;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn asum(
    n: isize,
    x: anytype,
    incx: isize,
    ctx: anytype,
) !Scalar(Child(@TypeOf(x))) {
    const X: type = types.Child(@TypeOf(x));

    var sum: Scalar(X) = try ops.init(Scalar(X), ctx);
    errdefer ops.deinit(&sum, ctx);

    try @import("asum_sub.zig").asum_sub(n, x, incx, &sum, ctx);

    return sum;
}
