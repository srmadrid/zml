const std = @import("std");

const types = @import("../../types.zig");
const Child = types.Child;
const Coerce = types.Coerce;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn dot(
    n: isize,
    x: anytype,
    incx: isize,
    y: anytype,
    incy: isize,
    ctx: anytype,
) !Coerce(Child(@TypeOf(x)), Child(@TypeOf(y))) {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = types.Coerce(X, Y);

    var sum: C = try ops.init(C, ctx);
    errdefer ops.deinit(&sum, ctx);

    try @import("dot_sub.zig").dot_sub(n, x, incx, y, incy, &sum, ctx);

    return sum;
}
