const std = @import("std");

const types = @import("../../types.zig");
const Child = types.Child;
const Coerce = types.Coerce;
const ops = @import("../../ops.zig");

const blas = @import("../blas.zig");

pub fn dotc(
    n: i32,
    x: anytype,
    incx: i32,
    y: anytype,
    incy: i32,
    ctx: anytype,
) !Coerce(Child(@TypeOf(x)), Child(@TypeOf(y))) {
    const X: type = types.Child(@TypeOf(x));
    const Y: type = types.Child(@TypeOf(y));
    const C: type = types.Coerce(X, Y);

    var sum: C = try ops.init(C, ctx);
    errdefer ops.deinit(&sum, ctx);

    try @import("dotc_sub.zig").dotc_sub(n, x, incx, y, incy, &sum, ctx);

    return sum;
}
