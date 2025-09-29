const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");

const blas = @import("../blas.zig");

pub inline fn vm(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = types.Coerce(types.Numeric(A), types.Numeric(B));

    const m: u32 = b.rows;
    const n: u32 = b.cols;

    var result: vector.Dense(C) = try .init(allocator, n);
    errdefer result.deinit(allocator);

    try blas.gemv(
        types.orderOf(B),
        .trans,
        types.scast(i32, m),
        types.scast(i32, n),
        1,
        b.data,
        types.scast(i32, b.ld),
        a.data,
        a.inc,
        0,
        result.data,
        result.inc,
        ctx,
    );

    return result;
}
