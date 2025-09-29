const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");

const blas = @import("../blas.zig");

pub inline fn mv(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = types.Coerce(types.Numeric(A), types.Numeric(B));

    const n: u32 = a.size;

    var result: vector.Dense(C) = try .init(allocator, n);
    errdefer result.deinit(allocator);

    try blas.hemv(
        types.orderOf(A),
        types.uploOf(A),
        types.scast(i32, n),
        1,
        a.data,
        types.scast(i32, a.ld),
        b.data,
        b.inc,
        0,
        result.data,
        result.inc,
        ctx,
    );

    return result;
}
