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

    const m: u32 = a.rows;
    const n: u32 = a.cols;

    var result: vector.Dense(C) = try .init(allocator, m);
    errdefer result.deinit(allocator);

    try blas.gbmv(
        types.orderOf(A),
        .no_trans,
        types.scast(i32, m),
        types.scast(i32, n),
        types.scast(i32, a.lower),
        types.scast(i32, a.upper),
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
