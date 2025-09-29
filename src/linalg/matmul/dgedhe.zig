const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");

const blas = @import("../blas.zig");

pub inline fn mm(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = types.Coerce(types.Numeric(A), types.Numeric(B));

    const m: u32 = a.rows;
    const n: u32 = b.size;

    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
    errdefer result.deinit(allocator);

    try blas.hemm(
        types.orderOf(A),
        .right,
        if (comptime types.orderOf(A) == types.orderOf(B)) types.uploOf(B) else types.uploOf(B).invert(),
        types.scast(i32, m),
        types.scast(i32, n),
        1,
        b.data,
        types.scast(i32, b.ld),
        a.data,
        types.scast(i32, a.ld),
        0,
        result.data,
        types.scast(i32, result.ld),
        ctx,
    );

    return result;
}
