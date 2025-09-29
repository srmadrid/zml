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
    const k: u32 = a.cols;
    const n: u32 = b.cols;

    var result: matrix.dense.General(C, types.orderOf(A)) = try .init(allocator, m, n);
    errdefer result.deinit(allocator);

    try blas.gemm(
        types.orderOf(A),
        .no_trans,
        if (comptime types.orderOf(A) == types.orderOf(B)) .no_trans else .trans,
        types.scast(i32, m),
        types.scast(i32, n),
        types.scast(i32, k),
        1,
        a.data,
        types.scast(i32, a.ld),
        b.data,
        types.scast(i32, b.ld),
        0,
        result.data,
        types.scast(i32, result.ld),
        ctx,
    );

    return result;
}
