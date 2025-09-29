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

    const n: u32 = b.size;

    var result: vector.Dense(C) = try .init(allocator, n);
    errdefer result.deinit(allocator);

    try blas.symv(
        types.orderOf(B).invert(),
        types.uploOf(B).invert(),
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
