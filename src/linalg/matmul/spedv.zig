const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const ops = @import("../../ops.zig");

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

    var i: u32 = 0;
    while (i < n) : (i += 1) {
        try ops.set(
            &result.data[i],
            b.get(a.data[i]) catch unreachable,
            ctx,
        );
    }

    return result;
}
