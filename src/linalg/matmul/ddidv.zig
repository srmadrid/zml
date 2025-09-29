const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const int = @import("../../int.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");

const blas = @import("../blas.zig");

pub inline fn mv(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = types.Coerce(types.Numeric(A), types.Numeric(B));

    const m: u32 = a.rows;
    const n: u32 = a.cols;

    var result: vector.Dense(C) = try .full(allocator, m, 0, ctx);
    errdefer result.deinit(allocator);

    try blas.copy(
        types.scast(i32, m),
        b.data,
        b.inc,
        result.data,
        result.inc,
        ctx,
    );

    var i: u32 = 0;
    while (i < int.min(m, n)) : (i += 1) {
        try ops.mul_( // result[i] *= a[i, i]
            &result.data[i],
            result.data[i],
            a.data[i],
            ctx,
        );
    }

    while (i < m) : (i += 1) {
        result.data[i] = try constants.zero(C, ctx);
    }

    return result;
}
