const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const int = @import("../../int.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");

const blas = @import("../blas.zig");

pub inline fn vm(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = types.Coerce(types.Numeric(A), types.Numeric(B));

    const m: u32 = b.rows;
    const n: u32 = b.cols;

    var result: vector.Dense(C) = try .full(allocator, n, 0, ctx);
    errdefer result.deinit(allocator);

    try blas.copy(
        types.scast(i32, m),
        a.data,
        a.inc,
        result.data,
        result.inc,
        ctx,
    );

    var i: u32 = 0;
    while (i < int.min(m, n)) : (i += 1) {
        try ops.mul_( // result[i] *= b[i, i]
            &result.data[i],
            result.data[i],
            b.data[i],
            ctx,
        );
    }

    while (i < n) : (i += 1) {
        result.data[i] = try constants.zero(C, ctx);
    }

    return result;
}
