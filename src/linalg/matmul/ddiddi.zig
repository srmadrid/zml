const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const int = @import("../../int.zig");
const ops = @import("../../ops.zig");

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");

const blas = @import("../blas.zig");

pub inline fn mm(allocator: std.mem.Allocator, a: anytype, b: anytype, ctx: anytype) !MulCoerce(@TypeOf(a), @TypeOf(b)) {
    const A: type = @TypeOf(a);
    const B: type = @TypeOf(b);
    const C: type = types.Coerce(types.Numeric(A), types.Numeric(B));

    const m: u32 = a.rows;
    const n: u32 = b.cols;

    var result: matrix.dense.Diagonal(C) = try .init(allocator, m, n);
    errdefer result.deinit(allocator);

    var i: u32 = 0;
    while (i < int.min(a.rows, b.cols)) : (i += 1) {
        result.data[i] = try ops.mul(
            a.data[i],
            b.data[i],
            ctx,
        );
    }

    return result;
}
