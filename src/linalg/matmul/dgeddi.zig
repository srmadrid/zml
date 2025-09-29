const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const int = @import("../../int.zig");
const constants = @import("../../constants.zig");

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

    var result: matrix.dense.General(C, types.orderOf(A)) = try .full(allocator, m, n, 0, ctx);
    errdefer result.deinit(allocator);

    // lapack.lascl2

    var j: u32 = 0;
    while (j < int.min(k, n)) : (j += 1) {
        try blas.axpy(
            types.scast(i32, m),
            b.data[j],
            a.data + if (comptime types.orderOf(A) == .col_major) j * a.ld else j,
            if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, a.ld),
            result.data + if (comptime types.orderOf(A) == .col_major) j * result.ld else j,
            if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
            ctx,
        );
    }

    while (j < n) : (j += 1) {
        var i: u32 = 0;
        while (i < m) : (i += 1) {
            if (comptime types.orderOf(A) == .col_major) {
                result.data[i + j * result.ld] = try constants.zero(C, ctx);
            } else {
                result.data[i * result.ld + j] = try constants.zero(C, ctx);
            }
        }
    }

    return result;
}
