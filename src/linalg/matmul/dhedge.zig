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
    const n: u32 = b.size;

    var result: matrix.General(C, types.orderOf(A)) = try .init(allocator, m, n);
    errdefer result.deinit(allocator);

    if (comptime types.orderOf(A) == types.orderOf(B)) {
        try blas.hemm(
            types.orderOf(B),
            .left,
            types.uploOf(A),
            types.scast(i32, m),
            types.scast(i32, n),
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
    } else {
        var j: u32 = 0;
        while (j < n) : (j += 1) {
            try blas.hemv(
                types.orderOf(A),
                types.uploOf(A),
                types.scast(i32, m),
                1,
                a.data,
                types.scast(i32, a.ld),
                b.data +
                    if (comptime types.orderOf(B) == .col_major)
                        j * b.ld
                    else
                        j,
                if (comptime types.orderOf(B) == .col_major) 1 else types.scast(i32, b.ld),
                0,
                result.data +
                    if (comptime types.orderOf(A) == .col_major)
                        j * result.ld
                    else
                        j,
                if (comptime types.orderOf(A) == .col_major) 1 else types.scast(i32, result.ld),
                ctx,
            );
        }
    }

    return result;
}
