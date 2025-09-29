const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const int = @import("../../int.zig");

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

    if (m == n) {
        try blas.copy(
            types.scast(i32, n),
            b.data,
            b.inc,
            result.data,
            result.inc,
            ctx,
        );

        try blas.trmv(
            types.orderOf(A),
            types.uploOf(A),
            .no_trans,
            types.diagOf(A),
            types.scast(i32, n),
            a.data,
            types.scast(i32, a.ld),
            result.data,
            result.inc,
            ctx,
        );
    } else {
        const min_dim: u32 = int.min(m, n);

        try blas.copy(
            types.scast(i32, min_dim),
            b.data,
            b.inc,
            result.data,
            result.inc,
            ctx,
        );

        try blas.trmv(
            types.orderOf(A),
            types.uploOf(A),
            .no_trans,
            types.diagOf(A),
            types.scast(i32, min_dim),
            a.data,
            types.scast(i32, a.ld),
            result.data,
            result.inc,
            ctx,
        );

        if (comptime types.uploOf(A) == .upper) {
            if (n > min_dim) { // extra rows (filled)
                try blas.gemv(
                    types.orderOf(A),
                    .no_trans,
                    types.scast(i32, m),
                    types.scast(i32, n - min_dim),
                    1,
                    a.data +
                        if (comptime types.orderOf(A) == .col_major)
                            min_dim * a.ld
                        else
                            min_dim,
                    types.scast(i32, a.ld),
                    b.data +
                        if (b.inc > 0)
                            min_dim * types.scast(u32, b.inc)
                        else
                            0,
                    b.inc,
                    1,
                    result.data,
                    result.inc,
                    ctx,
                );
            } else if (m > min_dim) { // extra rows (empty)
                try blas.scal(
                    types.scast(i32, m - min_dim),
                    0,
                    result.data + min_dim * types.scast(u32, result.inc),
                    result.inc,
                    ctx,
                );
            }
        } else {
            if (m > min_dim) { // extra rows (filled)
                try blas.gemv(
                    types.orderOf(A),
                    .no_trans,
                    types.scast(i32, m - min_dim),
                    types.scast(i32, n),
                    1,
                    a.data +
                        if (comptime types.orderOf(A) == .col_major)
                            min_dim
                        else
                            min_dim * a.ld,
                    types.scast(i32, a.ld),
                    b.data,
                    b.inc,
                    1,
                    result.data + min_dim * types.scast(u32, result.inc),
                    result.inc,
                    ctx,
                );
            }
        }
    }

    return result;
}
