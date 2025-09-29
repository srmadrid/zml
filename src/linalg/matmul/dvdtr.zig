const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const int = @import("../../int.zig");

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

    if (m == n) {
        try blas.copy(
            types.scast(i32, n),
            a.data,
            a.inc,
            result.data,
            result.inc,
            ctx,
        );

        try blas.trmv(
            types.orderOf(B),
            types.uploOf(B),
            .trans,
            types.diagOf(B),
            types.scast(i32, n),
            b.data,
            types.scast(i32, b.ld),
            result.data,
            result.inc,
            ctx,
        );
    } else {
        const min_dim: u32 = int.min(m, n);

        try blas.copy(
            types.scast(i32, min_dim),
            a.data,
            a.inc,
            result.data,
            result.inc,
            ctx,
        );

        try blas.trmv(
            types.orderOf(B),
            types.uploOf(B),
            .trans,
            types.diagOf(B),
            types.scast(i32, min_dim),
            b.data,
            types.scast(i32, b.ld),
            result.data,
            result.inc,
            ctx,
        );

        if (comptime types.uploOf(B) == .upper) {
            if (n > min_dim) {
                try blas.gemv(
                    types.orderOf(B),
                    .trans,
                    types.scast(i32, m),
                    types.scast(i32, n - min_dim),
                    1,
                    b.data +
                        if (comptime types.orderOf(B) == .col_major)
                            min_dim * b.ld
                        else
                            min_dim,
                    types.scast(i32, b.ld),
                    a.data,
                    a.inc,
                    0,
                    result.data + min_dim * types.scast(u32, result.inc),
                    result.inc,
                    ctx,
                );
            }
        } else {
            if (m > min_dim) {
                try blas.gemv(
                    types.orderOf(B),
                    .trans,
                    types.scast(i32, m - min_dim),
                    types.scast(i32, n),
                    1,
                    b.data +
                        if (comptime types.orderOf(B) == .col_major)
                            min_dim
                        else
                            min_dim * b.ld,
                    types.scast(i32, b.ld),
                    a.data +
                        if (a.inc > 0)
                            min_dim * types.scast(u32, a.inc)
                        else
                            0,
                    a.inc,
                    1,
                    result.data,
                    result.inc,
                    ctx,
                );
            } else if (n > min_dim) {
                try blas.scal(
                    types.scast(i32, n - min_dim),
                    0,
                    result.data + min_dim * types.scast(u32, result.inc),
                    result.inc,
                    ctx,
                );
            }
        }
    }

    return result;
}
