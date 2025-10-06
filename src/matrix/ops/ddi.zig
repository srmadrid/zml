const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = ReturnType2(op, Numeric(X), Numeric(Y));

    if (comptime !types.isDiagonalMatrix(@TypeOf(x))) {
        var result: matrix.Diagonal(R) = try .init(allocator, y.rows, y.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < int.min(result.rows, result.cols)) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x, y.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x, y.data[i], ctx);
            }
        }

        return result;
    } else if (comptime !types.isDiagonalMatrix(@TypeOf(y))) {
        var result: matrix.Diagonal(R) = try .init(allocator, x.rows, x.cols);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < int.min(result.rows, result.cols)) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[i], y);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[i], y, ctx);
            }
        }

        return result;
    }

    var result: matrix.Diagonal(R) = try .init(allocator, x.rows, x.cols);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    var i: u32 = 0;
    while (i < int.min(result.rows, result.cols)) : (i += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            result.data[i] = op(x.data[i], y.data[i]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            result.data[i] = try op(x.data[i], y.data[i], ctx);
        }
    }

    return result;
}
