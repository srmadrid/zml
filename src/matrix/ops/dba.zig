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
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y));

    if (comptime !types.isBandedMatrix(@TypeOf(x))) {
        var result: matrix.Banded(R, types.orderOf(Y)) = try .init(allocator, y.rows, y.cols, y.lower, y.upper);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(y)) == .col_major) { // c
            var j: u32 = 0;
            while (j < result.cols) : (j += 1) {
                var i: u32 = if (j < result.upper) 0 else j - result.upper;
                while (i <= int.min(result.rows - 1, j + result.lower)) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[(result.upper + i - j) + j * result.ld] = op(x, y.data[(y.upper + i - j) + j * y.ld]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[(result.upper + i - j) + j * result.ld] = try op(x, y.data[(y.upper + i - j) + j * y.ld], ctx);
                    }
                }
            }
        } else { // r
            var i: u32 = 0;
            while (i < result.rows) : (i += 1) {
                var j: u32 = if (i < result.lower) 0 else i - result.lower;
                while (j <= int.min(result.cols - 1, i + result.upper)) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + (result.lower + j - i)] = op(x, y.data[i * y.ld + (y.lower + j - i)]);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + (result.lower + j - i)] = try op(x, y.data[i * y.ld + (y.lower + j - i)], ctx);
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isBandedMatrix(@TypeOf(y))) {
        var result: matrix.Banded(R, types.orderOf(Y)) = try .init(allocator, x.rows, x.cols, x.lower, x.upper);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) { // c
            var j: u32 = 0;
            while (j < result.cols) : (j += 1) {
                var i: u32 = if (j < result.upper) 0 else j - result.upper;
                while (i <= int.min(result.rows - 1, j + result.lower)) : (i += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[(result.upper + i - j) + j * result.ld] = op(x.data[(x.upper + i - j) + j * x.ld], y);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[(result.upper + i - j) + j * result.ld] = try op(x.data[(x.upper + i - j) + j * x.ld], y, ctx);
                    }
                }
            }
        } else { // r
            var i: u32 = 0;
            while (i < result.rows) : (i += 1) {
                var j: u32 = if (i < result.lower) 0 else i - result.lower;
                while (j <= int.min(result.cols - 1, i + result.upper)) : (j += 1) {
                    if (comptime opinfo.@"fn".params.len == 2) {
                        result.data[i * result.ld + (result.lower + j - i)] = op(x.data[i * x.ld + (x.lower + j - i)], y);
                    } else if (comptime opinfo.@"fn".params.len == 3) {
                        result.data[i * result.ld + (result.lower + j - i)] = try op(x.data[i * x.ld + (x.lower + j - i)], y, ctx);
                    }
                }
            }
        }

        return result;
    }

    var result: matrix.Banded(R, types.orderOf(X)) = try .init(allocator, x.rows, x.cols, int.max(x.lower, y.lower), int.max(x.upper, y.upper));
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
        var j: u32 = 0;
        while (j < result.cols) : (j += 1) {
            var i: u32 = if (j < result.upper) 0 else j - result.upper;
            while (i <= int.min(result.rows - 1, j + result.lower)) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[(result.upper + i - j) + j * result.ld] = op(
                        x.get(i, j) catch unreachable,
                        y.get(i, j) catch unreachable,
                    );
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[(result.upper + i - j) + j * result.ld] = try op(
                        x.get(i, j) catch unreachable,
                        y.get(i, j) catch unreachable,
                        ctx,
                    );
                }
            }
        }
    } else {
        var i: u32 = 0;
        while (i < result.rows) : (i += 1) {
            var j: u32 = if (i < result.lower) 0 else i - result.lower;
            while (j <= int.min(result.cols - 1, i + result.upper)) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + (result.lower + j - i)] = op(
                        x.get(i, j) catch unreachable,
                        y.get(i, j) catch unreachable,
                    );
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + (result.lower + j - i)] = try op(
                        x.get(i, j) catch unreachable,
                        y.get(i, j) catch unreachable,
                        ctx,
                    );
                }
            }
        }
    }

    return result;
}
