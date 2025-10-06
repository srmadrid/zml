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

    var result: matrix.Banded(R, types.orderOf(X)) = try .init(allocator, x.rows, x.cols, x.lower, x.upper);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(x)) == .col_major) { // orderOf(result) == orderOf(x)
        var j: u32 = 0;
        while (j < x.cols) : (j += 1) {
            var i: u32 = if (j < result.upper) 0 else j - result.upper;
            while (i < int.min(j, x.rows)) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[(result.upper + i - j) + j * result.ld] = op(x.data[(x.upper + i - j) + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[(result.upper + i - j) + j * result.ld] = try op(x.data[(x.upper + i - j) + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                }
            }

            if (j < int.min(x.rows, x.cols)) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[result.upper + j * result.ld] = op(x.data[x.upper + j * x.ld], y.data[j]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[result.upper + j * result.ld] = try op(x.data[x.upper + j * x.ld], y.data[j], ctx);
                }
            }

            i = j + 1;
            while (i <= int.min(x.rows - 1, j + result.lower)) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[(result.upper + i - j) + j * result.ld] = op(x.data[(x.upper + i - j) + j * x.ld], constants.zero(Y, ctx) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[(result.upper + i - j) + j * result.ld] = try op(x.data[(x.upper + i - j) + j * x.ld], constants.zero(Y, ctx) catch unreachable, ctx);
                }
            }
        }
    } else {
        var i: u32 = 0;
        while (i < x.rows) : (i += 1) {
            var j: u32 = if (i < result.lower) 0 else i - result.lower;
            while (j < int.min(i, x.cols)) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + (result.lower + j - i)] = op(x.data[i * x.ld + (x.lower + j - i)], constants.zero(Y, ctx) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + (result.lower + j - i)] = try op(x.data[i * x.ld + (x.lower + j - i)], constants.zero(Y, ctx) catch unreachable, ctx);
                }
            }

            if (i < int.min(x.rows, x.cols)) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + result.lower] = op(x.data[i * x.ld + x.lower], y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + result.lower] = try op(x.data[i * x.ld + x.lower], y.data[i], ctx);
                }
            }

            j = i + 1;
            while (j <= int.min(x.cols - 1, i + result.upper)) : (j += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i * result.ld + (result.lower + j - i)] = op(x.data[i * x.ld + (x.lower + j - i)], constants.zero(Y, ctx) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i * result.ld + (result.lower + j - i)] = try op(x.data[i * x.ld + (x.lower + j - i)], constants.zero(Y, ctx) catch unreachable, ctx);
                }
            }
        }
    }

    return result;
}
