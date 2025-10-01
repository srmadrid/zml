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

    var result: R = try .init(allocator, x.rows, x.cols);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    var j: u32 = 0;
    while (j < x.cols) : (j += 1) {
        var i: u32 = 0;
        while (i < j) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.at(i, j));
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.at(i, j), ctx);
            }
        }

        if (comptime opinfo.@"fn".params.len == 2) {
            result.data[j + j * result.ld] = op(x.data[j], y.at(j, j));
        } else if (comptime opinfo.@"fn".params.len == 3) {
            result.data[j + j * result.ld] = try op(x.data[j], y.at(j, j), ctx);
        }

        i = j + 1;
        while (i < x.rows) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i + j * result.ld] = op(constants.zero(X, ctx) catch unreachable, y.at(i, j));
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i + j * result.ld] = try op(constants.zero(X, ctx) catch unreachable, y.at(i, j), ctx);
            }
        }
    }

    return result;
}
