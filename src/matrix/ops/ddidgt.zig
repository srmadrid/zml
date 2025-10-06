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

    var result: matrix.Tridiagonal(R) = try .init(allocator, x.rows);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    var rdata: [*]Numeric(R) = result.data;
    var ydata: [*]X = y.data + y.offset + y.osize + (y.osize - 1) - y.sdoffset;

    // Subdiagonal
    var j: u32 = 0;
    while (j < x.rows - 1) : (j += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            rdata[j] = op(constants.zero(X, ctx) catch unreachable, ydata[j]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            rdata[j] = try op(constants.zero(X, ctx) catch unreachable, ydata[j], ctx);
        }
    }

    rdata += result.size - 1;
    ydata = y.data + y.offset + y.osize - 1;

    // Diagonal
    j = 0;
    while (j < x.rows) : (j += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            rdata[j] = op(x.data[j], ydata[j]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            rdata[j] = try op(x.data[j], ydata[j], ctx);
        }
    }

    rdata += result.size;
    ydata = y.data + y.offset + y.sdoffset;

    // Superdiagonal
    j = 0;
    while (j < x.rows - 1) : (j += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            rdata[j] = op(constants.zero(Y, ctx) catch unreachable, ydata[j]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            rdata[j] = try op(constants.zero(Y, ctx) catch unreachable, ydata[j], ctx);
        }
    }

    return result;
}
