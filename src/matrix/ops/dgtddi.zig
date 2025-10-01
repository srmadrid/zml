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

    var result: R = try .init(allocator, x.size);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    var rdata: [*]Numeric(R) = result.data;
    var xdata: [*]X = x.data + x.offset + x.osize + (x.osize - 1) - x.sdoffset;

    // Subdiagonal
    var j: u32 = 0;
    while (j < x.size - 1) : (j += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            rdata[j] = op(xdata[j], constants.zero(Y, ctx) catch unreachable);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            rdata[j] = try op(xdata[j], constants.zero(Y, ctx) catch unreachable, ctx);
        }
    }

    rdata += result.size - 1;
    xdata = x.data + x.offset + x.osize - 1;

    // Diagonal
    j = 0;
    while (j < x.size) : (j += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            rdata[j] = op(xdata[j], y.data[j]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            rdata[j] = try op(xdata[j], y.data[j], ctx);
        }
    }

    rdata += result.size;
    xdata = x.data + x.offset + x.sdoffset;

    // Superdiagonal
    j = 0;
    while (j < x.size - 1) : (j += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            rdata[j] = op(xdata[j], constants.zero(Y, ctx) catch unreachable);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            rdata[j] = try op(xdata[j], constants.zero(Y, ctx) catch unreachable, ctx);
        }
    }

    return result;
}
