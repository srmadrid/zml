const std = @import("std");

const types = @import("../../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureVector = types.EnsureVector;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const ops = @import("../../ops.zig");
const constants = @import("../../../constants.zig");
const int = @import("../../../int.zig");

const vector = @import("../../../vector.zig");
const Dense = vector.Dense;

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !Dense(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = ReturnType2(op, Numeric(X), Numeric(Y));

    if (x.len != y.len)
        return vector.Error.DimensionMismatch;

    var result: Dense(R) = try .init(allocator, x.len);
    errdefer result.deinit(allocator);

    var i: u32 = 0;

    errdefer result._cleanup(
        i,
        types.renameStructFields(
            types.keepStructFields(
                ctx,
                &.{"allocator"},
            ),
            .{ .allocator = "element_allocator" },
        ),
    );

    const opinfo = @typeInfo(@TypeOf(op));
    if (y.inc == 1) {
        var ix: u32 = 0;
        while (i < result.len) : (i += 1) {
            if (x.idx[ix] == i) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[ix], y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[ix], y.data[i], ctx);
                }

                ix += 1;
            } else {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(constants.zero(types.Numeric(X), .{}) catch unreachable, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(constants.zero(types.Numeric(X), .{}) catch unreachable, y.data[i], ctx);
                }
            }
        }
    } else {
        var ix: u32 = 0;
        var iy: i32 = if (y.inc < 0) (-types.scast(i32, y.len) + 1) * y.inc else 0;
        while (i < result.len) : (i += 1) {
            if (x.idx[ix] == i) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[ix], y.data[types.scast(u32, iy)]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[ix], y.data[types.scast(u32, iy)], ctx);
                }

                ix += 1;
            } else {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(constants.zero(types.Numeric(X), .{}) catch unreachable, y.data[types.scast(u32, iy)]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(constants.zero(types.Numeric(X), .{}) catch unreachable, y.data[types.scast(u32, iy)], ctx);
                }
            }

            iy += y.inc;
        }
    }

    return result;
}
