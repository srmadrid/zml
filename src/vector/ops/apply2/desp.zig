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
    if (x.inc == 1) {
        var iy: u32 = 0;
        while (i < result.len) : (i += 1) {
            if (y.idx[iy] == i) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y.data[iy]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y.data[iy], ctx);
                }

                iy += 1;
            } else {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], constants.zero(types.Numeric(Y), .{}) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], constants.zero(types.Numeric(Y), .{}) catch unreachable, ctx);
                }
            }
        }
    } else {
        var ix: i32 = if (x.inc < 0) (-types.scast(i32, x.len) + 1) * x.inc else 0;
        var iy: u32 = 0;
        while (i < result.len) : (i += 1) {
            if (y.idx[iy] == i) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[types.scast(u32, ix)], y.data[iy]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[types.scast(u32, ix)], y.data[iy], ctx);
                }

                iy += 1;
            } else {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[types.scast(u32, ix)], constants.zero(types.Numeric(Y), .{}) catch unreachable);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[types.scast(u32, ix)], constants.zero(types.Numeric(Y), .{}) catch unreachable, ctx);
                }
            }

            ix += x.inc;
        }
    }

    return result;
}
