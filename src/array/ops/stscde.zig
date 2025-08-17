const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;

const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const array = @import("../../array.zig");
const Order = types.Order;
const Flags = array.Flags;
const Range = array.Range;

const Dense = array.Dense;
const Strided = array.Strided;

fn loop2_left_(
    result: anytype,
    x: anytype,
    y: anytype,
    comptime op_: anytype,
    depth: u32,
    order: types.IterationOrder,
    ir: i32,
    iy: u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op_));
        const idx: u32 = if (order == .left_to_right) 0 else result.ndim - 1;

        var jr: i32 = ir;
        var jy: u32 = iy;
        var j: u32 = 0;
        while (j < y.shape[idx]) : (j += 1) {
            if (comptime opinfo.@"fn".params.len == 3) {
                op_(&result.data[types.scast(u32, jr)], x, y.data[jy]);
            } else if (comptime opinfo.@"fn".params.len == 4) {
                try op_(&result.data[types.scast(u32, jr)], x, y.data[jy], ctx);
            }

            jr += result.strides[idx];
            jy += y.strides[idx];
        }
    } else {
        const idx: u32 = if (order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: i32 = ir;
        var jy: u32 = iy;
        var j: u32 = 0;
        while (j < y.shape[idx]) : (j += 1) {
            try loop2_left_(
                result,
                x,
                y,
                op_,
                depth - 1,
                order,
                jr,
                jy,
                ctx,
            );

            jr += result.strides[idx];
            jy += y.strides[idx];
        }
    }
}

pub fn apply2_(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    const Y: type = Numeric(@TypeOf(y));

    var yy: Dense(Y) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim])) {
        yy = y;
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], y.shape[0..y.ndim] });
        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
            return array.Error.NotBroadcastable;
        }

        yy = try @constCast(&y).broadcast(bct.shape[0..bct.ndim]);
    }

    try loop2_left_(
        o,
        x,
        &yy,
        op_,
        o.ndim - 1,
        y.flags.order.toIterationOrder(),
        types.scast(i32, o.offset),
        0,
        ctx,
    );

    return;
}
