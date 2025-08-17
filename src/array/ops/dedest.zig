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

fn loop2_(
    result: anytype,
    x: anytype,
    y: anytype,
    comptime op_: anytype,
    depth: u32,
    order: types.IterationOrder,
    ir: u32,
    ix: u32,
    iy: i32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op_));
        const idx: u32 = if (order == .left_to_right) 0 else result.ndim - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
        var jy: i32 = iy;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime opinfo.@"fn".params.len == 3) {
                op_(&result.data[jr], x.data[jx], y.data[types.scast(u32, jy)]);
            } else if (comptime opinfo.@"fn".params.len == 4) {
                try op_(&result.data[jr], x.data[jx], y.data[types.scast(u32, jy)], ctx);
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
            jy += y.strides[idx];
        }
    } else {
        const idx: u32 = if (order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
        var jy: i32 = iy;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            try loop2_(
                result,
                x,
                y,
                op_,
                depth - 1,
                order,
                jr,
                jx,
                jy,
                ctx,
            );

            jr += result.strides[idx];
            jx += x.strides[idx];
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
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));

    var xx: Dense(X) = undefined;
    var yy: Strided(Y) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
        std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim]))
    {
        xx = x;
        yy = y;
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim], y.shape[0..y.ndim] });

        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim]))
            return array.Error.NotBroadcastable;

        xx = try @constCast(&x).broadcast(bct.shape[0..bct.ndim]);
        yy = try @constCast(&y).broadcast(bct.shape[0..bct.ndim]);
    }

    try loop2_(
        o,
        &xx,
        &yy,
        op_,
        o.ndim - 1,
        o.flags.order.resolve2(xx.flags.order).toIterationOrder(),
        0,
        0,
        types.scast(i32, yy.offset),
        ctx,
    );

    return;
}
