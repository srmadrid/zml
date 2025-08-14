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

fn loop1_(
    result: anytype,
    x: anytype,
    comptime op_: anytype,
    depth: u32,
    order: types.IterationOrder,
    ir: u32,
    ix: i32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op_));
        const idx: u32 = if (order == .left_to_right) 0 else result.ndim - 1;

        var jr: u32 = ir;
        var jx: i32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                op_(&result.data[jr], x.data[types.scast(u32, jx)]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                try op_(&result.data[jr], x.data[types.scast(u32, jx)], ctx);
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    } else {
        const idx: u32 = if (order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: u32 = ir;
        var jx: i32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            try loop1_(
                result,
                x,
                op_,
                depth - 1,
                order,
                jr,
                jx,
                ctx,
            );

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    }
}

pub fn apply1_(
    o: anytype,
    x: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    const X: type = Numeric(@TypeOf(x));

    var xx: Strided(X) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
        // Same shape
        xx = x;
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });

        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim]))
            return array.Error.NotBroadcastable;

        xx = try @constCast(&x).broadcast(bct.shape[0..bct.ndim]);
    }

    try loop1_(
        o,
        &xx,
        op_,
        o.ndim - 1,
        o.flags.order.toIterationOrder(),
        0,
        types.scast(i32, xx.offset),
        ctx,
    );

    return;
}
