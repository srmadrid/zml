const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn sin(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or !types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.sin: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    var s: @TypeOf(z.re) = undefined;
    var c: @TypeOf(z.re) = undefined;
    if (ops.le(ops.abs(z.re, .{}) catch unreachable, 0.5, .{}) catch unreachable) {
        s = ops.sinh(z.im, .{}) catch unreachable;
        c = ops.cosh(z.im, .{}) catch unreachable;
    } else {
        var e: @TypeOf(z.re) = ops.exp(z.im, .{}) catch unreachable;
        const e_inv: @TypeOf(z.re) = ops.div(
            0.5,
            e,
            .{},
        ) catch unreachable;
        e = ops.mul(0.5, e, .{}) catch unreachable;
        s = ops.sub(e, e_inv, .{}) catch unreachable;
        c = ops.add(e, e_inv, .{}) catch unreachable;
    }

    return .{
        .re = ops.mul(ops.sin(z.re, .{}) catch unreachable, c, .{}) catch unreachable,
        .im = ops.mul(ops.cos(z.re, .{}) catch unreachable, s, .{}) catch unreachable,
    };
}
