const types = @import("../types.zig");
const cfloat = @import("../cfloat.zig");
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");

pub fn Pow(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.cfloat) or !types.numericType(Y).le(.cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("zml.cfloat.pow: at least one of x or y to be a cfloat, the other must be a bool, an int, a float or a cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub fn pow(x: anytype, y: anytype) Pow(@TypeOf(x), @TypeOf(y)) {
    const R: type = Pow(@TypeOf(x), @TypeOf(y));

    const xx: R = types.scast(R, x);
    const yy: R = types.scast(R, y);

    const absx: @TypeOf(xx.re) = cfloat.abs(xx);
    if (ops.eq(absx, 0, .{}) catch unreachable)
        return .{
            .re = constants.zero(@TypeOf(xx.re), .{}) catch unreachable,
            .im = constants.zero(@TypeOf(xx.im), .{}) catch unreachable,
        };

    const argx: @TypeOf(xx.re) = cfloat.arg(xx);
    var r: @TypeOf(xx.re) = ops.pow(
        absx,
        yy.re,
        .{},
    ) catch unreachable;
    var theta: @TypeOf(xx.re) = ops.mul(
        yy.re,
        argx,
        .{},
    ) catch unreachable;
    if (ops.ne(yy.im, 0, .{}) catch unreachable) {
        r = ops.mul(
            r,
            ops.exp(
                ops.mul(
                    ops.neg(yy.im, .{}) catch unreachable,
                    argx,
                    .{},
                ) catch unreachable,
                .{},
            ) catch unreachable,
            .{},
        ) catch unreachable;

        theta = ops.add(
            theta,
            ops.mul(
                yy.im,
                ops.log(absx, .{}) catch unreachable,
                .{},
            ) catch unreachable,
            .{},
        ) catch unreachable;
    }

    return .{
        .re = ops.mul(
            r,
            ops.cos(theta, .{}) catch unreachable,
            .{},
        ) catch unreachable,
        .im = ops.mul(
            r,
            ops.sin(theta, .{}) catch unreachable,
            .{},
        ) catch unreachable,
    };
}
