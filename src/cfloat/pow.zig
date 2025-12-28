const types = @import("../types.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const Coerce = types.Coerce;

pub fn pow(x: anytype, y: anytype) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.pow requires at least one of x or y to be a cfloat, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    const xx: C = types.scast(C, x);
    const yy: C = types.scast(C, y);

    const absx: @TypeOf(xx.re) = cfloat.abs(xx);
    if (absx == 0.0)
        return .{
            .re = 0.0,
            .im = 0.0,
        };

    const argx: @TypeOf(xx.re) = cfloat.arg(xx);
    var r: @TypeOf(xx.re) = float.pow(absx, yy.re);
    var theta: @TypeOf(xx.re) = yy.re * argx;
    if (yy.im != 0.0) {
        r *= float.exp(-yy.im * argx);
        theta += yy.im * float.log(absx);
    }

    return .{
        .re = r * float.cos(theta),
        .im = r * float.sin(theta),
    };
}
