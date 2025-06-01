const cfloat = @import("../cfloat.zig");
const types = @import("../types.zig");
const Coerce = types.Coerce;
const scast = types.scast;

pub fn pow(x: anytype, y: anytype) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.pow requires at least one of x or y to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    const l: C = scast(C, x);
    const r: C = scast(C, y);

    return cfloat.exp(r.mul(cfloat.log(l)));
}
