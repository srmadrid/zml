const cfloat = @import("../cfloat.zig");
const types = @import("../types.zig");
const Coerce = types.Coerce;
const scast = types.scast;

pub fn pow(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or
        (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.pow requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    const l: C = scast(C, left);
    const r: C = scast(C, right);

    return cfloat.exp(r.mul(cfloat.log(l)));
}
