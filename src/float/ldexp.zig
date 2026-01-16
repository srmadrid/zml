const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");

pub fn Ldexp(comptime X: type) type {
    comptime if (!types.isNumeric(X) or !types.numericType(X).le(.float))
        @compileError("zml.float.ldexp: x must be a bool, an int or a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return types.EnsureFloat(X);
}

pub inline fn ldexp(x: anytype, n: i32) Ldexp(@TypeOf(x)) {
    const xx: Ldexp(@TypeOf(x)) = types.scast(Ldexp(@TypeOf(x)), x);

    if (!std.math.isFinite(xx) or xx == 0)
        return xx + xx;

    return float.scalbn(xx, n);
}
