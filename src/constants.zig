const std = @import("std");
const types = @import("types.zig");
const cast = types.cast;
const Scalar = types.Scalar;
const Coerce = types.Coerce;
const CoerceToArray = types.CoerceToArray;
const Child = types.Child;
const needsAllocator = types.needsAllocator;

const int = @import("int.zig");
const float = @import("float.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");
const complex = @import("complex.zig");

pub inline fn pi(
    comptime T: type,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !T {
    _ = options;
    switch (types.numericType(T)) {
        .bool => @compileError("zml.pi not defined for " ++ @typeName(T)),
        .int => @compileError("zml.pi not defined for " ++ @typeName(T)),
        .float => return float.pi(T),
        .cfloat => return .{
            .re = float.pi(T.re),
            .im = 0,
        },
        else => @compileError("zml.pi not implemented for " ++ @typeName(T) ++ " yet"),
    }
}
