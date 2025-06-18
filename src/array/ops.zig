const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;

inline fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    options: struct {
        writeable: bool = true,
    },
) Array(ReturnType1(op, @TypeOf(x))) {
    comptime if (!types.isArray(x))
        @compileError("apply1: x must be an array, got " ++ @typeName(x));
}
