const std = @import("std");

const types = @import("../types.zig");
const Scalar = types.Scalar;
const ReturnType1 = types.ReturnType1;
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;

const dense = @import("dense.zig");
const strided = @import("strided.zig");

pub fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(ReturnType1(op, Scalar(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isSlice(X))
        @compileError("apply1: x must be an array or slice, got " ++ @typeName(X));

    switch (x.flags.storage) {
        .dense => return dense.apply1(allocator, Scalar(X), x, op, options.writeable),
        .strided => return strided.apply1(allocator, Scalar(X), x, op, options.writeable),
    }
}

// Only for in-pleace operations
pub fn apply1_(
    o: anytype,
    comptime op_: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply1: o must be an array or slice, got " ++ @typeName(O));

    if (o.flags.writeable == false)
        return error.ArrayNotWriteable;

    switch (o.flags.storage) {
        .dense => return dense.apply1_(Scalar(O), o, op_, options.allocator),
        .strided => return strided.apply1_(Scalar(O), o, op_, options.allocator),
    }
}

pub inline fn abs(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(Scalar(Scalar(@TypeOf(x)))) {
    return apply1(allocator, x, ops.abs, .{ .writeable = options.writeable });
}

pub inline fn abs_(
    o: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, ops.abs_, .{ .allocator = options.allocator });
}

pub inline fn ceil(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(Scalar(@TypeOf(x))) {
    return apply1(allocator, x, ops.ceil, .{ .writeable = options.writeable });
}

pub inline fn ceil_(
    o: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, ops.ceil_, .{ .allocator = options.allocator });
}
