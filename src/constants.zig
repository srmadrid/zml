const std = @import("std");
const types = @import("types.zig");
const cast = types.cast;
const Scalar = types.Scalar;
const Coerce = types.Coerce;
const CoerceToArray = types.CoerceToArray;
const Child = types.Child;
const isArbitraryPrecision = types.isArbitraryPrecision;
const validateContext = types.validateContext;

const int = @import("int.zig");
const float = @import("float.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");
const complex = @import("complex.zig");

pub inline fn zero(
    comptime T: type,
    ctx: anytype,
) !T {
    comptime if (types.isArbitraryPrecision(T)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    switch (types.numericType(T)) {
        .bool => return false,
        .int => return 0,
        .float => return 0,
        .cfloat => return .{
            .re = 0,
            .im = 0,
        },
        else => @compileError("zml.zero not implemented for " ++ @typeName(T) ++ " yet"),
    }
}

pub inline fn one(
    comptime T: type,
    ctx: anytype,
) !T {
    comptime if (types.isArbitraryPrecision(T)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    switch (types.numericType(T)) {
        .bool => return false,
        .int => return 1,
        .float => return 1,
        .cfloat => return .{
            .re = 1,
            .im = 0,
        },
        else => @compileError("zml.zero not implemented for " ++ @typeName(T) ++ " yet"),
    }
}

pub inline fn pi(
    comptime T: type,
    ctx: anytype,
) !T {
    comptime if (types.isArbitraryPrecision(T)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

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
