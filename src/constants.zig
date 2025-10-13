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
    switch (types.numericType(T)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return false;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 0.0;
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return .{ .re = 0.0, .im = 0.0 };
        },
        .integer => {
            comptime validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                },
            );

            if (types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null)) |allocator| {
                return .init(allocator, 2);
            } else {
                return .{
                    .limbs = &.{},
                    .size = 0,
                    ._llen = 0,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        else => @compileError("zml.one not implemented for " ++ @typeName(T) ++ " yet"),
    }
}

pub inline fn one(
    comptime T: type,
    ctx: anytype,
) !T {
    switch (types.numericType(T)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return true;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 1;
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 1.0;
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return .{ .re = 1.0, .im = 0.0 };
        },
        .integer => {
            comptime validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                },
            );

            if (types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null)) |allocator| {
                var result: integer.Integer = try .init(allocator, 2);
                result.limbs[0] = 1;
                result.size = 1;
                return result;
            } else {
                return .{
                    .limbs = @constCast((&[_]u32{1}).ptr),
                    .size = 1,
                    ._llen = 0,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        else => @compileError("zml.one not implemented for " ++ @typeName(T) ++ " yet"),
    }
}

pub inline fn two(
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
        .bool => return true,
        .int => return 2,
        .float => return 2,
        .cfloat => return .{
            .re = 2,
            .im = 0,
        },
        else => @compileError("zml.two not implemented for " ++ @typeName(T) ++ " yet"),
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
