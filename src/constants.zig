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
            const spec =
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                };

            comptime validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
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
        .rational => {
            const spec =
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                };

            comptime validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var num: integer.Integer = try .init(allocator, 2);
                errdefer num.deinit(allocator);

                var den: integer.Integer = try .init(allocator, 2);
                den.limbs[0] = 1;
                den.size = 1;

                return .{
                    .num = num,
                    .den = den,
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .num = zero(integer.Integer, .{}) catch unreachable,
                    .den = one(integer.Integer, .{}) catch unreachable,
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
    switch (comptime types.numericType(T)) {
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
            const spec =
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                };

            comptime validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
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
        .rational => {
            const spec =
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                };

            comptime validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var num: integer.Integer = try .init(allocator, 2);
                errdefer num.deinit(allocator);
                num.limbs[0] = 1;
                num.size = 1;

                var den: integer.Integer = try .init(allocator, 2);
                den.limbs[0] = 1;
                den.size = 1;

                return .{
                    .num = num,
                    .den = den,
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .num = one(integer.Integer, .{}) catch unreachable,
                    .den = one(integer.Integer, .{}) catch unreachable,
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
    switch (comptime types.numericType(T)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return true;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 2;
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 2.0;
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return .{ .re = 2.0, .im = 0.0 };
        },
        .integer => {
            const spec =
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                };

            comptime validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var result: integer.Integer = try .init(allocator, 2);
                result.limbs[0] = 2;
                result.size = 1;
                return result;
            } else {
                return .{
                    .limbs = @constCast((&[_]u32{2}).ptr),
                    .size = 1,
                    ._llen = 0,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .rational => {
            const spec =
                .{
                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                };

            comptime validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var num: integer.Integer = try .init(allocator, 2);
                errdefer num.deinit(allocator);
                num.limbs[0] = 2;
                num.size = 1;

                var den: integer.Integer = try .init(allocator, 2);
                den.limbs[0] = 1;
                den.size = 1;

                return .{
                    .num = num,
                    .den = den,
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .num = two(integer.Integer, .{}) catch unreachable,
                    .den = one(integer.Integer, .{}) catch unreachable,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
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
