const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

/// Sets the value of `o` to `x`.
///
/// When `o` is of an arbitrary precision type, its already allocated memory is
/// used, evading a new allocation (reallocation may be needed if more space is
/// needed). For `o` a fixed precision type, this is equivalent to `cast`.
pub inline fn set(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.set requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isNumeric(O))
        @compileError("zml.set requires the output to be a pointer to a numeric type, got " ++ @typeName(@TypeOf(o)));

    comptime if (!types.isNumeric(X))
        @compileError("zml.set requires the input to be a numeric type, got " ++ @typeName(X));

    switch (comptime types.numericType(O)) {
        .bool => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x != 0;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x != 0.0;
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x.re != 0.0 or x.im != 0.0;
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .int => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = if (x) 1 else 0;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = @intCast(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = @intFromFloat(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = @intFromFloat(x.re);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .float => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = if (x) 1.0 else 0.0;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = @floatFromInt(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = @floatCast(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = @floatCast(x.re);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .cfloat => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = .{
                    .re = if (x) 1.0 else 0.0,
                    .im = 0.0,
                };
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = .{
                    .re = @floatFromInt(x),
                    .im = 0.0,
                };
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = .{
                    .re = @floatCast(x),
                    .im = 0.0,
                };
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = .{
                    .re = @floatCast(x.re),
                    .im = @floatCast(x.im),
                };
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}
