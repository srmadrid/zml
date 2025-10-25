const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const real = @import("../real.zig");
const complex = @import("../complex.zig");

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
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = integer.ne(x, 0);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = rational.ne(x, 0);
            },
            .real => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = complex.ne(x, 0);
            },
            .expression => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
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
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x.toInt(O);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x.toInt(O);
            },
            .real => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x.re.toInt(O);
            },
            .expression => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
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
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = try x.toFloat(O);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = try x.toFloat(O);
            },
            .real => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x.re.toFloat(O);
            },
            .expression => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
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
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = .{
                    .re = try x.toFloat(types.Scalar(O)),
                    .im = 0.0,
                };
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = .{
                    .re = try x.toFloat(types.Scalar(O)),
                    .im = 0.0,
                };
            },
            .real => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = .{
                    .re = x.re.toFloat(types.Scalar(O)),
                    .im = x.im.toFloat(types.Scalar(O)),
                };
            },
            .expression => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .integer => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );

            try o.set(ctx.allocator, x);
        },
        .rational => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );

            try o.set(ctx.allocator, x, 1);
        },
        .real => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        .complex => {
            if (comptime types.isComplex(X)) {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try o.set(ctx.allocator, x.re, x.im);
            } else {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try o.set(ctx.allocator, x, 0);
            }
        },
        .expression => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}
