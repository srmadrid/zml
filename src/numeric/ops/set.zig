const std = @import("std");

const types = @import("../../types.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const cfloat = @import("../../cfloat.zig");
const integer = @import("../../integer.zig");
const rational = @import("../../rational.zig");
const real = @import("../../real.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

/// Sets the value of `o` to `x`.
///
/// When `o` is of an arbitrary precision type, its already allocated memory is
/// used, evading a new allocation (reallocation may be needed if more space is
/// needed). For `o` a fixed precision type, this is equivalent to `cast`.
pub inline fn set(o: anytype, x: anytype, ctx: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O) or
        !types.isNumeric(types.Child(O)) or
        !types.isNumeric(X))
        @compileError("zml.numeric.set: o must be a mutable one-itme pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = types.Child(O);

    switch (comptime types.numericType(O)) {
        .bool => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .custom => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
        },
        .int => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
        },
        .float => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
        },
        .dyadic => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
        },
        .cfloat => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
        },
        .integer => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the integer's memory deallocation. Must be the same allocator used to initialize it.",
                    },
                },
            );

            try o.set(ctx.allocator, x);
        },
        .rational => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the rational's memory deallocation. Must be the same allocator used to initialize it.",
                    },
                },
            );

            try o.set(ctx.allocator, x, 1);
        },
        .real => @compileError("zml.numeric.set: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.set: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(O)) {
                comptime if (!types.hasMethod(O, "set", fn (*O, std.mem.Allocator, X) anyerror!void, &.{ *O, std.mem.Allocator, X }))
                    @compileError("zml.numeric.set: " ++ @typeName(O) ++ " must implement `fn set(*" ++ @typeName(O) ++ ", std.mem.Allocator, " ++ @typeName(X) ++ ") !void`");

                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the custom numeric's memory allocation. Must be the same allocator used to initialize it.",
                        },
                    },
                );

                try o.set(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(O, "set", fn (*O, X) void, &.{ *O, X }))
                    @compileError("zml.numeric.set: " ++ @typeName(O) ++ " must implement `fn set(*" ++ @typeName(O) ++ ", " ++ @typeName(X) ++ ") void`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                o.set(x);
            }
        },
    }
}
