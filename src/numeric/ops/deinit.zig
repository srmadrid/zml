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

pub inline fn deinit(x: anytype, ctx: anytype) void {
    comptime var X: type = @TypeOf(x);

    comptime if (!types.isPointer(X) or types.isConstPointer(X) or
        !types.isNumeric(types.Child(X)))
        @compileError("zml.numeric.deinit: x must be a mutable one-itme pointer to a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    X = types.Child(X);

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .dyadic, .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            // No deinitialization needed for fixed precision types, this is a no-op.
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

            x.deinit(ctx.allocator);
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

            x.deinit(ctx.allocator);
        },
        .real => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the real's memory deallocation. Must be the same allocator used to initialize it.",
                    },
                },
            );

            x.deinit(ctx.allocator);
        },
        .complex => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the complex's memory deallocation. Must be the same allocator used to initialize it.",
                    },
                },
            );

            x.deinit(ctx.allocator);
        },
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "deinit", fn (*X, std.mem.Allocator) void, &.{}))
                    @compileError("zml.deinit: " ++ @typeName(X) ++ " must implement `fn deinit(*" ++ @typeName(X) ++ ", std.mem.Allocator) void`");

                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the custom numeric's memory deallocation. Must be the same allocator used to initialize it.",
                        },
                    },
                );

                x.deinit(ctx.allocator);
            } else {
                comptime types.validateContext(@TypeOf(ctx), .{});

                // No deinitialization needed for non-allocated custom types, this is a no-op.
            }
        },
    }
}
