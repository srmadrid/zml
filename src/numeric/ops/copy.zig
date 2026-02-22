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

pub inline fn copy(x: anytype, ctx: anytype) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.copy: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .dyadic, .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .integer => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the integer's memory allocation.",
                    },
                },
            );

            return x.copy(ctx.allocator);
        },
        .rational => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the rational's memory allocation.",
                    },
                },
            );

            return x.copy(ctx.allocator);
        },
        .real => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the real's memory allocation.",
                    },
                },
            );

            return x.copy(ctx.allocator);
        },
        .complex => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the complex's memory allocation.",
                    },
                },
            );

            return x.copy(ctx.allocator);
        },
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "copy", fn (X, std.mem.Allocator) anyerror!X, &.{ X, std.mem.Allocator }))
                    @compileError("zml.numeric.copy: " ++ @typeName(X) ++ " must implement `fn copy(" ++ @typeName(X) ++ ", std.mem.Allocator) !" ++ @typeName(X) ++ "`");

                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the custom numeric's memory allocation.",
                        },
                    },
                );

                return x.copy(ctx.allocator);
            } else {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return x;
            }
        },
    }
}
