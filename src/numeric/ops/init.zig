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

const constants = @import("../../constants.zig");

pub inline fn init(comptime N: type, ctx: anytype) !N {
    comptime if (!types.isNumeric(N))
        @compileError("zml.numeric.init: N must be a numeric type, got \n\tN: " ++ @typeName(N) ++ "\n");

    switch (comptime types.numericType(N)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return false;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0.0;
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return .zero;
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return constants.zero(N, .{}) catch unreachable;
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

            return .init(ctx.allocator, 0);
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

            return .init(ctx.allocator, 0, 0);
        },
        .real => @compileError("zml.numeric.init: not implemented for " ++ @typeName(N) ++ " yet."),
        .complex => @compileError("zml.numeric.init: not implemented for " ++ @typeName(N) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(N)) {
                comptime if (!types.hasMethod(N, "init", fn (std.mem.Allocator) anyerror!N, &.{std.mem.Allocator}))
                    @compileError("zml.numeric.init: " ++ @typeName(N) ++ " must implement `fn init(std.mem.Allocator) !" ++ @typeName(N) ++ "`");

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

                return N.init(ctx.allocator);
            } else {
                comptime if (!types.hasMethod(N, "init", fn () N, &.{}))
                    @compileError("zml.numeric.init: " ++ @typeName(N) ++ " must implement `fn init() " ++ @typeName(N) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return N.init();
            }
        },
    }
}
