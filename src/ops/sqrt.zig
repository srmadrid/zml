const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");

///
pub inline fn sqrt(
    x: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.sqrt not defined for " ++ @typeName(X));

    switch (comptime types.domainType(X)) {
        .array => {
            comptime if (types.isArbitraryPrecision(Numeric(X))) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            } else {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            };

            return array.sqrt(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sqrt not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sqrt(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sqrt(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sqrt(x);
            },
            else => @compileError("zml.sqrt for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
