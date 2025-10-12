const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const array = @import("../array.zig");

///
pub inline fn abs(
    x: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Scalar(Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.abs not defined for " ++ @typeName(X));

    switch (comptime types.domainType(X)) {
        .array => {
            comptime if (types.isArbitraryPrecision(Numeric(X))) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .element_allocator = .{ .type = ?std.mem.Allocator, .required = false },
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

            return array.abs(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.abs not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.abs(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.abs(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.abs(x);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                return integer.abs(types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null), x);
            },
            else => @compileError("zml.abs for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
