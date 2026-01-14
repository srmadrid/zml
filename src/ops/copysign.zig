const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

///
pub inline fn copysign(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Numeric(@TypeOf(x))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.copysign not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domain(X)) {
        .array => switch (comptime types.domain(Y)) {
            .array, .numeric => { // copysign(array, array), copysign(array, numeric)
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
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

                return array.copysign(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.copysign not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domain(Y)) {
            .numeric => { // copysign(numeric, numeric)
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.copysign not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.copysign(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.copysign(x, y);
                    },
                    .cfloat => @compileError("zml.copysign not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .complex => @compileError("zml.copysign not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    else => @compileError("zml.copysign between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.copysign not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.copysign not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
