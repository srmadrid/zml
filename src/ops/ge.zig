const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const Coerce = types.Coerce;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

///
pub inline fn ge(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .array, .numeric => { // array >= array, array >= numeric
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return array.ge(
                    ctx.array_allocator,
                    x,
                    y,
                );
            },
            else => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .array => { // numeric >= array
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return array.ge(
                    ctx.array_allocator,
                    x,
                    y,
                );
            },
            .numeric => { // numeric >= numeric
                switch (comptime types.numericType(C)) {
                    .bool => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return x >= y;
                    },
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.ge(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.ge(x, y);
                    },
                    .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
