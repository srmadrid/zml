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

pub inline fn tocompleteabs_(o: anytype, x: anytype, ctx: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O) or
        !types.isNumeric(types.Child(O)) or
        !types.isNumeric(X))
        @compileError("zml.numeric.abs_: o must be a mutable one-itme pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = types.Child(O);
    const R: type = numeric.Abs(X);

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .dyadic, .cfloat => switch (comptime types.numericType(X)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational, .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    numeric.abs(x, ctx) catch unreachable,
                    ctx,
                ) catch unreachable;
            },
            .complex => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
            .custom => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
        },
        .integer => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.abs_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    int.abs(x),
                    ctx,
                );
            },
            .float => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    float.abs(x),
                    ctx,
                );
            },
            .cfloat => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    cfloat.abs(x),
                    ctx,
                );
            },
            .integer => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (o.limbs == x.limbs) {
                    o.positive = true;
                } else {
                    if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                        try ops.set(
                            o,
                            integer.abs(null, x) catch unreachable,
                            .{ .allocator = allocator },
                        );
                    } else {
                        return error.AllocatorRequired;
                    }
                }
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    rational.abs(null, x) catch unreachable,
                    ctx,
                );
            },
            .real => @compileError("zml.abs_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .complex => @compileError("zml.abs_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .rational => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.abs_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    int.abs(x),
                    ctx,
                );
            },
            .float => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    float.abs(x),
                    ctx,
                );
            },
            .cfloat => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    cfloat.abs(x),
                    ctx,
                );
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                try ops.set(
                    o,
                    integer.abs(null, x) catch unreachable,
                    ctx,
                );
            },
            .rational => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (o.num.limbs == x.num.limbs and
                    o.den.limbs == x.den.limbs)
                {
                    o.num.positive = true;
                } else {
                    if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                        try ops.set(
                            o,
                            rational.abs(null, x) catch unreachable,
                            .{ .allocator = allocator },
                        );
                    } else {
                        return error.AllocatorRequired;
                    }
                }
            },
            .real => @compileError("zml.abs_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .complex => @compileError("zml.abs_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .real => @compileError("zml.abs_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .complex => @compileError("zml.abs_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
