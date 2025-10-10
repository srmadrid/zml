const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const array = @import("../array.zig");

///
pub inline fn sinh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sinh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.sinh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domainType(O)) {
        .array => switch (comptime types.domainType(X)) {
            .array, .numeric => { // array = sinh(numeric), array = sinh(array)
                comptime if (types.isArbitraryPrecision(types.Numeric(O))) {
                    // To be though about: when O == X it is trivial. If O != X, we
                    // must reason about wether a view of type O can be created from
                    // X.
                } else {
                    types.validateContext(@TypeOf(ctx), .{});
                };

                return array.sinh_(
                    o,
                    x,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => unreachable,
        },
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => { // numeric = sinh(numeric)
                switch (comptime types.numericType(X)) {
                    .bool => @compileError("zml.sinh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    .int => {
                        comptime if (@typeInfo(X).int.signedness == .unsigned)
                            @compileError("zml.sinh_ not defined for unsigned integers, got " ++ @typeName(X));

                        comptime if (types.isArbitraryPrecision(O)) {
                            // To be though about
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        } else {
                            types.validateContext(@TypeOf(ctx), .{});
                        };

                        try ops.set(
                            o,
                            float.sinh(x),
                            ctx,
                        );
                    },
                    .float => {
                        comptime if (types.isArbitraryPrecision(O)) {
                            // To be though about
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        } else {
                            types.validateContext(@TypeOf(ctx), .{});
                        };

                        try ops.set(
                            o,
                            float.sinh(x),
                            ctx,
                        );
                    },
                    .cfloat => {
                        comptime if (types.isArbitraryPrecision(O)) {
                            // To be though about
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        } else {
                            types.validateContext(@TypeOf(ctx), .{});
                        };

                        try ops.set(
                            o,
                            cfloat.sinh(x),
                            ctx,
                        );
                    },
                    else => @compileError("zml.sinh_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.sinh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.sinh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
