const std = @import("std");

const types = @import("../../types.zig");
const Coerce = types.Coerce;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const arrops = @import("../ops.zig");

///
pub inline fn min(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y))
        @compileError("zml.array.min requires at least one array argument");

    const C: type = types.Coerce(types.Numeric(X), types.Numeric(Y));

    if (comptime types.isArbitraryPrecision(C)) {
        comptime if (types.numericType(C) == .rational and (types.numericType(X) == .float or types.numericType(Y) == .float)) {
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .element_allocator = .{ .type = ?std.mem.Allocator, .required = false },
                },
            );
        } else {
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        };

        return arrops.apply2(
            allocator,
            x,
            y,
            ops.min,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        comptime types.validateContext(@TypeOf(ctx), .{});

        return arrops.apply2(
            allocator,
            x,
            y,
            ops.min,
            ctx,
        );
    }
}
