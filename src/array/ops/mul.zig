const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const arrops = @import("../ops.zig");

///
pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !MulCoerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y))
        @compileError("zml.array.mul requires at least one array argument");

    const C: type = types.MulCoerce(types.Numeric(X), types.Numeric(Y));

    if (comptime types.isArbitraryPrecision(C)) {
        comptime types.validateContext(
            @TypeOf(ctx),
            .{
                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
            },
        );

        return arrops.apply2(
            allocator,
            x,
            y,
            ops.mul,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        if (comptime types.numericType(C) == .int) {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .mode = .{ .type = int.Mode, .required = false },
                },
            );
        } else {
            comptime types.validateContext(@TypeOf(ctx), .{});
        }

        return arrops.apply2(
            allocator,
            x,
            y,
            ops.mul,
            ctx,
        );
    }
}
