const std = @import("std");

const types = @import("../../types.zig");
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const vecops = @import("../ops.zig");

pub inline fn add(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(Numeric(X), Numeric(Y));

    comptime if (!types.isVector(@TypeOf(x)) or !types.isVector(@TypeOf(y)))
        @compileError("Both arguments to add must be vector types");

    if (comptime types.isArbitraryPrecision(C)) {
        comptime types.validateContext(
            @TypeOf(ctx),
            .{ .element_allocator = .{ .type = std.mem.Allocator, .required = true } },
        );

        return vecops.apply2(
            allocator,
            x,
            y,
            ops.add,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        if (comptime types.numericType(C) == .int) {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            comptime types.validateContext(@TypeOf(ctx), .{});
        }

        return vecops.apply2(
            allocator,
            x,
            y,
            ops.add,
            ctx,
        );
    }
}
