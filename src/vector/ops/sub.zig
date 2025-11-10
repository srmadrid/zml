const std = @import("std");

const types = @import("../../types.zig");
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const vecops = @import("../ops.zig");

pub inline fn sub(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (!types.isVector(@TypeOf(x)) or !types.isVector(@TypeOf(y)))
        @compileError("Both arguments to sub must be vector types");

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .mode = .{ .type = int.Mode, .required = false, .default = .default },
                },
            );
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return vecops.apply2(
        allocator,
        x,
        y,
        ops.sub,
        ctx,
    );
}
