const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");

const arrops = @import("../ops.zig");

pub inline fn gt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.array.gt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime !types.isArray(O))
        @compileError("zml.array.gt_ requires the output to be an array, got " ++ @typeName(O));

    if (comptime types.isArbitraryPrecision(types.Numeric(O))) {
        comptime types.validateContext(
            @TypeOf(ctx),
            .{
                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
            },
        );

        return arrops.apply2_(
            o,
            x,
            y,
            ops.gt_,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        comptime types.validateContext(@TypeOf(ctx), .{});

        return arrops.apply2_(
            o,
            x,
            y,
            ops.gt_,
            ctx,
        );
    }
}
