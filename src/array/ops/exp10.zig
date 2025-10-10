const std = @import("std");

const types = @import("../../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");

const arrops = @import("../ops.zig");

///
pub inline fn exp10(
    allocator: std.mem.Allocator,
    x: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X))
        @compileError("zml.array.exp10 requires x to be an array, got " ++ @typeName(X));

    if (comptime types.isArbitraryPrecision(X)) {
        comptime types.validateContext(
            @TypeOf(ctx),
            .{
                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
            },
        );

        return arrops.apply1(
            allocator,
            x,
            ops.exp10,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        comptime types.validateContext(@TypeOf(ctx), .{});

        return arrops.apply1(
            allocator,
            x,
            ops.exp10,
            ctx,
        );
    }
}
