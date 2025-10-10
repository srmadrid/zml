const std = @import("std");

const types = @import("../../types.zig");
const EnsureArray = types.EnsureArray;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");

const arrops = @import("../ops.zig");

///
pub inline fn neg(
    allocator: std.mem.Allocator,
    x: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Numeric(@TypeOf(x))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X))
        @compileError("zml.array.neg requires x to be an array, got " ++ @typeName(X));

    if (comptime types.isArbitraryPrecision(X)) {
        comptime types.validateContext(
            @TypeOf(ctx),
            .{
                .element_allocator = .{ .type = ?std.mem.Allocator, .required = false },
            },
        );

        return arrops.apply1(
            allocator,
            x,
            ops.neg,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        comptime types.validateContext(@TypeOf(ctx), .{});

        return arrops.apply1(
            allocator,
            x,
            ops.neg,
            ctx,
        );
    }
}
