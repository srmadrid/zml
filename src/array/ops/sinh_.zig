const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");

const arrops = @import("../ops.zig");

///
pub inline fn sinh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.array.sinh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime !types.isArray(O))
        @compileError("zml.array.sinh_ requires the output to be an array, got " ++ @typeName(O));

    if (comptime types.isArbitraryPrecision(types.Numeric(O))) {
        // To be thought about
        types.validateContext(
            @TypeOf(ctx),
            .{
                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
            },
        );

        return arrops.apply1_(
            o,
            x,
            ops.sinh_,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        comptime if (types.isArbitraryPrecision(X)) {
            // To be thought about
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        };

        return arrops.apply1_(
            o,
            x,
            ops.sinh_,
            ctx,
        );
    }
}
