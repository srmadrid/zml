const std = @import("std");

const types = @import("../../types.zig");
const int = @import("../../int.zig");
const ops = @import("../../ops.zig");

const arrops = @import("../ops.zig");

pub inline fn mul_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.array.mul_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime !types.isArray(O))
        @compileError("zml.array.mul_ requires the output to be an array, got " ++ @typeName(O));

    const C: type = types.Coerce(types.Numeric(X), types.Numeric(Y));

    if (comptime types.isArbitraryPrecision(types.Numeric(O))) {
        comptime if (types.isArbitraryPrecision(C)) {
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        };

        return arrops.apply2_(
            o,
            x,
            y,
            ops.mul_,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    } else {
        comptime if (types.isArbitraryPrecision(C)) {
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(types.Numeric(C)) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            }
        };

        return arrops.apply2_(
            o,
            x,
            y,
            ops.mul_,
            ctx,
        );
    }
}
