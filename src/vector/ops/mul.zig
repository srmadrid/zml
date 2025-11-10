const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const vecops = @import("../ops.zig");

const linalg = @import("../../linalg.zig");

pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !MulCoerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = MulCoerce(Numeric(X), Numeric(Y));

    comptime if (!types.isVector(X) and !types.isVector(Y))
        @compileError("At least one of the arguments must be a vector type");

    if (comptime types.isVector(X) and types.isVector(Y)) { // vector * vector
        comptime if (types.isArbitraryPrecision(C)) {
            @compileError("Arbitrary precision types not implemented yet");
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        };

        return linalg.dot(x, y, ctx);
    } else {
        comptime if (types.isArbitraryPrecision(C)) { // scalar * vector  or  vector * scalar
            @compileError("Arbitrary precision types not implemented yet");
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
            ops.mul,
            ctx,
        );
    }
}
