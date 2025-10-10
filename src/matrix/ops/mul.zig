const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const matops = @import("../ops.zig");

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

    comptime if (!types.isMatrix(X) and !types.isMatrix(Y))
        @compileError("At least one of the arguments must be a matrix type");

    if (comptime (types.isMatrix(X) and types.isMatrix(Y)) or
        types.isVector(X) or types.isVector(Y))
    { // matrix * matrix  or  vector * matrix  or  matrix * vector
        comptime if (types.isArbitraryPrecision(C)) {
            @compileError("Arbitrary precision types not implemented yet");
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        };

        return linalg.matmul(allocator, x, y, ctx);
    } else {
        comptime if (types.isArbitraryPrecision(C)) { // scalar * matrix  or  matrix * scalar
            @compileError("Arbitrary precision types not implemented yet");
        } else {
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            }
        };

        return matops.apply2(
            allocator,
            x,
            y,
            ops.mul,
            ctx,
        );
    }
}
