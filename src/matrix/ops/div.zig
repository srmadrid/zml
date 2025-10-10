const std = @import("std");

const types = @import("../../types.zig");
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const matops = @import("../ops.zig");

pub inline fn div(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(Numeric(X), Numeric(Y));

    comptime if (!types.isMatrix(X) and types.isMatrix(Y))
        @compileError("First argument must be a matrix type and second argument must be a scalar type");

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("Arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return matops.apply2(
        allocator,
        x,
        y,
        ops.div,
        ctx,
    );
}
