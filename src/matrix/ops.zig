const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const EnsureMatrix = types.EnsureMatrix;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const matrix = @import("../array.zig");

const general = @import("general.zig");
const symmetric = @import("symmetric.zig");
const hermitian = @import("hermitian.zig");
const triangular = @import("triangular.zig");
const diagonal = @import("diagonal.zig");
const banded = @import("banded.zig");
const tridiagonal = @import("tridiagonal.zig");
const sparse = @import("sparse.zig");

const linalg = @import("../linalg.zig");

// Example
pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    _ = allocator;
    _ = opts;
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(Numeric(X), Numeric(Y));

    comptime if (!types.isMatrix(X) and !types.isMatrix(Y))
        @compileError("At least one of the arguments must be a matrix type");

    if (comptime types.isMatrix(X) and types.isMatrix(Y)) { // matrix * matrix
        comptime if (types.isArbitraryPrecision(C)) {
            @compileError("Arbitrary precision types not implemented yet");
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        };

        // return linalg.matmul(...);
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

        // return apply2(...);
    }
}
