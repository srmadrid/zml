const std = @import("std");

const types = @import("../../types.zig");
const MulCoerce = types.MulCoerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const vecops = @import("../ops.zig");

const linalg = @import("../../linalg.zig");

/// Performs multiplication between two vectors, or between a vector and a
/// scalar, automatically handling any combination of dense and sparse vectors.
///
/// Signature
/// ---------
/// ```zig
/// fn mul(x: X, y: Y, ctx: anytype) !MulCoerce(X, Y)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// `ctx` (`anytype`):
/// A context struct providing necessary resources and configuration for the
/// operation. The required fields depend on the operand types. If the context
/// is missing required fields or contains unnecessary or wrongly typed fields,
/// the compiler will emit a detailed error message describing the expected
/// structure.
///
/// Returns
/// -------
/// `MulCoerce(@TypeOf(x), @TypeOf(y))`:
/// The result of the multiplication.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
///
/// `vector.Error.DimensionMismatch`:
/// If the two vectors do not have the same length. Can only happen when both
/// operands are vectors.
pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !MulCoerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isVector(@TypeOf(x)) and !types.isVector(@TypeOf(y)))
        @compileError("vector.mul: at least one of x or y must be a vector, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (comptime types.isVector(X) and types.isVector(Y)) { // vector * vector
        @compileError("vector * vector not implemented yet");

        // return linalg.dot(x, y, ctx);
    } else {
        comptime switch (types.numericType(types.Numeric(C))) {
            .bool => @compileError("vector.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .int, .float, .cfloat => {
                types.validateContext(@TypeOf(ctx), .{});
            },
            .integer, .rational, .real, .complex => {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            },
        };

        return vecops.apply2(
            allocator,
            x,
            y,
            ops.mul,
            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
        );
    }
}
