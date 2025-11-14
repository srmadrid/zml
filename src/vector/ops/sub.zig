const std = @import("std");

const types = @import("../../types.zig");
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const vecops = @import("../ops.zig");

/// Performs subtraction between two vectors, automatically handling any
/// combination of dense and sparse vectors.
///
/// Signature
/// ---------
/// ```zig
/// fn sub(x: X, y: Y, ctx: anytype) !Coerce(X, Y)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The left vector operand.
///
/// `y` (`anytype`):
/// The right vector operand.
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
/// `Coerce(@TypeOf(x), @TypeOf(y))`:
/// The result of the subtraction.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
///
/// `vector.Error.DimensionMismatch`:
/// If the two vectors do not have the same length.
pub inline fn sub(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isVector(@TypeOf(x)) or !types.isVector(@TypeOf(y)))
        @compileError("vector.sub: both x and y must be vectors, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime switch (types.numericType(types.Numeric(C))) {
        .bool => @compileError("vector.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
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
        ops.sub,
        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
    );
}
