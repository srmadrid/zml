const std = @import("std");

const types = @import("../../types.zig");
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const vecops = @import("../ops.zig");

/// Performs division of a vector by a scalar, automatically handling any
/// combination of dense and sparse vectors.
///
/// Signature
/// ---------
/// ```zig
/// fn div(x: X, y: Y, ctx: anytype) !Coerce(X, Y)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The left vector operand.
///
/// `y` (`anytype`):
/// The right scalar operand.
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
/// The result of the division.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails.
pub inline fn div(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(Numeric(X), Numeric(Y));

    comptime if (!(types.isVector(@TypeOf(x)) and types.isNumeric(@TypeOf(y))))
        @compileError("vector.div: x must be a vector and y must be a numeric, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime switch (types.numericType(types.Numeric(C))) {
        .bool => @compileError("vector.div not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
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
        ops.div,
        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
    );
}
