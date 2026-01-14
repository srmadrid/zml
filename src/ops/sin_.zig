const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// Performs in-place computation of the sine of `x`.
///
/// The `sin_` routine computes the sine of its input `x`, and stores the
/// result directly into `o`, automatically validating the provided context. The
/// operation is performed in the input's precision, and the result is then cast
/// to the output type. It supports both fixed-precision and arbitrary-precision
/// arithmetic, as well as structured data domains. The supported type
/// combinations are:
/// - **Numeric = sin(Numeric)**: scalar sine.
/// - **Matrix = sin(Matrix)**: matrix sine (not implemented yet).
/// - **Array = sin(Array)**: element-wise sine.
/// - **Expression = sin(Expression)**: symbolic sine.
///
/// Signature
/// ---------
/// ```zig
/// fn sin_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The operand to compute the sine of.
///
/// `ctx` (`anytype`):
/// A context struct providing necessary resources and configuration for the
/// operation. The required fields depend on the output and operand types. If
/// the context is missing required fields or contains unnecessary or wrongly
/// typed fields, the compiler will emit a detailed error message describing the
/// expected structure.
///
/// Returns
/// -------
/// `void`
///
/// Errors
/// ------
/// ``:
///
/// Notes
/// -----
/// When the output and input types are the same, aliasing is allowed.
///
/// When the type of the operand is of arbitrary precision, the context may
/// provide an optional pre-allocated buffer to store intermediate results,
/// avoiding repeated allocations in scenarios where `sin_` is called multiple
/// times. If no buffer is provided, the operation will allocate a temporary
/// buffer internally, using the allocator specified in the context. Aliasing
/// between `o` and the buffer is not checked, and will lead to extra
/// allocations.
pub inline fn sin_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sin_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.sin_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domain(O)) {
        .array => switch (comptime types.domain(X)) {
            .array => { // array = sin(array)
                comptime switch (types.numericType(types.Numeric(O))) {
                    .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int, .float, .cfloat => {
                            types.validateContext(@TypeOf(ctx), .{});
                        },
                        else => @compileError("zml.sin_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                    },
                    else => @compileError("zml.sin_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                };

                return array.sin_(
                    o,
                    x,
                    ctx,
                );
            },
            else => @compileError("zml.sin_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .numeric => switch (comptime types.domain(X)) {
            .numeric => { // numeric = sin(numeric)
                switch (comptime types.numericType(O)) {
                    .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.sin_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int, .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                float.sin(x),
                                .{},
                            ) catch unreachable;
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                cfloat.sin(x),
                                .{},
                            ) catch unreachable;
                        },
                        .integer => @compileError("zml.sin_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .rational => @compileError("zml.sin_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .real => @compileError("zml.sin_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.sin_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    else => @compileError("zml.sin_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.sin_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.sin_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
