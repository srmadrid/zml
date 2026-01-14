const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// Performs in-place computation of the base-2 logarithm of `x`.
///
/// The `log2_` routine computes the base-2 logarithm of its input `x`, and
/// stores the result directly into `o`, automatically validating the provided
/// context. The operation is performed in the input's precision, and the result
/// is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = log2(Numeric)**: scalar base-2 logarithm.
/// - **Matrix = log2(Matrix)**: matrix base-2 logarithm (not implemented yet).
/// - **Array = log2(Array)**: element-wise base-2 logarithm.
/// - **Expression = log2(Expression)**: symbolic base-2 logarithm.
///
/// Signature
/// ---------
/// ```zig
/// fn log2_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The operand to compute the base-2 logarithm of.
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
/// avoiding repeated allocations in scenarios where `log2_` is called multiple
/// times. If no buffer is provided, the operation will allocate a temporary
/// buffer internally, using the allocator specified in the context. Aliasing
/// between `o` and the buffer is not checked, and will lead to extra
/// allocations.
pub inline fn log2_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.log2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domain(O)) {
        .array => switch (comptime types.domain(X)) {
            .array => { // array = log2(array)
                comptime switch (types.numericType(types.Numeric(O))) {
                    .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int, .float, .cfloat => {
                            types.validateContext(@TypeOf(ctx), .{});
                        },
                        else => @compileError("zml.log2_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                    },
                    else => @compileError("zml.log2_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                };

                return array.log2_(
                    o,
                    x,
                    ctx,
                );
            },
            else => @compileError("zml.log2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .numeric => switch (comptime types.domain(X)) {
            .numeric => { // numeric = log2(numeric)
                switch (comptime types.numericType(O)) {
                    .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.log2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int, .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                float.log2(x),
                                .{},
                            ) catch unreachable;
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                cfloat.log2(x),
                                .{},
                            ) catch unreachable;
                        },
                        .integer => @compileError("zml.log2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .rational => @compileError("zml.log2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .real => @compileError("zml.log2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.log2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    else => @compileError("zml.log2_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.log2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.log2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
