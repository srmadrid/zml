const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// Performs in-place computation of the hyperbolic arccosine of `x`.
///
/// The `acosh_` routine computes the hyperbolic arccosine of its input `x`, and
/// stores the result directly into `o`, automatically validating the provided
/// context. The operation is performed in the input's precision, and the result
/// is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = acosh(Numeric)**: scalar hyperbolic arccosine.
/// - **Matrix = acosh(Matrix)**: matrix hyperbolic arccosine (not implemented
///   yet).
/// - **Array = acosh(Array)**: element-wise hyperbolic arccosine.
/// - **Expression = acosh(Expression)**: symbolic hyperbolic arccosine.
///
/// Signature
/// ---------
/// ```zig
/// fn acosh_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The operand to compute the hyperbolic arccosine of.
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
/// avoiding repeated allocations in scenarios where `acosh_` is called multiple
/// times. If no buffer is provided, the operation will allocate a temporary
/// buffer internally, using the allocator specified in the context. Aliasing
/// between `o` and the buffer is not checked, and will lead to extra
/// allocations.
pub inline fn acosh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.acosh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.acosh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domain(O)) {
        .array => switch (comptime types.domain(X)) {
            .array => { // array = acos(array)
                comptime switch (types.numericType(types.Numeric(O))) {
                    .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int, .float, .cfloat => {
                            types.validateContext(@TypeOf(ctx), .{});
                        },
                        else => @compileError("zml.acosh_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                    },
                    else => @compileError("zml.acosh_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                };

                return array.acosh_(
                    o,
                    x,
                    ctx,
                );
            },
            else => @compileError("zml.acosh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .numeric => switch (comptime types.domain(X)) {
            .numeric => { // numeric = acos(numeric)
                switch (comptime types.numericType(O)) {
                    .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.acosh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int, .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                float.acosh(x),
                                .{},
                            ) catch unreachable;
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                cfloat.acosh(x),
                                .{},
                            ) catch unreachable;
                        },
                        .integer => @compileError("zml.acosh_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .rational => @compileError("zml.acosh_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .real => @compileError("zml.acosh_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.acosh_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    else => @compileError("zml.acosh_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.acosh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.acosh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
