const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");

/// The return type of the `log` routine for an input of type `X`.
pub fn Log(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Log(types.Numeric(X))),
        .matrix => @compileError("zml.Log not implemented for matrices yet"),
        .vector => @compileError("zml.Log not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the natural logarithm of `x`.
///
/// The `log` routine computes the natural logarithm of its input `x`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar natural logarithm.
/// - **Matrix**: matrix natural logarithm (not implemented yet).
/// - **Array**: element-wise natural logarithm.
///
/// Signature
/// ---------
/// ```zig
/// fn log(x: X, ctx: anytype) !Log(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the natural logarithm of.
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
/// `Log(@TypeOf(x))`:
/// The natural logarithm of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn log(
    x: anytype,
    ctx: anytype,
) !Log(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.log not defined for " ++ @typeName(X));

    switch (comptime types.domainType(X)) {
        .array => {
            comptime switch (types.numericType(types.Numeric(X))) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
                else => @compileError("zml.log for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.log(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.log(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.log(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.log(x);
            },
            else => @compileError("zml.log for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
