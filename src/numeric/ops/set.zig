const std = @import("std");

const types = @import("../../types.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const cfloat = @import("../../cfloat.zig");
const integer = @import("../../integer.zig");
const rational = @import("../../rational.zig");
const real = @import("../../real.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

/// Sets the value of `o` to `x`.
///
/// ## Signature
/// ```zig
/// numeric.set(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The input operand.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X` and
///   `Y`. If the context is missing required fields or contains unnecessary or
///   wrongly typed fields, the compiler will emit a detailed error message
///   describing the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `O`.
///
/// #### `O` is not allocated
/// The context must be empty.
///
/// #### `O` is allocated
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `O` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` must implement the required `set` method. The expected signature
/// and behavior of `set` are as follows:
/// * `O` is not allocated: `fn set(*O, X) void`: Sets the value of `o` to `x`.
/// * `O` is allocated: `fn set(std.mem.Allocator, *O, X) !void`: Sets the value
///   of `o` to `x`, using the provided allocator for any necessary memory
///   management.
pub inline fn set(o: anytype, x: anytype, ctx: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O) or
        !types.isNumeric(types.Child(O)) or
        !types.isNumeric(X))
        @compileError("zml.numeric.set: o must be a mutable one-itme pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = types.Child(O);

    if (comptime types.isCustomType(O)) {
        if (comptime types.isCustomType(X)) { // O and X both custom
            if (comptime types.isAllocated(O)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ O, X },
                    "set",
                    fn (std.mem.Allocator, *O, X) anyerror!void,
                    &.{ std.mem.Allocator, *O, X },
                ) orelse
                    @compileError("zml.numeric.set: " ++ @typeName(O) ++ " or " ++ @typeName(X) ++ " must implement `fn set(std.mem.Allocator, *" ++ @typeName(O) ++ ", " ++ @typeName(X) ++ ") !void`");

                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the custom numeric's memory allocation.",
                        },
                    },
                );

                return Impl.set(ctx.allocator, o, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ O, X },
                    "set",
                    fn (*O, X) void,
                    &.{ *O, X },
                ) orelse
                    @compileError("zml.numeric.set: " ++ @typeName(O) ++ " or " ++ @typeName(X) ++ " must implement `fn set(*" ++ @typeName(O) ++ ", " ++ @typeName(X) ++ ") void`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.set(o, x);
            }
        } else { // only O custom
            if (comptime types.isAllocated(O)) {
                comptime if (!types.hasMethod(O, "set", fn (std.mem.Allocator, *O, X) anyerror!void, &.{ std.mem.Allocator, *O, X }))
                    @compileError("zml.numeric.set: " ++ @typeName(O) ++ " must implement `fn set(std.mem.Allocator, *" ++ @typeName(O) ++ ", " ++ @typeName(X) ++ ") !void`");

                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the custom numeric's memory allocation.",
                        },
                    },
                );

                return O.set(ctx.allocator, o, x);
            } else {
                comptime if (!types.hasMethod(O, "set", fn (*O, X) void, &.{ *O, X }))
                    @compileError("zml.numeric.set: " ++ @typeName(O) ++ " must implement `fn set(*" ++ @typeName(O) ++ ", " ++ @typeName(X) ++ ") void`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return O.set(o, x);
            }
        }
    } else if (comptime types.isCustomType(X)) { // only X custom
        if (comptime types.isAllocated(O)) {
            comptime if (!types.hasMethod(X, "set", fn (std.mem.Allocator, *O, X) anyerror!void, &.{ std.mem.Allocator, *O, X }))
                @compileError("zml.numeric.set: " ++ @typeName(X) ++ " must implement `fn set(std.mem.Allocator, *" ++ @typeName(O) ++ ", " ++ @typeName(X) ++ ") !void`");

            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the custom numeric's memory allocation.",
                    },
                },
            );

            return X.set(ctx.allocator, o, x);
        } else {
            comptime if (!types.hasMethod(X, "set", fn (*O, X) void, &.{ *O, X }))
                @compileError("zml.numeric.set: " ++ @typeName(X) ++ " must implement `fn set(*" ++ @typeName(O) ++ ", " ++ @typeName(X) ++ ") void`");

            comptime types.validateContext(@TypeOf(ctx), .{});

            return X.set(o, x);
        }
    }

    switch (comptime types.numericType(O)) {
        .bool => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(bool, x);
            },
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => unreachable,
        },
        .cfloat => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .real => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .complex => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = types.scast(O, x);
            },
            .custom => unreachable,
        },
        .integer => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the integer's memory deallocation. Must be the same allocator used to initialize it.",
                    },
                },
            );

            try integer.Integer.set(o, ctx.allocator, x);
        },
        .rational => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the rational's memory deallocation. Must be the same allocator used to initialize it.",
                    },
                },
            );

            try rational.Rational.set(o, ctx.allocator, x, 1);
        },
        .real => @compileError("zml.numeric.set: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.set: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
        .custom => unreachable,
    }
}
