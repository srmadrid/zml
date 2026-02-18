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

pub fn Conj(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.conj: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return X,
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return X,
        .custom => {
            if (comptime !types.hasMethod(X, "Conj", fn (type) type, &.{}))
                @compileError("zml.numeric.conj: " ++ @typeName(X) ++ " must implement `fn Conj(type) type`");

            return X.Conj(X);
        },
    }
}

/// Returns the complex conjugate of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.conj(x: X, ctx: anytype) !numeric.Conj(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the complex conjugate of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `X`.
///
/// #### `X` is not allocated
/// The context must be empty.
///
/// #### `X` is allocated
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// ## Returns
/// `numeric.Conj(@TypeOf(x))`: The complex conjugate of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Conj` and `conj` methods. The expected
/// signature and behavior of `Conj` and `conj` are as follows:
/// * `fn Conj(type) type`: Returns the return type of `conj` for the custom
///   numeric type.
/// * Non-allocated: `fn conj(X) X.Conj(X)`: Returns the complex conjugate of
///   `x`.
/// * Allocated: `fn conj(?std.mem.Allocator, X) !X.Conj(X)`: Returns the
///   complex conjugate of `x` as a newly allocated value, if the allocator is
///   provided, or a read-only view if not.
pub inline fn conj(x: anytype, ctx: anytype) !numeric.Conj(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Conj(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x.conj();
        },
        .integer => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = false,
                        .description = "The allocator to use for the integer's memory allocation. If not provided, a view will be returned.",
                    },
                },
            );

            return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                x.copy(ctx.allocator)
            else
                x.view();
        },
        .rational => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = false,
                        .description = "The allocator to use for the rational's memory allocation. If not provided, a view will be returned.",
                    },
                },
            );

            return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                x.copy(ctx.allocator)
            else
                x.view();
        },
        .real => @compileError("zml.numeric.conj: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.conj: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "conj", fn (?std.mem.Allocator, X) anyerror!R, &.{}))
                    @compileError("zml.numeric.conj: " ++ @typeName(X) ++ " must implement `fn conj(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the custom numeric's memory allocation. If not provided, a view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    X.conj(ctx.allocator, x)
                else
                    X.conj(null, x);
            } else {
                comptime if (!types.hasMethod(X, "conj", fn (X) R, &.{}))
                    @compileError("zml.numeric.conj: " ++ @typeName(X) ++ " must implement `fn conj(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.conj(x);
            }
        },
    }
}
