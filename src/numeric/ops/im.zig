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

const constants = @import("../../constants.zig");

pub fn Im(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.im: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return types.Scalar(X),
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return types.Scalar(X),
        .custom => {
            if (comptime !types.hasMethod(X, "Im", fn (type) type, &.{}))
                @compileError("zml.numeric.im: " ++ @typeName(X) ++ " must implement `fn Im(type) type`");

            return X.Im(X);
        },
    }
}

/// Returns the imaginary part of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.im(x: X, ctx: anytype) !numeric.Im(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the imaginary part of.
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
/// `numeric.Im(@TypeOf(x))`: The imaginary part of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Im` and `im` methods. The expected
/// signature and behavior of `Im` and `im` are as follows:
/// * `fn Im(type) type`: Returns the return type of `im` for the custom
///   numeric type.
/// * Non-allocated: `fn im(X) X.Im(X)`: Returns the imaginary part of `x`.
/// * Allocated: `fn im(?std.mem.Allocator, X) !X.Im(X)`: Returns the imaginary
///   part of `x` as a newly allocated value, if the allocator is provided, or a
///   read-only view if not.
pub inline fn im(x: anytype, ctx: anytype) !numeric.Im(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Im(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return false;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0.0;
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return .zero;
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x.im;
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

            return constants.zero(R, ctx);
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

            return constants.zero(R, ctx);
        },
        .real => @compileError("zml.numeric.im: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.im: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "im", fn (?std.mem.Allocator, X) anyerror!R, &.{}))
                    @compileError("zml.numeric.im: " ++ @typeName(X) ++ " must implement `fn im(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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
                    X.im(ctx.allocator, x)
                else
                    X.im(null, x);
            } else {
                comptime if (!types.hasMethod(X, "im", fn (X) R, &.{}))
                    @compileError("zml.numeric.im: " ++ @typeName(X) ++ " must implement `fn im(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.im(x);
            }
        },
    }
}
