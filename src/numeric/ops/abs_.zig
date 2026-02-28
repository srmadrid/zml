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

/// Performs in-place computation of the absolute value of a numeric `x` into
/// a numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.abs_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the absolute value of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X` and
///   `Y`. If the context is missing required fields or contains unnecessary or
///   wrongly typed fields, the compiler will emit a detailed error message
///   describing the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `O` and `X`.
///
/// #### `O` is not allocated
/// The context must be empty.
///
/// #### `O` is allocated, `O == X`, and `X.has_simple_abs` exists and is true
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, the operation will return an error unless
///   `o` and `x` refer to the same instance, in which case the operation will
///   be performed in-place without allocation.
///
/// #### `O` is allocated, and `O != X` or `X.has_simple_abs` does not exist or is false
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `O` is allocated and an allocator is provided.
/// * `integer.Error.AllocatorRequired`: If `O == X == Integer`, `o` and `x`
///   refer to different instances, and no allocator is provided in the context.
/// * `rational.Error.AllocatorRequired`: If `O == X == Rational`, `o` and `x`
///   refer to different instances, and no allocator is provided in the context.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `abs_` method. The expected
/// signature and behavior of `abs_` are as follows:
/// * `O` is not allocated: `fn abs_(*O, X) void`: Computes the absolute value
///   of `x` and stores it in `o`.
/// * `O` is allocated, `O == X` and `X.has_simple_abs` exists and is true: `fn abs_(?std.mem.Allocator, *O, X) !void`:
///   Computes the absolute value of `x` and stores it in `o`. If an allocator
///   is not provided, the operation must be performed in-place without
///   allocation or return an error.
/// * `O` is allocated, `O != X` or `X.has_simple_abs` does not exist or is false: `fn abs_(std.mem.Allocator, *O, X) !void`:
///   Computes the absolute value of `x` and stores it in `o`.
///
/// If neither `O` nor `X` implement the required `abs_` method, the function
/// will fall back to using `numeric.set` with the result of `numeric.abs`. In
/// this case, `O`, `X` and `ctx`  must adhere to the requirements of these
/// functions.
pub inline fn abs_(o: anytype, x: anytype, ctx: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O) or
        !types.isNumeric(types.Child(O)) or
        !types.isNumeric(X))
        @compileError("zml.numeric.abs_: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = types.Child(O);

    if (comptime types.isCustomType(O)) {
        if (comptime types.isCustomType(X)) { // O and X both custom
            if (comptime types.isAllocated(O)) {
                if (comptime O == X and @hasDecl(X, "has_simple_abs") and X.has_simple_abs) {
                    if (comptime !types.hasMethod(O, "abs_", fn (?std.mem.Allocator, *O, X) anyerror!void, &.{ std.mem.Allocator, *O, O })) {
                        try numeric.set(
                            o,
                            numeric.abs(x, .{}) catch unreachable,
                            ctx,
                        );

                        return;
                    }

                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{
                                .type = std.mem.Allocator,
                                .required = false,
                                .description = "The allocator to use for the custom numeric's memory allocation. If not provided, the operation will return an error unless o and x refer to the same instance, in which case the operation will be performed in-place without allocation.",
                            },
                        },
                    );

                    if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                        try O.abs_(ctx.allocator, o, x)
                    else
                        O.abs_(null, o, x) catch unreachable;

                    return;
                } else {
                    const Impl: type = comptime types.anyHasMethod(
                        &.{ O, X },
                        "abs_",
                        fn (std.mem.Allocator, *O, X) anyerror!void,
                        &.{ std.mem.Allocator, *O, X },
                    ) orelse {
                        var abs = try numeric.abs(x, if (types.isAllocated(numeric.Abs(X))) ctx else .{});
                        defer numeric.deinit(&abs, if (types.isAllocated(numeric.Abs(X))) ctx else .{});

                        try numeric.set(
                            o,
                            abs,
                            ctx,
                        );

                        return;
                    };

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

                    try Impl.abs_(ctx.allocator, o, x);

                    return;
                }
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ O, X },
                    "abs_",
                    fn (*O, X) void,
                    &.{ *O, X },
                ) orelse {
                    var abs = try numeric.abs(x, ctx);
                    defer numeric.deinit(&abs, ctx);

                    numeric.set(
                        o,
                        abs,
                        .{},
                    ) catch unreachable;

                    return;
                };

                comptime types.validateContext(@TypeOf(ctx), .{});

                Impl.abs_(o, x);

                return;
            }
        } else { // only O custom
            if (comptime types.isAllocated(O)) {
                if (comptime !types.hasMethod(O, "abs_", fn (std.mem.Allocator, *O, X) anyerror!void, &.{ std.mem.Allocator, *O, X })) {
                    var abs = try numeric.abs(x, if (types.isAllocated(numeric.Abs(X))) ctx else .{});
                    defer numeric.deinit(&abs, if (types.isAllocated(numeric.Abs(X))) ctx else .{});

                    try numeric.set(
                        o,
                        abs,
                        ctx,
                    );

                    return;
                }

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

                try O.abs_(ctx.allocator, o, x);

                return;
            } else {
                if (comptime !types.hasMethod(O, "abs_", fn (*O, X) void, &.{ *O, X })) {
                    var abs = try numeric.abs(x, ctx);
                    defer numeric.deinit(&abs, ctx);

                    numeric.set(
                        o,
                        abs,
                        .{},
                    ) catch unreachable;

                    return;
                }

                comptime types.validateContext(@TypeOf(ctx), .{});

                O.abs_(o, x);

                return;
            }
        }
    } else if (comptime types.isCustomType(X)) { // only X custom
        if (comptime types.isAllocated(O)) {
            if (comptime !types.hasMethod(X, "abs_", fn (std.mem.Allocator, O) anyerror!void, &.{ std.mem.Allocator, O })) {
                var abs = try numeric.abs(x, if (types.isAllocated(numeric.Abs(X))) ctx else .{});
                defer numeric.deinit(&abs, if (types.isAllocated(numeric.Abs(X))) ctx else .{});

                try numeric.set(
                    o,
                    abs,
                    ctx,
                );

                return;
            }

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

            try X.abs_(ctx.allocator, o, x);

            return;
        } else {
            if (comptime !types.hasMethod(X, "abs_", fn (X) O, &.{X})) {
                var abs = try numeric.abs(x, ctx);
                defer numeric.deinit(&abs, ctx);

                numeric.set(
                    o,
                    abs,
                    .{},
                ) catch unreachable;

                return;
            }

            comptime types.validateContext(@TypeOf(ctx), .{});

            X.abs_(o, x);

            return;
        }
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .dyadic, .cfloat => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    x,
                    ctx,
                ) catch unreachable;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    int.abs(x),
                    ctx,
                ) catch unreachable;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    float.abs(x),
                    ctx,
                ) catch unreachable;
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    dyadic.abs(x),
                    ctx,
                ) catch unreachable;
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    cfloat.abs(x),
                    ctx,
                ) catch unreachable;
            },
            .integer => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    integer.abs(null, x) catch unreachable,
                    ctx,
                ) catch unreachable;
            },
            .rational => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                numeric.set(
                    o,
                    rational.abs(null, x) catch unreachable,
                    ctx,
                ) catch unreachable;
            },
            .real => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
            .complex => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
            .custom => unreachable,
        },
        .integer => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    x,
                    ctx,
                );
            },
            .int => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    int.abs(x),
                    ctx,
                );
            },
            .float => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    float.abs(x),
                    ctx,
                );
            },
            .dyadic => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    dyadic.abs(x),
                    ctx,
                );
            },
            .cfloat => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    cfloat.abs(x),
                    ctx,
                );
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator for the integer's memory allocation. If not provided, the operation will return an error unless o and x refer to the same integer.",
                        },
                    },
                );

                if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    try integer.abs_(ctx.allocator, o, x)
                else
                    integer.abs_(null, o, x) catch unreachable;
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    rational.abs(null, x) catch unreachable,
                    ctx,
                );
            },
            .real => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
            .complex => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(X)) {
            .bool => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    x,
                    ctx,
                );
            },
            .int => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    int.abs(x),
                    ctx,
                );
            },
            .float => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    float.abs(x),
                    ctx,
                );
            },
            .dyadic => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    dyadic.abs(x),
                    ctx,
                );
            },
            .cfloat => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    cfloat.abs(x),
                    ctx,
                );
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                try numeric.set(
                    o,
                    integer.abs(null, x) catch unreachable,
                    ctx,
                );
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator for the rational's memory allocation. If not provided, the operation will return an error unless o and x refer to the same rational.",
                        },
                    },
                );

                if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    try rational.abs_(ctx.allocator, o, x)
                else
                    try rational.abs_(null, o, x);
            },
            .real => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
            .complex => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
            .custom => unreachable,
        },
        .real => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.abs_: not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet."),
        .custom => unreachable,
    }
}
