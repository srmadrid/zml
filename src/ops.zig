const std = @import("std");

const types = @import("types.zig");
const cast = types.cast;
const scast = types.scast;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Coerce = types.Coerce;
const CoerceToArray = types.CoerceToArray;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const needsAllocator = types.needsAllocator;
const validateContext = types.validateContext;
const getFieldOrDefault = types.getFieldOrDefault;
const mixStructs = types.mixStructs;
const stripStruct = types.stripStruct;

const int = @import("int.zig");
const float = @import("float.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");
const complex = @import("complex.zig");
const array = @import("array.zig");

/// Adds two values of any two supported types.
///
/// The function supports addition for values of any combination of supported
/// numeric types, `Array`s and slices. The operation performed is
///
/// ```zig
///     x + y
/// ```
///
/// Parameters
/// ----------
/// x (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The left-hand side operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// ctx (`anytype`): Context for the addition operation.
///
/// Returns
/// -------
/// `Coerce(@TypeOf(x), @TypeOf(y))`: The result of the addition operation.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision or an `Array`.
///
/// Raises
/// ------
/// `@compileError`: If the addition operation is not defined for the types of
/// the inputs, if the types of the inputs cannot be coerced to a common type,
/// or if the types of the inputs are not supported numeric types, or `Array`s
/// or slices of unsupported types.
///
/// See Also
/// --------
/// `add_`: For in-place addition of two values.
///
/// `int.add`: For addition of an `int` and another `int` or `bool`.
///
/// `float.add`: For addition of a `float` and another `float`, `int` or `bool`.
///
/// `cfloat.add`: For addition of a `cfloat` and another `cfloat`, `float`,
/// `int` or `bool`.
///
/// `integer.add`: For addition of an `integer` and another `integer`, `bool`
/// or `int`.
///
/// `rational.add`: For addition of a `rational` and another `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `real.add`: For addition of a `real` and another `real`, `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `complex.add`: For addition of a `complex` and another `complex`, `real`,
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.add`: For addition of an `expression` and another  `expression`,
/// `complex`, `real`, `rational`, `integer`, `float`, `int` or `bool`.
///
/// `array.add`: For addition of two `Array`s or slices.
///
/// Notes
/// -----
/// The addition is performed with the precision of the coerced type of the
/// inputs.
pub inline fn add(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            if (types.numericType(Numeric(C)) == .int) {
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .mode = .{ .type = int.Mode, .required = false },
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .order = .{ .type = ?array.Order, .required = false },
                    },
                );
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .order = .{ .type = ?array.Order, .required = false },
                    },
                );
            }
        };

        return array.add(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => {
            comptime validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );

            return int.add(x, y, getFieldOrDefault(ctx, "mode", int.Mode, .default));
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.add(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.add(x, y);
        },
        else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

/// Adds two values of any two supported types in-place.
///
/// The function supports addition for values of any combination of supported
/// numeric types, `Array`s and slices, and stores the result in the output
/// pointer. The operation performed is
///
/// ```zig
///     o = x + y
/// ```
///
/// Parameters
/// ----------
/// o (pointer to `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The output pointer where the
/// result of the addition operation will be stored. If the output type is of
/// arbitrary precision or an `Array`, it must be initialized.
///
/// x (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The left-hand side operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the addition operation.
/// - `mode` (`int.Mode`): The mode of the addition operation. Only needed when
/// adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision. May not be used if the output has enouph memory allocated
/// already.
///
/// Returns
/// -------
/// `void`: The result of the addition operation is stored in the output
/// pointer.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision.
///
/// Raises
/// ------
/// `@compileError`: If the addition operation is not defined for the types of
/// the inputs, if the types of the inputs cannot be coerced to a common type,
/// if the types of the inputs are not supported numeric types, `Array`s or
/// slices, if the output pointer is not a mutable pointer, if the output's
/// child type is not a supported numeric type, `Array` or slice, or if the
/// coerced type is an `Array` and the output's child type is not an `Array`.
///
/// See Also
/// --------
/// `add`: For addition of two values and returning the result.
///
/// Notes
/// -----
/// The addition is performed with the precision of the coerced type of the
/// inputs, and the result is cast to the output type if necessary. This cast
/// is not checked for safety.
///
/// Aliasing is allowed.
pub inline fn add_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.add_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting
                if (types.numericType(Numeric(C)) == .int) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .mode = .{ .type = int.Mode, .required = false },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                }
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                if (types.numericType(Numeric(C)) == .int) {
                    validateContext(
                        @TypeOf(ctx),
                        .{ .mode = .{ .type = int.Mode, .required = false } },
                    );
                } else {
                    validateContext(@TypeOf(ctx), .{});
                }
            }
        };

        return array.add_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.add_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.add_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => {
                comptime validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );

                o.* = scast(O, int.add(x, y, getFieldOrDefault(ctx, "mode", int.Mode, .default)));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.add(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.add(x, y));
            },
            else => @compileError("zml.add_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.add_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

/// Subtracts two values of any two supported types.
///
/// The function supports subtraction for values of any combination of supported
/// numeric types, `Array`s and slices. The operation performed is
///
/// ```zig
///     x - y
/// ```
///
/// Parameters
/// ----------
/// x (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The left-hand side operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the subtraction operation.
/// - `mode` (`int.Mode`): The mode of the subtraction operation. Only needed
/// when adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision.
/// - `writeable` (`bool`): Whether the output should be writeable. Only needed
/// if the output type is an `Array`.
///
/// Returns
/// -------
/// `Coerce(@TypeOf(x), @TypeOf(y))`: The result of the subtraction operation.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision.
///
/// Raises
/// ------
/// `@compileError`: If the subtraction operation is not defined for the types
/// of the inputs, if the types of the inputs cannot be coerced to a common
/// type, or if the types of the inputs are not supported numeric types,
/// `Array`s or slices.
///
/// See Also
/// --------
/// `sub_`: For in-place subtraction of two values.
///
/// `int.sub`: For subtraction of an `int` and another `int` or `bool`.
///
/// `float.sub`: For subtraction of a `float` and another `float`, `int` or
/// `bool`.
///
/// `cfloat.sub`: For subtraction of a `cfloat` and another `cfloat`, `float`,
/// `int` or `bool`.
///
/// `integer.sub`: For subtraction of an `integer` and another `integer`, `bool`
/// or `int`.
///
/// `rational.sub`: For subtraction of a `rational` and another `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `real.sub`: For subtraction of a `real` and another `real`, `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `complex.sub`: For subtraction of a `complex` and another `complex`, `real`,
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.sub`: For subtraction of an `expression` and another
/// `expression`, `complex`, `real`, `rational`, `integer`, `float`, `int` or
/// `bool`.
///
/// `array.sub`: For subtraction of two `Array`s or slices.
///
/// Notes
/// -----
/// The subtraction is performed with the precision of the coerced type of the
/// inputs.
pub inline fn sub(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            if (types.numericType(Numeric(C)) == .int) {
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .mode = .{ .type = int.Mode, .required = false },
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .order = .{ .type = ?array.Order, .required = false },
                    },
                );
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .order = .{ .type = ?array.Order, .required = false },
                    },
                );
            }
        };

        return array.sub(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => {
            comptime validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );

            return int.sub(x, y, getFieldOrDefault(ctx, "mode", int.Mode, .default));
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.sub(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.sub(x, y);
        },
        else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

/// Sibtracts two values of any two supported types in-place.
///
/// The function supports subtraction for values of any combination of supported
/// numeric types, `Array`s and slices, and stores the result in the output
/// pointer. The operation performed is
///
/// ```zig
///     o = x - y
/// ```
///
/// Parameters
/// ----------
/// o (pointer to `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The output pointer where the
/// result of the subtraction operation will be stored. If the output type is of
/// arbitrary precision, an `Array` or a slice, it must be initialized.
///
/// x (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The left-hand side operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the subtraction operation.
/// - `mode` (`int.Mode`): The mode of the subtraction operation. Only needed
/// when adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for reallocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision. May not be used if the output has enouph memory allocated
/// already.
///
/// Returns
/// -------
/// `void`: The result of the addition operation is stored in the output
/// pointer.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision.
/// `array.Error.Bla`: Put errors when `array.add_` is implemented.
///
/// Raises
/// ------
/// `@compileError`: If the addition operation is not defined for the types of
/// the inputs, if the types of the inputs cannot be coerced to a common type,
/// if the types of the inputs are not supported numeric types, `Array`s or
/// slices, if the output pointer is not a mutable pointer, if the output's
/// child type is not a supported numeric type, `Array` or slice, or if the
/// coerced type is an `Array` and the output's child type is not an `Array`.
///
/// See Also
/// --------
/// `add`: For addition of two values and returning the result.
///
/// Notes
/// -----
/// The addition is performed with the precision of the coerced type of the
/// inputs, and the result is cast to the output type if necessary. This cast
/// is not checked for safety.
pub inline fn sub_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sub_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting
                if (types.numericType(Numeric(C)) == .int) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .mode = .{ .type = int.Mode, .required = false },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                }
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                if (types.numericType(Numeric(C)) == .int) {
                    validateContext(
                        @TypeOf(ctx),
                        .{ .mode = .{ .type = int.Mode, .required = false } },
                    );
                } else {
                    validateContext(@TypeOf(ctx), .{});
                }
            }
        };

        return array.sub_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.sub_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.sub_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => {
                comptime validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );

                o.* = scast(O, int.sub(x, y, getFieldOrDefault(ctx, "mode", int.Mode, .default)));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.sub(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.sub(x, y));
            },
            else => @compileError("zml.sub_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.sub_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

/// Multiplies two values of any two supported types.
///
/// The function supports multiplication for values of any combination of
/// supported numeric types, `Array`s and slices. The operation performed is
///
/// \begin{equation*}
///     x \cdot y
/// \end{equation*}
///
/// Parameters
/// ----------
/// x (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The left-hand side operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the multiplication operation.
/// - `mode` (`int.Mode`): The mode of the multiplication operation. Only needed
/// when adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision.
///
/// Returns
/// -------
/// `Coerce(@TypeOf(x), @TypeOf(y))`: The result of the multiplication
/// operation.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision.
///
/// Raises
/// ------
/// `@compileError`: If the multiplication operation is not defined for the
/// types of the inputs, if the types of the inputs cannot be coerced to a
/// common type, or if the types of the inputs are not supported numeric types,
/// `Array`s or slices.
///
/// See Also
/// --------
/// `mul_`: For in-place multiplication of two values.
///
/// `int.mul`: For multiplication of an `int` and another `int` or `bool`.
///
/// `float.mul`: For multiplication of a `float` and another `float`, `int` or
/// `bool`.
///
/// `cfloat.mul`: For multiplication of a `cfloat` and another `cfloat`,
/// `float`, `int` or `bool`.
///
/// `integer.mul`: For multiplication of an `integer` and another `integer`,
/// `bool` or `int`.
///
/// `rational.mul`: For multiplication of a `rational` and another `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `real.mul`: For multiplication of a `real` and another `real`, `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `complex.mul`: For multiplication of a `complex` and another `complex`,
/// `real`, `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.mul`: For multiplication of an `expression` and another
/// `expression`, `complex`, `real`, `rational`, `integer`, `float`, `int` or
/// `bool`.
///
/// `array.mul`: For multiplication of two `Array`s or slices.
///
/// Notes
/// -----
/// The multiplication is performed with the precision of the coerced type of
/// the inputs.
pub inline fn mul(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            if (types.numericType(Numeric(C)) == .int) {
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .mode = .{ .type = int.Mode, .required = false },
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .order = .{ .type = ?array.Order, .required = false },
                    },
                );
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .order = .{ .type = ?array.Order, .required = false },
                    },
                );
            }
        };

        return array.mul(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => {
            comptime validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );

            return int.mul(x, y, getFieldOrDefault(ctx, "mode", int.Mode, .default));
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.mul(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.mul(x, y);
        },
        else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn mul_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.mul_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting
                if (types.numericType(Numeric(C)) == .int) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .mode = .{ .type = int.Mode, .required = false },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                }
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                if (types.numericType(Numeric(C)) == .int) {
                    validateContext(
                        @TypeOf(ctx),
                        .{ .mode = .{ .type = int.Mode, .required = false } },
                    );
                } else {
                    validateContext(@TypeOf(ctx), .{});
                }
            }
        };

        return array.mul_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.mul_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.mul_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => {
                comptime validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );

                o.* = scast(O, int.mul(x, y, getFieldOrDefault(ctx, "mode", int.Mode, .default)));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.mul(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.mul(x, y));
            },
            else => @compileError("zml.mul_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.mul_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

/// Divides two values of any two supported types.
///
/// The function supports division for values of any combination of supported
/// numeric types, `Array`s and slices. The operation performed is
///
/// \begin{equation*}
///     x / y
/// \end{equation*}
///
/// Parameters
/// ----------
/// x (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The left-hand side operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the division operation.
/// - `mode` (`int.Mode`): The mode of the division operation. Only needed when
/// adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision.
///
/// Returns
/// -------
/// `Coerce(@TypeOf(x), @TypeOf(y))`: The result of the division operation.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision.
///
/// Raises
/// ------
/// `@compileError`: If the division operation is not defined for the types of
/// the inputs, if the types of the inputs cannot be coerced to a common type,
/// or if the types of the inputs are not supported numeric types, `Array`s or
/// slices.
///
/// See Also
/// --------
/// `div_`: For in-place division of two values.
///
/// `int.div`: For division of an `int` and another `int` or `bool`.
///
/// `float.div`: For division of a `float` and another `float`, `int` or `bool`.
///
/// `cfloat.div`: For division of a `cfloat` and another `cfloat`, `float`,
/// `int` or `bool`.
///
/// `integer.div`: For division of an `integer` and another `integer`, `bool`
/// or `int`.
///
/// `rational.div`: For division of a `rational` and another `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `real.div`: For division of a `real` and another `real`, `rational`,
/// `integer`, `float`, `int` or `bool`.
///
/// `complex.div`: For division of a `complex` and another `complex`, `real`,
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.div`: For division of an `expression` and another  `expression`,
/// `complex`, `real`, `rational`, `integer`, `float`, `int` or `bool`.
///
/// `array.div`: For division of two `Array`s or slices.
///
/// Notes
/// -----
/// The division is performed with the precision of the coerced type of the
/// inputs.
pub inline fn div(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.mul(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.mul(x, y, getFieldOrDefault(ctx, "mode", int.Mode, .default));
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.mul(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.mul(x, y);
        },
        else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn div_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.div_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting

                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.div_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.div_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.div_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.div(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.div(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.div(x, y));
            },
            else => @compileError("zml.div_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.div_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn eq(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime validateContext(
            @TypeOf(ctx),
            .{
                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                .order = .{ .type = ?array.Order, .required = false },
            },
        );

        return array.eq(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
        );
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.eq does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return x == y;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.eq(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.eq(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.eq(x, y);
        },
        else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn eq_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.eq_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.eq_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.eq_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.eq_ does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x == y);
            },
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.eq(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.eq(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.eq(x, y));
            },
            else => @compileError("zml.eq_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.eq_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn ne(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime validateContext(
            @TypeOf(ctx),
            .{
                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                .order = .{ .type = ?array.Order, .required = false },
            },
        );

        return array.ne(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
        );
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.ne does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return x == y;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.ne(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.ne(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.ne(x, y);
        },
        else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn ne_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ne_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.ne_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.ne_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.ne_ does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x == y);
            },
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.ne(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.ne(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.ne(x, y));
            },
            else => @compileError("zml.ne_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.ne_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn lt(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime validateContext(
            @TypeOf(ctx),
            .{
                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                .order = .{ .type = ?array.Order, .required = false },
            },
        );

        return array.lt(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
        );
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.lt is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => unreachable,
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.lt(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.lt(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.lt(x, y);
        },
        else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn lt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.lt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.lt_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.lt_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.lt_ is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => unreachable,
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.lt(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.lt(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.lt(x, y));
            },
            else => @compileError("zml.lt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.lt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn le(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime validateContext(
            @TypeOf(ctx),
            .{
                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                .order = .{ .type = ?array.Order, .required = false },
            },
        );

        return array.le(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
        );
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.le is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => unreachable,
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.le(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.le(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.le(x, y);
        },
        else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn le_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.le_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.le_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.le_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.le_ is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => unreachable,
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.le(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.le(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.le(x, y));
            },
            else => @compileError("zml.le_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.le_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn gt(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime validateContext(
            @TypeOf(ctx),
            .{
                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                .order = .{ .type = ?array.Order, .required = false },
            },
        );

        return array.gt(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
        );
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.gt is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => unreachable,
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.gt(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.gt(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.gt(x, y);
        },
        else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn gt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.gt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.gt_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.gt_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.lt is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => unreachable,
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.gt(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.gt(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.gt(x, y));
            },
            else => @compileError("zml.gt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.gt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn ge(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime validateContext(
            @TypeOf(ctx),
            .{
                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                .order = .{ .type = ?array.Order, .required = false },
            },
        );

        return array.ge(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
        );
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.lt is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => unreachable,
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.ge(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.ge(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.ge(x, y);
        },
        else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn ge_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ge_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.ge_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.ge_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if (types.numericType(X) == .bool or types.numericType(Y) == .bool)
        @compileError("zml.lt is not defined for `bool` types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => unreachable,
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.ge(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.ge(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.ge(x, y));
            },
            else => @compileError("zml.ge_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.ge_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn max(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.max(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.max does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return x or y;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.max(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.max(x, y);
        },
        .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn max_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.max_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.max_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.max_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.max_ does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x or y);
            },
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.max(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.max(x, y));
            },
            .cfloat => @compileError("zml.max_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .complex => @compileError("zml.max_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            else => @compileError("zml.max_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
        },
        else => @compileError("zml.max_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn min(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.min(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.min does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(C)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return x or y;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.min(x, y);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.min(x, y);
        },
        .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn min_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.min_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        };

        return array.min_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.min_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    comptime if ((types.numericType(X) == .bool and !types.numericType(Y) == .bool) or
        (types.numericType(Y) == .bool and !types.numericType(X) == .bool))
        @compileError("zml.min_ does not support comparing `bool` with other numeric types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x or y);
            },
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.min(x, y));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.min(x, y));
            },
            .cfloat => @compileError("zml.min_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .complex => @compileError("zml.min_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            else => @compileError("zml.min_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
        },
        else => @compileError("zml.min_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

// Basic operations
pub inline fn abs(
    x: anytype,
    ctx: anytype,
) !CoerceToArray(@TypeOf(x), Scalar(Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                    .copy = .{ .type = bool, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.abs(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs not defined for " ++ @typeName(X)),
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.abs(x);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.abs(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.abs(x);
        },
        else => @compileError("zml.abs not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn abs_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .copy = .{ .type = bool, .required = false },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.abs_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.abs_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.abs_ not defined for " ++ @typeName(X) ++ " input type"),
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.abs(x));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.abs(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.abs(x));
            },
            else => @compileError("zml.abs_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.abs_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn abs2(
    x: anytype,
    ctx: anytype,
) !CoerceToArray(@TypeOf(x), Scalar(Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.abs2(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs2 not defined for " ++ @typeName(X)),
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return int.mul(x, x, .default);
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.pow(x, 2);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.abs2(x);
        },
        else => @compileError("zml.abs2 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn abs2_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.abs2_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.abs2_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(X) ++ " input type"),
            .int => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, int.mul(x, x, .default));
            },
            .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.pow(x, 2));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.abs2(x));
            },
            else => @compileError("zml.abs2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.abs2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

// Exponential functions
pub inline fn exp(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.exp(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.exp(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.exp(x);
        },
        else => @compileError("zml.exp not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.exp_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.exp_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.exp(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.exp(x));
            },
            else => @compileError("zml.exp_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.exp_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn exp10(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.exp10(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp10 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.exp10(x);
        },
        else => @compileError("zml.exp10 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp10_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp10_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.exp10_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.exp10_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp10_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.exp10(x));
            },
            else => @compileError("zml.exp10_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.exp10_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn exp2(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.exp2(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp2 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.exp2(x);
        },
        else => @compileError("zml.exp2 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp2_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.exp2_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.exp2_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp2_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.exp2(x));
            },
            else => @compileError("zml.exp2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.exp2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn exp10m1(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.exp10m1(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp10m1 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.exp10m1(x);
        },
        else => @compileError("zml.exp10m1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp10m1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp10m1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.exp10m1_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.exp10m1_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp10m1_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.exp10m1(x));
            },
            else => @compileError("zml.exp10m1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.exp10m1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn exp2m1(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.exp2m1(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp2m1 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.exp2m1(x);
        },
        else => @compileError("zml.exp2m1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp2m1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp2m1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.exp2m1_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.exp2m1_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp2m1_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.exp2m1(x));
            },
            else => @compileError("zml.exp2m1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.exp2m1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn expm1(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.expm1(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.expm1 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.expm1(x);
        },
        else => @compileError("zml.expm1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn expm1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.expm1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.expm1_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.expm1_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.expm1_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.expm1(x));
            },
            else => @compileError("zml.expm1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.expm1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn log(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.log(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.log(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.log(x);
        },
        else => @compileError("zml.log not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.log_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.log_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.log(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.log(x));
            },
            else => @compileError("zml.log_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.log_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn log10(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.log10(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log10 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.log10(x);
        },
        else => @compileError("zml.log10 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log10_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log10_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.log10_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.log10_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log10_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.log10(x));
            },
            else => @compileError("zml.log10_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.log10_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn log2(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.log2(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log2 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.log2(x);
        },
        else => @compileError("zml.log2 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

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

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.log2_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.log2_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log2_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.log2(x));
            },
            else => @compileError("zml.log2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.log2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn log10p1(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.log10p1(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log10p1 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.log10p1(x);
        },
        else => @compileError("zml.log10p1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log10p1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log10p1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.log10p1_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.log10p1_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log10p1_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.log10p1(x));
            },
            else => @compileError("zml.log10p1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.log10p1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn log2p1(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.log2p1(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log2p1 not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.log2p1(x);
        },
        else => @compileError("zml.log2p1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log2p1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log2p1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.log2p1_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.log2p1_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log2p1_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.log2p1(x));
            },
            else => @compileError("zml.log2p1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.log2p1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn log1p(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.log1p(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log1p not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.log1p(x);
        },
        else => @compileError("zml.log1p not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log1p_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log1p_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.log1p_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.log1p_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log1p_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.log1p(x));
            },
            else => @compileError("zml.log1p_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.log1p_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

// Power functions
pub inline fn pow(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.pow(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => @compileError("zml.pow between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.pow(x, y);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.pow(x, y);
        },
        else => @compileError("zml.pow between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn pow_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.pow_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting

                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.pow_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.pow_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.pow_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.pow(x, y));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.pow(x, y));
            },
            else => @compileError("zml.pow_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.pow_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn sqrt(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.sqrt(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sqrt not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.sqrt(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.sqrt(x);
        },
        else => @compileError("zml.sqrt not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn sqrt_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sqrt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.sqrt_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.sqrt_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sqrt_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.sqrt(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.sqrt(x));
            },
            else => @compileError("zml.sqrt_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.sqrt_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn cbrt(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.cbrt(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cbrt not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.cbrt(x);
        },
        else => @compileError("zml.cbrt not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cbrt_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cbrt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.cbrt_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.cbrt_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.cbrt_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.cbrt(x));
            },
            else => @compileError("zml.cbrt_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.cbrt_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn hypot(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.hypot(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.hypot(x, y);
        },
        else => @compileError("zml.hypot between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn hypot_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.hypot_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting

                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.hypot_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.hypot_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.hypot_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.hypot(x, y));
            },
            else => @compileError("zml.hypot_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.hypot_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

// Trigonometric functions
pub inline fn sin(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.sin(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sin not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.sin(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.sin(x);
        },
        else => @compileError("zml.sin not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

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

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.sin_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.sin_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sin_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.sin(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.sin(x));
            },
            else => @compileError("zml.sin_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.sin_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn cos(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.cos(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cos not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.cos(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.cos(x);
        },
        else => @compileError("zml.cos not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cos_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cos_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.cos_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.cos_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.cos_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.cos(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.cos(x));
            },
            else => @compileError("zml.cos_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.cos_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn tan(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.tan(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tan not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.tan(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.tan(x);
        },
        else => @compileError("zml.tan not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tan_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.tan_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.tan_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.tan_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.tan_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.tan(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.tan(x));
            },
            else => @compileError("zml.tan_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.tan_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn asin(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.asin(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asin not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.asin(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.asin(x);
        },
        else => @compileError("zml.asin not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asin_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.asin_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.asin_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.asin_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.asin_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.asin(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.asin(x));
            },
            else => @compileError("zml.asin_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.asin_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn acos(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.acos(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acos not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.acos(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.acos(x);
        },
        else => @compileError("zml.acos not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acos_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.acos_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.acos_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.acos_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.acos_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.acos(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.acos(x));
            },
            else => @compileError("zml.acos_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.acos_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn atan(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.atan(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atan not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.atan(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.atan(x);
        },
        else => @compileError("zml.atan not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atan_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atan_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.atan_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.atan_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.atan_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.atan(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.atan(x));
            },
            else => @compileError("zml.atan_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.atan_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn atan2(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.atan2(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.atan2(x, y);
        },
        else => @compileError("zml.atan2 between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn atan2_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atan2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting

                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.atan2_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.atan2_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.atan2_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.atan2(x, y));
            },
            else => @compileError("zml.atan2_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.atan2_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn sinpi(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.sinpi(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sinpi not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.sinpi(x);
        },
        else => @compileError("zml.sinpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn sinpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sinpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.sinpi_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.sinpi_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sinpi_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.sinpi(x));
            },
            else => @compileError("zml.sinpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.sinpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn cospi(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.cospi(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cospi not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.cospi(x);
        },
        else => @compileError("zml.cospi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cospi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cospi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.cospi_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.cospi_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.cospi_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.cospi(x));
            },
            else => @compileError("zml.cospi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.cospi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn tanpi(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.tanpi(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tanpi not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.tanpi(x);
        },
        else => @compileError("zml.tanpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tanpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.tanpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.tanpi_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.tanpi_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.tanpi_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.tanpi(x));
            },
            else => @compileError("zml.tanpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.tanpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn asinpi(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.asinpi(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asinpi not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.asinpi(x);
        },
        else => @compileError("zml.asinpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asinpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.asinpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.asinpi_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.asinpi_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.asinpi_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.asinpi(x));
            },
            else => @compileError("zml.asinpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.asinpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn acospi(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.acospi(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acospi not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.acospi(x);
        },
        else => @compileError("zml.acospi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acospi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.acospi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.acospi_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.acospi_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.acospi_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.acospi(x));
            },
            else => @compileError("zml.acospi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.acospi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn atanpi(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.atanpi(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atanpi not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.atanpi(x);
        },
        else => @compileError("zml.atanpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atanpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atanpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.atanpi_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.atanpi_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.atanpi_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.atanpi(x));
            },
            else => @compileError("zml.atanpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.atanpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn atan2pi(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
    {
        comptime if (types.isArbitraryPrecision(Numeric(C))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.atan2pi(
            ctx.array_allocator,
            x,
            y,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.atan2pi not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.atan2pi(x, y);
        },
        else => @compileError("zml.atan2pi between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn atan2pi_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atan2pi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Both types are arbitrary precision
                if (Numeric(O) == Numeric(C)) {
                    // Equal types: output can be used for the operations, needing
                    // only the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                    );
                } else {
                    // Different types: internal allocator is required to perform
                    // the operation at `x`'s precision, and then cast the result to
                    // the output with the output's allocator
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                // Only the output is arbitrary precision, so we need the output's
                // allocator to perform the casting

                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(C))) {
                // Only the input is arbitrary precision, so we need the internal
                // allocator to perform the operation at `x`'s precision
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.atan2pi_(o, x, y, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X) or types.isArray(Y) or types.isSlice(Y)) {
        @compileError("zml.atan2pi_: o must be an `Array` or slice if x or y is an `Array` or slice, got " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y));
    }

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.atan2pi_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.atan2pi(x, y));
            },
            else => @compileError("zml.atan2pi_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.atan2pi_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

// Hyperbolic functions
pub inline fn sinh(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.sinh(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sinh not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.sinh(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.sinh(x);
        },
        else => @compileError("zml.sinh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn sinh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sinh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.sinh_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.sinh_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sinh_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.sinh(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.sinh(x));
            },
            else => @compileError("zml.sinh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.sinh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn cosh(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.cosh(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cosh not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.cosh(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.cosh(x);
        },
        else => @compileError("zml.cosh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cosh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cosh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.cosh_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.cosh_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.cosh_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.cosh(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.cosh(x));
            },
            else => @compileError("zml.cosh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.cosh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn tanh(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.tanh(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tanh not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.tanh(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.tanh(x);
        },
        else => @compileError("zml.tanh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tanh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.tanh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.tanh_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.tanh_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.tanh_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.tanh(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.tanh(x));
            },
            else => @compileError("zml.tanh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.tanh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn asinh(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.asinh(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asinh not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.asinh(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.asinh(x);
        },
        else => @compileError("zml.asinh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asinh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.asinh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.asinh_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.asinh_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.asinh_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.asinh(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.asinh(x));
            },
            else => @compileError("zml.asinh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.asinh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn acosh(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.acosh(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acosh not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.acosh(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.acosh(x);
        },
        else => @compileError("zml.acosh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

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

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.acosh_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.acosh_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.acosh_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.acosh(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.acosh(x));
            },
            else => @compileError("zml.acosh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.acosh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn atanh(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.atanh(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atanh not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.atanh(x);
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return cfloat.atanh(x);
        },
        else => @compileError("zml.atanh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atanh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atanh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.atanh_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.atanh_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.atanh_ not defined for " ++ @typeName(X) ++ " input type"),
            .int, .float => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.atanh(x));
            },
            .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, cfloat.atanh(x));
            },
            else => @compileError("zml.atanh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        },
        else => @compileError("zml.atanh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

// Error and gamma functions
pub inline fn erf(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.erf(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.erf not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.erf(x);
        },
        else => @compileError("zml.erf not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn erf_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.erf_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.erf_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.erf_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool => @compileError("zml.erf_ not defined for " ++ @typeName(X) ++ " input type"),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            o.* = scast(O, float.erf(x));
        },
        else => @compileError("zml.erf_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn erfc(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.erfc(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.erfc not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.erfc(x);
        },
        else => @compileError("zml.erfc not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn erfc_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.erfc_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.erfc_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.erfc_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool => @compileError("zml.erfc_ not defined for " ++ @typeName(X) ++ " input type"),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            o.* = scast(O, float.erfc(x));
        },
        else => @compileError("zml.erfc_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn gamma(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.gamma(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.gamma not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.gamma(x);
        },
        else => @compileError("zml.gamma not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn gamma_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.gamma_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.gamma_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.gamma_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool => @compileError("zml.gamma_ not defined for " ++ @typeName(X) ++ " input type"),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            o.* = scast(O, float.gamma(x));
        },
        else => @compileError("zml.gamma_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn lgamma(
    x: anytype,
    ctx: anytype,
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.lgamma(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.lgamma not defined for " ++ @typeName(X)),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.lgamma(x);
        },
        else => @compileError("zml.lgamma not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn lgamma_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.lgamma_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.lgamma_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.lgamma_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool => @compileError("zml.lgamma_ not defined for " ++ @typeName(X) ++ " input type"),
        .int, .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            o.* = scast(O, float.lgamma(x));
        },
        else => @compileError("zml.lgamma_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn conjugate(
    x: anytype,
    ctx: anytype,
) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                    .copy = .{ .type = bool, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.conjugate(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.conjugate not defined for " ++ @typeName(X)),
        .int => return x,
        .float => return x,
        .cfloat => return x.conjugate(),
        else => @compileError("zml.conjugate not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

// Nearest integer operations
pub inline fn ceil(
    x: anytype,
    ctx: anytype,
) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    .order = .{ .type = ?array.Order, .required = false },
                },
            );
        };

        return array.ceil(
            ctx.array_allocator,
            x,
            .{ .order = getFieldOrDefault(ctx, "order", ?array.Order, null) },
            stripStruct(ctx, &.{ "array_allocator", "order" }),
        );
    }

    switch (types.numericType(X)) {
        .bool => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        .int => return x,
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return float.ceil(x);
        },
        .cfloat => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        else => @compileError("zml.ceil not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn ceil_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ceil_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.ceil_(o, x, ctx);
    } else if (comptime types.isArray(X) or types.isSlice(X)) {
        @compileError("zml.ceil_: o must be an `Array` or slice if x is an `Array` or slice, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    }

    switch (comptime types.numericType(X)) {
        .bool => @compileError("zml.ceil_ not defined for " ++ @typeName(X) ++ " input type"),
        .int => o.* = x,
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            o.* = scast(O, float.ceil(x));
        },
        .cfloat => @compileError("zml.ceil_ not defined for " ++ @typeName(X) ++ " input type"),
        else => @compileError("zml.ceil_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
    }
}

pub inline fn init(
    comptime T: type,
    ctx: anytype,
) !T {
    if (comptime types.isArray(T) or types.isSlice(T))
        @compileError("zml.init not implemented for arrays or slices yet.");

    switch (comptime types.numericType(T)) {
        .bool => {
            comptime validateContext(@TypeOf(ctx), .{});

            return false;
        },
        .int => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .float => {
            comptime validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return .{ .re = 0, .im = 0 };
        },
        else => @compileError("zml.init not implemented for " ++ @typeName(T) ++ " yet"),
    }
}

/// Sets the value of `o` to `x`.
///
/// When `o` is of an arbitrary precision type, its already allocated memory is
/// used, evading a new allocation (reallocation may be needed if more space is
/// needed). For `o` a fixed precision type, this is equivalent to `cast`.
pub inline fn set(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.set requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        @compileError("zml.set not implemented for arrays or slices yet.");

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
            .bool, .int, .float, .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .integer => switch (comptime types.numericType(X)) {
            .bool, .int, .float, .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like integer.set
        },
        .rational => switch (comptime types.numericType(X)) {
            .bool, .int, .float, .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like rational.set
        },
        .real => switch (comptime types.numericType(X)) {
            .bool, .int, .float, .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like real.set
        },
        .complex => switch (comptime types.numericType(X)) {
            .bool, .int, .float, .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like complex.set
        },
        .expression => switch (comptime types.numericType(X)) {
            .bool, .int, .float, .cfloat => {
                comptime validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, x);
            },
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like expression.set
        },
    }
}

pub inline fn copy(
    x: anytype,
    ctx: anytype,
) @TypeOf(x) {
    comptime if (!types.isPointer(@TypeOf(x)) or types.isConstPointer(@TypeOf(x)))
        @compileError("zml.copy requires x a mutable pointer, got " ++ @typeName(@TypeOf(x)));

    const X: type = Child(@TypeOf(x));

    if (comptime types.isArray(X) or types.isSlice(X))
        @compileError("zml.copy not implemented for arrays or slices yet.");

    switch (types.numericType(X)) {
        .bool, .int, .float, .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            return x.*;
        },
        else => @compileError("zml.copy not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn deinit(
    x: anytype,
    ctx: anytype,
) void {
    comptime if (!types.isPointer(@TypeOf(x)) or types.isConstPointer(@TypeOf(x)))
        @compileError("zml.deinit requires x a mutable pointer, got " ++ @typeName(@TypeOf(x)));

    const X: type = Child(@TypeOf(x));

    if (comptime types.isArray(X) or types.isSlice(X))
        @compileError("zml.deinit not implemented for arrays or slices yet.");

    switch (types.numericType(X)) {
        .bool, .int, .float, .cfloat => {
            comptime validateContext(@TypeOf(ctx), .{});

            // No deinitialization needed for fixed precision types, this is a no-op.
            return;
        },
        else => @compileError("zml.deinit not implemented for " ++ @typeName(X) ++ " yet"),
    }
}
