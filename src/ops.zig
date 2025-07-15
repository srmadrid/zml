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
/// options (`struct`): Options for the addition operation.
/// - mode (`int.Mode`): The mode of the addition operation. Only needed when
/// adding two `int` values.
/// - allocator (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision or an `Array`.
/// - writeable (`bool`): Whether the output should be writeable. Only needed if
/// the output type is an `Array`.
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
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.add(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.add(x, y, .{ .mode = options.mode }),
        .float => return float.add(x, y),
        .cfloat => return cfloat.add(x, y),
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
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.add_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.add_(o, x, y, .{ .allocator = options.allocator, .mode = options.mode });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.add_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.add(x, y, .{ .mode = options.mode })),
            .float => o.* = scast(O, float.add(x, y)),
            .cfloat => o.* = scast(O, cfloat.add(x, y)),
            else => @compileError("zml.add_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        //.integer => switch (comptime types.numericType(C)) {
        //    .bool => @compileError("zml.add_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
        //    .int => try set(o, int.add(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
        //    .float => try set(o, float.add(x, y), .{ .allocator = options.allocator }),
        //    .cfloat => try set(o, cfloat.add(x, y), .{ .allocator = options.allocator }),
        //    else => @compileError("zml.add_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        //},
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
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.sub(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.sub(x, y, .{ .mode = options.mode }),
        .float => return float.sub(x, y),
        .cfloat => return cfloat.sub(x, y),
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
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sub_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.sub_(o, x, y, .{ .allocator = options.allocator, .mode = options.mode });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.sub_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.sub(x, y, .{ .mode = options.mode })),
            .float => o.* = scast(O, float.sub(x, y)),
            .cfloat => o.* = scast(O, cfloat.sub(x, y)),
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
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.mul(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.mul(x, y, .{ .mode = options.mode }),
        .float => return float.mul(x, y),
        .cfloat => return cfloat.mul(x, y),
        else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn mul_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.mul_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.mul_(o, x, y, .{ .allocator = options.allocator, .mode = options.mode });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.mul_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.mul(x, y, .{ .mode = options.mode })),
            .float => o.* = scast(O, float.mul(x, y)),
            .cfloat => o.* = scast(O, cfloat.mul(x, y)),
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
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.div(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.div not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.div(x, y),
        .float => return float.div(x, y),
        .cfloat => return cfloat.div(x, y),
        else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn div_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.div_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.div_(o, x, y, .{ .allocator = options.allocator });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.div_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
        .int => try set(o, int.div(x, y), .{ .allocator = options.allocator }),
        .float => try set(o, float.div(x, y), .{ .allocator = options.allocator }),
        .cfloat => try set(o, cfloat.div(x, y), .{ .allocator = options.allocator }),
        else => @compileError("zml.div_ not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn eq(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.eq(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => return x == y,
        .int => return int.eq(x, y),
        .float => return float.eq(x, y),
        .cfloat => return cfloat.eq(x, y),
        else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn eq_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.eq_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.eq_(o, x, y, .{ .allocator = options.allocator });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.eq_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.eq(x, y)),
            .float => o.* = scast(O, float.eq(x, y)),
            .cfloat => o.* = scast(O, cfloat.eq(x, y)),
            else => @compileError("zml.eq_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.eq_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn ne(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.ne(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => return x != y,
        .int => return int.ne(x, y),
        .float => return float.ne(x, y),
        .cfloat => return cfloat.ne(x, y),
        else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn ne_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ne_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.ne_(o, x, y, .{ .allocator = options.allocator });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.ne_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.ne(x, y)),
            .float => o.* = scast(O, float.ne(x, y)),
            .cfloat => o.* = scast(O, cfloat.ne(x, y)),
            else => @compileError("zml.ne_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.ne_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn lt(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.lt(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.lt(x, y),
        .float => return float.lt(x, y),
        .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn lt_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.lt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.lt_(o, x, y, .{ .allocator = options.allocator });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.lt_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.lt(x, y)),
            .float => o.* = scast(O, float.lt(x, y)),
            .cfloat => o.* = scast(O, cfloat.lt(x, y)),
            else => @compileError("zml.lt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.lt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn le(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.le(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.le(x, y),
        .float => return float.le(x, y),
        .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn le_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.le_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.le_(o, x, y, .{ .allocator = options.allocator });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.le_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.le(x, y)),
            .float => o.* = scast(O, float.le(x, y)),
            .cfloat => o.* = scast(O, cfloat.le(x, y)),
            else => @compileError("zml.le_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.le_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn gt(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.gt(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.gt(x, y),
        .float => return float.gt(x, y),
        .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn gt_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.gt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.gt_(o, x, y, .{ .allocator = options.allocator });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.gt_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.gt(x, y)),
            .float => o.* = scast(O, float.gt(x, y)),
            .cfloat => o.* = scast(O, cfloat.gt(x, y)),
            else => @compileError("zml.gt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.gt_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn ge(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.ge(options.allocator.?, x, y, .{ .writeable = true });

    switch (comptime types.numericType(C)) {
        .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => return int.ge(x, y),
        .float => return float.ge(x, y),
        .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn ge_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ge_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.ge_(o, x, y, .{ .allocator = options.allocator });

    switch (comptime types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
            .bool => @compileError("zml.ge_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = scast(O, int.ge(x, y)),
            .float => o.* = scast(O, float.ge(x, y)),
            .cfloat => o.* = scast(O, cfloat.ge(x, y)),
            else => @compileError("zml.ge_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.ge_ not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn max(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.max not implemented for arrays or slices yet.");
    // return array.max(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.max not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @max(x, y),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @max(x, y),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
    }
}

pub inline fn min(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.min not implemented for arrays or slices yet.");
    // return array.min(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.min not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @min(x, y),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @min(x, y),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
    }
}

// Basic operations
pub inline fn abs(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
        copy: bool = true, // If arbitrary precision, copy the value instead of returning a view (just editing `.positive = true`). For copy = false, allocator is not needed.
    },
) !CoerceToArray(@TypeOf(x), Scalar(Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.abs(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs not defined for " ++ @typeName(X)),
        .int => return int.abs(x),
        .float => return float.abs(x),
        .cfloat => return cfloat.abs(x),
        else => @compileError("zml.abs not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn abs_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.abs_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => try set(o, int.abs(x), .{ .allocator = options.allocator }),
        .float => try set(o, float.abs(x), .{ .allocator = options.allocator }),
        .cfloat => try set(o, cfloat.abs(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.abs_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn abs2(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.abs2(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs2 not defined for " ++ @typeName(X)),
        .int => return int.mul(x, x),
        .float => return float.pow(x, 2),
        .cfloat => return cfloat.abs2(x),
        else => @compileError("zml.abs2 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn abs2_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.abs2_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, int.mul(x, x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.pow(x, 2), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.abs2(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.abs2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

// Exponential functions
pub inline fn exp(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.exp(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp not defined for " ++ @typeName(X)),
        .int, .float => return float.exp(x),
        .cfloat => return cfloat.exp(x),
        else => @compileError("zml.exp not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.exp_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.exp(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.exp(x), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.exp(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.exp_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp10(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.exp10(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp10 not defined for " ++ @typeName(X)),
        .int, .float => return float.exp10(x),
        else => @compileError("zml.exp10 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp10_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp10_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.exp10_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp10_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.exp10(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.exp10(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.exp10_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp2(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.exp2(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp2 not defined for " ++ @typeName(X)),
        .int, .float => return float.exp2(x),
        else => @compileError("zml.exp2 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp2_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.exp2_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.exp2(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.exp2(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.exp2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp10m1(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.exp10m1(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp10m1 not defined for " ++ @typeName(X)),
        .int, .float => return float.exp10m1(x),
        else => @compileError("zml.exp10m1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp10m1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp10m1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.exp10m1_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp10m1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.exp10m1(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.exp10m1(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.exp10m1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp2m1(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.exp2m1(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp2m1 not defined for " ++ @typeName(X)),
        .int, .float => return float.exp2m1(x),
        else => @compileError("zml.exp2m1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn exp2m1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.exp2m1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.exp2m1_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.exp2m1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.exp2m1(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.exp2m1(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.exp2m1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn expm1(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.expm1(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.expm1 not defined for " ++ @typeName(X)),
        .int, .float => return float.expm1(x),
        else => @compileError("zml.expm1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn expm1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.expm1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.expm1_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.expm1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.expm1(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.expm1(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.expm1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.log(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log not defined for " ++ @typeName(X)),
        .int, .float => return float.log(x),
        .cfloat => return cfloat.log(x),
        else => @compileError("zml.log not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.log_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.log(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.log(x), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.log(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.log_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log10(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.log10(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log10 not defined for " ++ @typeName(X)),
        .int, .float => return float.log10(x),
        else => @compileError("zml.log10 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log10_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log10_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.log10_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log10_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.log10(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.log10(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.log10_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log2(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.log2(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log2 not defined for " ++ @typeName(X)),
        .int, .float => return float.log2(x),
        else => @compileError("zml.log2 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log2_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.log2_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.log2(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.log2(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.log2_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log10p1(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.log10p1(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log10p1 not defined for " ++ @typeName(X)),
        .int, .float => return float.log10p1(x),
        else => @compileError("zml.log10p1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log10p1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log10p1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.log10p1_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log10p1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.log10p1(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.log10p1(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.log10p1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log2p1(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.log2p1(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log2p1 not defined for " ++ @typeName(X)),
        .int, .float => return float.log2p1(x),
        else => @compileError("zml.log2p1 not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log2p1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log2p1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.log2p1_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log2p1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.log2p1(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.log2p1(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.log2p1_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log1p(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.log1p(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log1p not defined for " ++ @typeName(X)),
        .int, .float => return float.log1p(x),
        else => @compileError("zml.log1p not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn log1p_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.log1p_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.log1p_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.log1p_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.log1p(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.log1p(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.log1p_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

// Power functions
pub inline fn pow(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.pow(options.allocator.?, x, y, .{ .writeable = true });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int, .float => return float.pow(x, y),
        .cfloat => return cfloat.pow(x, y),
        else => @compileError("zml.pow between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn pow_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.pow_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.pow_(o, x, y, .{ .allocator = options.allocator });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.pow_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
        .int, .float => o.* = try cast(O, float.pow(x, y), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.pow(x, y), .{ .allocator = options.allocator }),
        else => @compileError("zml.pow_ not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn sqrt(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.sqrt(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sqrt not defined for " ++ @typeName(X)),
        .int, .float => return float.sqrt(x),
        .cfloat => return cfloat.sqrt(x),
        else => @compileError("zml.sqrt not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn sqrt_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sqrt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.sqrt_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sqrt_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.sqrt(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.sqrt(x), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.sqrt(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.sqrt_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cbrt(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.cbrt(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cbrt not defined for " ++ @typeName(X)),
        .int, .float => return float.cbrt(x),
        else => @compileError("zml.cbrt not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cbrt_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cbrt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.cbrt_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cbrt_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.cbrt(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.cbrt(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.cbrt_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn hypot(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.hypot(options.allocator.?, x, y, .{ .writeable = options.writeable });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int, .float => return float.hypot(x, y),
        else => @compileError("zml.hypot between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn hypot_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.hypot_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.hypot_(o, x, y, .{ .allocator = options.allocator });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.hypot_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
        .int, .float => o.* = try cast(O, float.hypot(x, y), .{ .allocator = options.allocator }),
        else => @compileError("zml.hypot_ not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

// Trigonometric functions
pub inline fn sin(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.sin(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sin not defined for " ++ @typeName(X)),
        .int, .float => return float.sin(x),
        .cfloat => return cfloat.sin(x),
        else => @compileError("zml.sin not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn sin_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sin_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.sin_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sin_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.sin(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.sin(x), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.sin(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.sin_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cos(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.cos(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cos not defined for " ++ @typeName(X)),
        .int, .float => return float.cos(x),
        .cfloat => return cfloat.cos(x),
        else => @compileError("zml.cos not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cos_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cos_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.cos_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cos_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.cos(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.cos(x), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.cos(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.cos_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tan(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.tan(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tan not defined for " ++ @typeName(X)),
        .int, .float => return float.tan(x),
        .cfloat => return cfloat.tan(x),
        else => @compileError("zml.tan not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tan_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.tan_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.tan_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tan_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.tan(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.tan(x), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.tan(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.tan_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asin(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.asin(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asin not defined for " ++ @typeName(X)),
        .int, .float => return float.asin(x),
        else => @compileError("zml.asin not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asin_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.asin_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.asin_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asin_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.asin(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.asin(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.asin_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acos(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.acos(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acos not defined for " ++ @typeName(X)),
        .int, .float => return float.acos(x),
        else => @compileError("zml.acos not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acos_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.acos_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.acos_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acos_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.acos(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.acos(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.acos_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atan(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.atan(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atan not defined for " ++ @typeName(X)),
        .int, .float => return float.atan(x),
        else => @compileError("zml.atan not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atan_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atan_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.atan_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atan_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.atan(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.atan(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.atan_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atan2(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.atan2(options.allocator.?, x, y, .{ .writeable = options.writeable });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int, .float => return float.atan2(x, y),
        else => @compileError("zml.atan2 between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn atan2_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atan2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.atan2_(o, x, y, .{ .allocator = options.allocator });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.atan2_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
        .int, .float => o.* = try cast(O, float.atan2(x, y), .{ .allocator = options.allocator }),
        else => @compileError("zml.atan2_ not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

pub inline fn sinpi(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.sinpi(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sinpi not defined for " ++ @typeName(X)),
        .int, .float => return float.sinpi(x),
        else => @compileError("zml.sinpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn sinpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sinpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.sinpi_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sinpi_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.sinpi(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.sinpi(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.sinpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cospi(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.cospi(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cospi not defined for " ++ @typeName(X)),
        .int, .float => return float.cospi(x),
        else => @compileError("zml.cospi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cospi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cospi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.cospi_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cospi_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.cospi(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.cospi(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.cospi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tanpi(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.tanpi(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tanpi not defined for " ++ @typeName(X)),
        .int, .float => return float.tanpi(x),
        else => @compileError("zml.tanpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tanpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.tanpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.tanpi_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tanpi_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.tanpi(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.tanpi(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.tanpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asinpi(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.asinpi(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asinpi not defined for " ++ @typeName(X)),
        .int, .float => return float.asinpi(x),
        else => @compileError("zml.asinpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asinpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.asinpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.asinpi_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asinpi_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.asinpi(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.asinpi(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.asinpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acospi(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.acospi(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acospi not defined for " ++ @typeName(X)),
        .int, .float => return float.acospi(x),
        else => @compileError("zml.acospi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acospi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.acospi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.acospi_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acospi_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.acospi(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.acospi(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.acospi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atanpi(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.atanpi(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atanpi not defined for " ++ @typeName(X)),
        .int, .float => return float.atanpi(x),
        else => @compileError("zml.atanpi not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atanpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atanpi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.atanpi_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atanpi_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.atanpi(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.atanpi(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.atanpi_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atan2pi(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.atan2pi(options.allocator.?, x, y, .{ .writeable = options.writeable });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.atan2pi not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int, .float => return float.atan2pi(x, y),
        else => @compileError("zml.atan2pi between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

pub inline fn atan2pi_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atan2pi_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        return array.atan2pi_(o, x, y, .{ .allocator = options.allocator });

    switch (types.numericType(C)) {
        .bool => @compileError("zml.atan2pi_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
        .int, .float => o.* = try cast(O, float.atan2pi(x, y), .{ .allocator = options.allocator }),
        else => @compileError("zml.atan2pi_ not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
    }
}

// Hyperbolic functions
pub inline fn sinh(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.sinh(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sinh not defined for " ++ @typeName(X)),
        .int, .float => return float.sinh(x),
        else => @compileError("zml.sinh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn sinh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sinh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.sinh_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.sinh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.sinh(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.sinh(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.sinh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cosh(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.cosh(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cosh not defined for " ++ @typeName(X)),
        .int, .float => return float.cosh(x),
        else => @compileError("zml.cosh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn cosh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.cosh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.cosh_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.cosh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.cosh(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.cosh(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.cosh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tanh(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.tanh(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tanh not defined for " ++ @typeName(X)),
        .int, .float => return float.tanh(x),
        else => @compileError("zml.tanh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn tanh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.tanh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.tanh_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.tanh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.tanh(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.tanh(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.tanh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asinh(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.asinh(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asinh not defined for " ++ @typeName(X)),
        .int, .float => return float.asinh(x),
        else => @compileError("zml.asinh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn asinh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.asinh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.asinh_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.asinh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.asinh(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.asinh(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.asinh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acosh(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.acosh(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acosh not defined for " ++ @typeName(X)),
        .int, .float => return float.acosh(x),
        else => @compileError("zml.acosh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn acosh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.acosh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.acosh_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.acosh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.acosh(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.acosh(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.acosh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atanh(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.atanh(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atanh not defined for " ++ @typeName(X)),
        .int, .float => return float.atanh(x),
        else => @compileError("zml.atanh not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn atanh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atanh_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.atanh_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.atanh_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.atanh(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.atanh(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.atanh_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

// Error and gamma functions
pub inline fn erf(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.erf(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.erf not defined for " ++ @typeName(X)),
        .int, .float => return float.erf(x),
        else => @compileError("zml.erf not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn erf_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.erf_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.erf_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.erf_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.erf(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.erf(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.erf_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn erfc(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.erfc(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.erfc not defined for " ++ @typeName(X)),
        .int, .float => return float.erfc(x),
        else => @compileError("zml.erfc not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn erfc_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.erfc_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.erfc_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.erfc_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.erfc(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.erfc(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.erfc_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn gamma(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.gamma(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.gamma not defined for " ++ @typeName(X)),
        .int, .float => return float.gamma(x),
        else => @compileError("zml.gamma not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn gamma_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.gamma_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.gamma_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.gamma_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.gamma(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.gamma(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.gamma_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn lgamma(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !EnsureFloat(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.lgamma(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.lgamma not defined for " ++ @typeName(X)),
        .int, .float => return float.lgamma(x),
        else => @compileError("zml.lgamma not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn lgamma_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.lgamma_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.lgamma_(o, x, .{ .allocator = options.allocator });

    switch (types.numericType(X)) {
        .bool => @compileError("zml.lgamma_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, float.lgamma(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.lgamma(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.lgamma_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn conjugate(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
        copy: bool = true, // If arbitrary precision, copy the value instead of returning a view (just editing `.positive = true`). For copy = false, allocator is not needed.
    },
) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.conjugate(options.allocator.?, x, .{ .writeable = options.writeable });

    switch (types.numericType(X)) {
        .bool, .int, .float => return x,
        .cfloat => return x.conjugate(),
        else => @compileError("zml.conjugate not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

// Nearest integer operations
pub inline fn ceil(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        return array.ceil(options.allocator.?, x, .{});

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        .int => return x,
        .float => return float.ceil(x),
        .cfloat => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        .complex => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        else => @compileError("zml.ceil not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn ceil_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ceil_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.ceil_(o, x, .{ .allocator = options.allocator });

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ceil_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        .int => o.* = try cast(O, int.ceil(x), .{ .allocator = options.allocator }),
        .float => o.* = try cast(O, float.ceil(x), .{ .allocator = options.allocator }),
        .cfloat => o.* = try cast(O, cfloat.ceil(x), .{ .allocator = options.allocator }),
        else => @compileError("zml.ceil_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn init(
    comptime T: type,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !T {
    if (comptime types.isArray(T) or types.isSlice(T))
        @compileError("zml.init not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (comptime types.numericType(T)) {
        .bool => return false,
        .int => return 0,
        .float => return 0,
        .cfloat => return .{ .re = 0, .im = 0 },
        else => @compileError("zml.init not implemented for " ++ @typeName(T) ++ " yet"),
    }
}

/// Sets the value of `o` to `x`, using the specified allocator if needed.
///
/// When `o` is of an arbitrary precision type, its already allocated memory is
/// used, evading a new allocation (reallocation may be needed if more space is
/// needed). For `o` a fixed precision type, this is equivalent to `cast`.
pub inline fn set(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.set requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        @compileError("zml.set not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(O)) {
        .bool, .int, .float, .cfloat => switch (types.numericType(X)) {
            .bool, .int, .float, .cfloat => o.* = try cast(O, x, .{ .allocator = options.allocator }),
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .integer => switch (types.numericType(X)) {
            .bool, .int, .float, .cfloat => o.* = try cast(O, x, .{ .allocator = options.allocator }),
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like integer.set
        },
        .rational => switch (types.numericType(X)) {
            .bool, .int, .float, .cfloat => o.* = try cast(O, x, .{ .allocator = options.allocator }),
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like rational.set
        },
        .real => switch (types.numericType(X)) {
            .bool, .int, .float, .cfloat => o.* = try cast(O, x, .{ .allocator = options.allocator }),
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like real.set
        },
        .complex => switch (types.numericType(X)) {
            .bool, .int, .float, .cfloat => o.* = try cast(O, x, .{ .allocator = options.allocator }),
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like complex.set
        },
        .expression => switch (types.numericType(X)) {
            .bool, .int, .float, .cfloat => o.* = try cast(O, x, .{ .allocator = options.allocator }),
            else => @compileError("zml.set not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"), // Something like expression.set
        },
    }
}

pub inline fn copy(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) @TypeOf(x) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X) or types.isSlice(X))
        @compileError("zml.copy not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool, .int, .float, .cfloat => return x,
        else => @compileError("zml.copy not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn deinit(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) void {
    comptime if (!types.isPointer(@TypeOf(x)) or types.isConstPointer(@TypeOf(x)))
        @compileError("zml.deinit requires x a mutable pointer, got " ++ @typeName(@TypeOf(x)));

    const X: type = Child(@TypeOf(x));

    if (comptime types.isArray(X) or types.isSlice(X))
        @compileError("zml.deinit not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool, .int, .float, .cfloat => {},
        else => @compileError("zml.deinit not implemented for " ++ @typeName(X) ++ " yet"),
    }
}
