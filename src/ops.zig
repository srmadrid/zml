const std = @import("std");
const types = @import("types.zig");
const cast = types.cast;
const Scalar = types.Scalar;
const Coerce = types.Coerce;
const CoerceToArray = types.CoerceToArray;
const Child = types.Child;
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
/// - `mode` (`int.Mode`): The mode of the addition operation. Only needed when
/// adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision.
///
/// Returns
/// -------
/// `Coerce(@TypeOf(x), @TypeOf(y))`: The result of the addition operation.
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
/// or if the types of the inputs are not supported numeric types, `Array`s or
/// slices.
///
/// See Also
/// --------
/// `add_`: For in-place addition of two values.
///
/// `add_to`: For addition of two values and storing the result in the output
/// pointer.
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
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = types.Coerce(X, Y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.add not implemented for arrays or slices yet.");
    // return array.add(x, y, options.allocator.?);

    switch (types.numericType(C)) {
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
///     o += y
/// ```
///
/// Parameters
/// ----------
/// o (pointer to `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The output pointer where the
/// result of the addition operation will be stored, and the left-hand side
/// operand.
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
/// `add_to`: For addition of two values and storing the result in the output
/// pointer.
///
/// `int.add_`: For in-place addition of an `int` and another `int` or `bool`.
///
/// `float.add_`: For in-place addition of a `float` and another `float`, `int`
/// or `bool`.
///
/// `cfloat.add_`: For in-place addition of a `cfloat` and another `cfloat`,
/// `float`, `int` or `bool`.
///
/// `integer.add_`: For in-place addition of an `integer` and another `integer`,
/// `bool` or `int`.
///
/// `rational.add_`: For in-place addition of a `rational` and another
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `real.add_`: For in-place addition of a `real` and another `real`,
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `complex.add_`: For in-place addition of a `complex` and another `complex`,
/// `real`, `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.add_`: For in-place addition of an `expression` and another
/// `expression`, `complex`, `real`, `rational`, `integer`, `float`, `int` or
/// `bool`.
///
/// `array.add_`: For in-place addition of two `Array`s or slices.
///
/// Notes
/// -----
/// The addition is performed with the precision of the output type.
pub inline fn add_(
    o: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.add_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.add_ not implemented for arrays or slices yet.");
    // return array.add(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ " output type"),
        .int => return int.add_(o, y, .{ .mode = options.mode }),
        .float => return float.add_(o, y),
        .cfloat => return cfloat.add_(o, y),
        else => @compileError("zml.add_ not implemented for " ++ @typeName(O) ++ " output type"),
    }
}

/// Adds two values of any two supported types storing the result in the output
/// pointer in-place.
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
/// result of the addition operation will be stored.
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
/// `add_`: For in-place addition of two values.
///
/// Notes
/// -----
/// The addition is performed with the precision of the coerced type of the
/// inputs, and the result is cast to the output type if necessary. This cast
/// is not checked for safety.
pub inline fn add_to(
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
        @compileError("zml.add_to requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.add_to not implemented for arrays or slices yet.");
    // return array.add_to(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => switch (types.numericType(C)) {
            .bool => @compileError("zml.add_to not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = try cast(O, int.add(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.add(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.add(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.add_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .int => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.add(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.add(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.add(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.add_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .float => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.add(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.add(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.add(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.add_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .cfloat => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.add(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.add(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.add(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.add_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.add_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
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
/// `sub_to`: For subtraction of two values and storing the result in the
/// output pointer.
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
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.sub not implemented for arrays or slices yet.");
    // return array.sub(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.sub(x, y, .{ .mode = options.mode }),
                .float => return float.sub(x, y),
                .cfloat => return cfloat.sub(x, y),
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.sub(x, y, .{ .mode = options.mode }),
                .float => return float.sub(x, y),
                .cfloat => return cfloat.sub(x, y),
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.sub(x, y),
                .cfloat => return cfloat.sub(x, y),
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.sub(x, y),
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

/// Subtracts two values of any two supported types in-place.
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
/// result of the addition operation will be stored, and the left-hand side
/// operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the subtraction operation.
/// - `mode` (`int.Mode`): The mode of the subtraction operation. Only needed
/// when adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision. May not be used if the output has enouph memory allocated
/// already.
///
/// Returns
/// -------
/// `void`: The result of the subtraction operation is stored in the output
/// pointer.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision.
/// `array.Error.Bla`: Put errors when `array.sub_` is implemented.
///
/// Raises
/// ------
/// `@compileError`: If the subtraction operation is not defined for the types
/// of the inputs, if the types of the inputs cannot be coerced to a common
/// type, if the types of the inputs are not supported numeric types, `Array`s
/// or slices, if the output pointer is not a mutable pointer, if the output's
/// child type is not a supported numeric type, `Array` or slice, or if the
/// coerced type is an `Array` and the output's child type is not an `Array`.
///
/// See Also
/// --------
/// `sub`: For subtraction of two values and returning the result.
///
/// `sub_to`: For subtraction of two values and storing the result in the
/// output pointer.
///
/// `int.sub_`: For in-place subtraction of an `int` and another `int` or
/// `bool`.
///
/// `float.sub_`: For in-place subtraction of a `float` and another `float`,
/// `int` or `bool`.
///
/// `cfloat.sub_`: For in-place subtraction of a `cfloat` and another `cfloat`,
/// `float`, `int` or `bool`.
///
/// `integer.sub_`: For in-place subtraction of an `integer` and another
/// `integer`, `bool` or `int`.
///
/// `rational.sub_`: For in-place subtraction of a `rational` and another
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `real.sub_`: For in-place subtraction of a `real` and another `real`,
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `complex.sub_`: For in-place subtraction of a `complex` and another
/// `complex`, `real`, `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.sub_`: For in-place subtraction of an `expression` and another
/// `expression`, `complex`, `real`, `rational`, `integer`, `float`, `int` or
/// `bool`.
///
/// `array.sub_`: For in-place subtraction of two `Array`s or slices.
///
/// Notes
/// -----
/// The subtraction is performed with the precision of the output type.
pub inline fn sub_(
    o: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sub_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.sub_ not implemented for arrays or slices yet.");
    // return array.add(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => @compileError("zml.sub_ not defined for " ++ @typeName(O) ++ " output type"),
        .int => return int.sub_(o, y, .{ .mode = options.mode }),
        .float => return float.sub_(o, y),
        .cfloat => return cfloat.sub_(o, y),
        else => @compileError("zml.sub_ not implemented for " ++ @typeName(O) ++ " output type"),
    }
}

/// Sibtracts two values of any two supported types storing the result in the
/// output pointer in-place.
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
/// result of the subtraction operation will be stored.
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
/// `add_`: For in-place addition of two values.
///
/// Notes
/// -----
/// The addition is performed with the precision of the coerced type of the
/// inputs, and the result is cast to the output type if necessary. This cast
/// is not checked for safety.
pub inline fn sub_to(
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
        @compileError("zml.sub_to requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.sub_to not implemented for arrays or slices yet.");
    // return array.sub_to(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => switch (types.numericType(C)) {
            .bool => @compileError("zml.sub_to not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = try cast(O, int.sub(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.sub(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.sub(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.sub_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .int => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.sub(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.sub(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.sub(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.sub_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .float => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.sub(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.sub(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.sub(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.sub_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .cfloat => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.sub(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.sub(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.sub(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.sub_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.sub_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
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
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.mul not implemented for arrays or slices yet.");
    // return array.mul(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.mul(x, y, .{ .mode = options.mode }),
                .float => return float.mul(x, y),
                .cfloat => return cfloat.mul(x, y),
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.mul(x, y, .{ .mode = options.mode }),
                .float => return float.mul(x, y),
                .cfloat => return cfloat.mul(x, y),
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.mul(x, y),
                .cfloat => return cfloat.mul(x, y),
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.mul(x, y),
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

/// Multiplies two values of any two supported types in-place.
///
/// The function supports multiplication for values of any combination of
/// supported numeric types, `Array`s and slices, and stores the result in the
/// output pointer. The operation performed is
///
/// ```zig
///     o = x * y
/// ```
///
/// Parameters
/// ----------
/// o (pointer to `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The output pointer where the
/// result of the addition operation will be stored, and the left-hand side
/// operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the multiplication operation.
/// - `mode` (`int.Mode`): The mode of the multiplication operation. Only needed
/// when adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision. May not be used if the output has enouph memory allocated
/// already.
///
/// Returns
/// -------
/// `void`: The result of the multiplication operation is stored in the output
/// pointer.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value. Only occurs if the output type is of arbitrary
/// precision.
/// `array.Error.Bla`: Put errors when `array.mul_` is implemented.
///
/// Raises
/// ------
/// `@compileError`: If the multiplication operation is not defined for the
/// types of the inputs, if the types of the inputs cannot be coerced to a
/// common type, if the types of the inputs are not supported numeric types,
/// `Array`s or slices, if the output pointer is not a mutable pointer, if the
/// output's child type is not a supported numeric type, `Array` or slice, or if
/// the coerced type is an `Array` and the output's child type is not an
/// `Array`.
///
/// See Also
/// --------
/// `mul`: For multiplication of two values and returning the result.
///
/// `int.mul_`: For in-place multiplication of an `int` and another `int` or
/// `bool`.
///
/// `float.mul_`: For in-place multiplication of a `float` and another `float`,
/// `int` or `bool`.
///
/// `cfloat.mul_`: For in-place multiplication of a `cfloat` and another
/// `cfloat`, `float`, `int` or `bool`.
///
/// `integer.mul_`: For in-place multiplication of an `integer` and another
/// `integer`, `bool` or `int`.
///
/// `rational.mul_`: For in-place multiplication of a `rational` and another
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `real.mul_`: For in-place multiplication of a `real` and another `real`,
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `complex.mul_`: For in-place multiplication of a `complex` and another
/// `complex`, `real`, `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.mul_`: For in-place multiplication of an `expression` and
/// another `expression`, `complex`, `real`, `rational`, `integer`, `float`,
/// `int` or `bool`.
///
/// `array.mul_`: For in-place multiplication of two `Array`s or slices.
///
/// Notes
/// -----
/// The multiplication is performed with the precision of the output type.
pub inline fn mul_(
    o: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.mul_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.mul_ not implemented for arrays or slices yet.");
    // return array.add(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ " output type"),
        .int => return int.mul_(o, y, .{ .mode = options.mode }),
        .float => return float.mul_(o, y),
        .cfloat => return cfloat.mul_(o, y),
        else => @compileError("zml.mul_ not implemented for " ++ @typeName(O) ++ " output type"),
    }
}

pub inline fn mul_to(
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
        @compileError("zml.mul_to requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.mul_to not implemented for arrays or slices yet.");
    // return array.mul_to(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => switch (types.numericType(C)) {
            .bool => @compileError("zml.mul_to not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = try cast(O, int.mul(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.mul(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.mul(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.mul_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .int => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.mul(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.mul(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.mul(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.mul_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .float => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.mul(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.mul(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.mul(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.mul_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .cfloat => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.mul(x, y, .{ .mode = options.mode }), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.mul(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.mul(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.mul_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.mul_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
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
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.div not implemented for arrays or slices yet.");
    // return array.div(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.div not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.div(x, y),
                .float => return float.div(x, y),
                .cfloat => return cfloat.div(x, y),
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.div(x, y),
                .float => return float.div(x, y),
                .cfloat => return cfloat.div(x, y),
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.div(x, y),
                .cfloat => return cfloat.div(x, y),
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.div(x, y),
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
    }
}

/// Divides two values of any two supported types in-place.
///
/// The function supports division for values of any combination of supported
/// numeric types, `Array`s and slices, and stores the result in the output
/// pointer. The operation performed is
///
/// ```zig
///     o = x / y
/// ```
///
/// Parameters
/// ----------
/// o (pointer to `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The output pointer where the
/// result of the addition operation will be stored, and the left-hand side
/// operand.
///
/// y (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`, `expression`, `Array` or slice): The right-hand side operand.
///
/// options (`struct`): Options for the division operation.
/// - `mode` (`int.Mode`): The mode of the division operation. Only needed when
/// adding two `int` values.
/// - `allocator` (`std.mem.Allocator`): An allocator to use for allocating
/// memory for the output value. Only needed if the output type is of arbitrary
/// precision. May not be used if the output has enouph memory allocated
/// already.
///
/// Returns
/// -------
/// `void`: The result of the division operation is stored in the output
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
/// `@compileError`: If the division operation is not defined for the types of
/// the inputs, if the types of the inputs cannot be coerced to a common type,
/// if the types of the inputs are not supported numeric types, `Array`s or
/// slices, if the output pointer is not a mutable pointer, if the output's
/// child type is not a supported numeric type, `Array` or slice, or if the
/// coerced type is an `Array` and the output's child type is not an `Array`.
///
/// See Also
/// --------
/// `div`: For division of two values and returning the result.
///
/// `int.div_`: For in-place division of an `int` and another `int` or `bool`.
///
/// `float.div_`: For in-place division of a `float` and another `float`, `int`
/// or `bool`.
///
/// `cfloat.div_`: For in-place division of a `cfloat` and another `cfloat`,
/// `float`, `int` or `bool`.
///
/// `integer.div_`: For in-place division of an `integer` and another `integer`,
/// `bool` or `int`.
///
/// `rational.div_`: For in-place division of a `rational` and another
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `real.div_`: For in-place division of a `real` and another `real`,
/// `rational`, `integer`, `float`, `int` or `bool`.
///
/// `complex.div_`: For in-place division of a `complex` and another `complex`,
/// `real`, `rational`, `integer`, `float`, `int` or `bool`.
///
/// `expression.div_`: For in-place division of an `expression` and another
/// `expression`, `complex`, `real`, `rational`, `integer`, `float`, `int` or
/// `bool`.
///
/// `array.div_`: For in-place division of two `Array`s or slices.
///
/// Notes
/// -----
/// The division is performed with the precision of the output type.
pub inline fn div_(
    o: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.div_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.div_ not implemented for arrays or slices yet.");
    // return array.add(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ " output type"),
        .int => return int.div_(o, y, .{ .mode = options.mode }),
        .float => return float.div_(o, y),
        .cfloat => return cfloat.div_(o, y),
        else => @compileError("zml.div_ not implemented for " ++ @typeName(O) ++ " output type"),
    }
}

pub inline fn div_to(
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
        @compileError("zml.div_to requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.div_to not implemented for arrays or slices yet.");
    // return array.div_to(x, y, options.allocator);

    switch (types.numericType(O)) {
        .bool => switch (types.numericType(C)) {
            .bool => @compileError("zml.div_to not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types"),
            .int => o.* = try cast(O, int.div(x, y), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.div(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.div(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.div_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .int => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.div(x, y), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.div(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.div(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.div_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .float => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.div(x, y), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.div(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.div(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.div_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        .cfloat => switch (types.numericType(C)) {
            .bool, .int => o.* = try cast(O, int.div(x, y), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.div(x, y), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.div(x, y), .{ .allocator = options.allocator }),
            else => @compileError("zml.div_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
        },
        else => @compileError("zml.div_to not implemented for " ++ @typeName(O) ++ " output type, and " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " input types yet"),
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

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.zml.eq not implemented for arrays or slices yet.");
    // return array.eq(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.eq not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.eq(x, y),
                .float => return float.eq(x, y),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.eq(x, y),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
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

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.ne not implemented for arrays or slices yet.");
    // return array.ne(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.ne(x, y),
                .float => return float.ne(x, y),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.ne(x, y),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
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

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.lt not implemented for arrays or slices yet.");
    // return array.lt(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.lt not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.lt(x, y),
                .float => return float.lt(x, y),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.lt(x, y),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
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

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.le not implemented for arrays or slices yet.");
    // return array.le(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.le not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.le(x, y),
                .float => return float.le(x, y),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.le(x, y),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
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

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.gt not implemented for arrays or slices yet.");
    // return array.gt(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.gt not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.gt(x, y),
                .float => return float.gt(x, y),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.gt(x, y),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
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

    if (comptime types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.ge not implemented for arrays or slices yet.");
    // return array.ge(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ge not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.ge(x, y),
                .float => return float.ge(x, y),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.ge(x, y),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
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

// Include:
// - sub
// - mul
// - div
// - mod
// - pow
// - neg
// - abs
// - sqrt
// - sin
// - cos
// - tan
// - asin
// - acos
// - atan
// - sinh
// - cosh
// - tanh
// - asinh
// - acosh
// - atanh
// - sec
// - csc
// - cot
// - asec
// - acsc
// - acot
// - sech
// - csch
// - coth
// - asech
// - acsch
// - acoth
// - log
// - log10
// - log2
// - log1p
// - exp
// - any other basic math functions

pub inline fn abs(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !CoerceToArray(@TypeOf(x), Scalar(Scalar(@TypeOf(x)))) {
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
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.abs_(o, .{ .allocator = options.allocator });

    _ = options.allocator;
    switch (types.numericType(O)) {
        .bool => @compileError("zml.abs_ not defined for " ++ @typeName(O) ++ " output type"),
        .int => o.* = int.abs(o.*),
        .float => o.* = float.abs(o.*),
        .cfloat => o.* = @compileError("zml.abs_ not defined for " ++ @typeName(O) ++ " output type"),
        else => @compileError("zml.abs_ not implemented for " ++ @typeName(O) ++ " yet"),
    }
}

pub inline fn abs_to(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    if (comptime !types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_to requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.abs_to(o, x, .{ .allocator = options.allocator });

    _ = options.allocator;
    switch (types.numericType(O)) {
        .bool => switch (types.numericType(X)) {
            .bool => @compileError("zml.abs_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.abs(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.abs(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.abs(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.abs_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .int => switch (types.numericType(X)) {
            .bool => @compileError("zml.abs_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.abs(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.abs(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.abs(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.abs_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .float => switch (types.numericType(X)) {
            .bool => @compileError("zml.abs_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.abs(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.abs(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.abs(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.abs_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .cfloat => switch (types.numericType(X)) {
            .bool => @compileError("zml.abs_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.abs(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.abs(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.abs(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.abs_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        // The idea is that integer, rational, etc. all ap types have special functions
        // that can use a preinitialized element
        else => @compileError("zml.abs_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

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
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ceil_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.ceil_(o, .{ .allocator = options.allocator });

    _ = options.allocator;
    switch (types.numericType(O)) {
        .bool => @compileError("zml.ceil_ not defined for " ++ @typeName(O) ++ " output type"),
        .int => return,
        .float => o.* = float.ceil(o.*),
        .cfloat => @compileError("zml.ceil_ not defined for " ++ @typeName(O) ++ " output type"),
        .complex => @compileError("zml.ceil_ not defined for " ++ @typeName(O) ++ " output type"),
        else => @compileError("zml.ceil_ not implemented for " ++ @typeName(O) ++ " yet"),
    }
}

pub inline fn ceil_to(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ceil_to requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O) or types.isSlice(O))
        return array.ceil_to(o, x, .{ .allocator = options.allocator });

    _ = options.allocator;
    switch (types.numericType(O)) {
        .bool => switch (types.numericType(X)) {
            .bool => @compileError("zml.ceil_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.ceil(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.ceil(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.ceil(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.ceil_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .int => switch (types.numericType(X)) {
            .bool => @compileError("zml.ceil_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.ceil(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.ceil(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.ceil(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.ceil_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .float => switch (types.numericType(X)) {
            .bool => @compileError("zml.ceil_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.ceil(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.ceil(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.ceil(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.ceil_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        .cfloat => switch (types.numericType(X)) {
            .bool => @compileError("zml.ceil_to not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
            .int => o.* = try cast(O, int.ceil(x), .{ .allocator = options.allocator }),
            .float => o.* = try cast(O, float.ceil(x), .{ .allocator = options.allocator }),
            .cfloat => o.* = try cast(O, cfloat.ceil(x), .{ .allocator = options.allocator }),
            else => @compileError("zml.ceil_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
        },
        // The idea is that integer, rational, etc. all ap types have special functions
        // that can use a preinitialized element
        else => @compileError("zml.ceil_to not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " yet"),
    }
}

/// Sets the value of `o` to `x`, using the specified allocator if provided.
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
