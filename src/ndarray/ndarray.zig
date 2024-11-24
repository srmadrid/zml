// MAKE A `squish` FUNCTION THAT GIVEN AN ARRAY "SQUISHES" the dimensions with 1,
// i.e., (1, 1, 10, 1, 1, 9, 1) would become (10, 9). It returns a view.

// MAKE A TRANSFER OWNERSHIP, which given an array and a view of it (with the same
// shape since if the view has smaller shape then it cant happen, unless it is a
// squished shape, as explained in the above paragraph) changes the view to the
// owner and the owner to the view.

// Create a subarray view creator. Useful for i.e. LU decomposition

// FOR MOST FUNCTIONS, MAKE TWO VERSIONS (NAMES NOT FINAL)?? MAYBE JUST OPTIONAL
// PARAMETERS (?T)?
// - add, addInPlace: (self: *Self, left: Self, right: Self) void: adds left +
//   right and stores the result in self
// - addNew, add: (allocator: std.mem.Allocator, left: Self, right: Self) Self:
//   adds left + right and returns a newly allocated array with the result
const std = @import("std");
const zml = @import("../zml.zig");
const core = @import("../core/core.zig");

pub const Iterator = @import("iterators.zig").Iterator;
pub const MultiIterator = @import("iterators.zig").MultiIterator;

const ndarray = @This();

const _add = core.supported._add;
const _sub = core.supported._sub;
const _mul = core.supported._mul;
const _div = core.supported._div;

/// Maximal number of dimensions.
pub const MaxDimensions = 32;

/// An n-dimensional array of elements of type `T`. Initialize with `init` and
/// deinitialize with `deinit`.
///
/// Only supports the following types for `T`:
/// - Integers: u8, u16, u32, u64, u128, i8, i16, i32, i64, i128
/// - Floating-point numbers: f16, f32, f64, f80, f128
/// - C compatibility types
/// - Booleans: bool
/// - Custom types: cf16, cf32, cf64, cf128, BigInt, Fraction, Complex,
///                 Expression
///
/// Stores an allocator for internal memory management.
pub fn NDArray(comptime T: type) type {
    // Catch any attempt to create with unsupported type.
    _ = core.supported.whatSupportedNumericType(T);
    return struct {
        const Self: type = @This();
        /// The data of the array.
        data: []T,
        /// The number dimensions of the array.
        ndim: usize,
        /// The shape of the array, i.e., the dimensions of the array.
        shape: [MaxDimensions]usize,
        /// The strides of the array. These are the number of elements to skip
        /// to get the next element in each dimension.
        strides: [MaxDimensions]usize,
        /// Total number of elements in the array.
        size: usize,
        /// Owner of the data. Set to `null` if itself is the owner.
        base: ?*NDArray(T),
        /// Flags holding info on the storage of the array.
        flags: Flags,
        /// The allocator used for internal memory management.
        allocator: ?std.mem.Allocator,

        /// Initializes an array.
        ///
        /// **Description**:
        ///
        /// Initializes a matrix with the given shape and flags. It is highly
        /// recommended to initialize all elements of the array after calling
        /// this function using the `set` or `setAll` functions. Failure to do
        /// so may result in undefined behavior if the uninitialized elements
        /// are accessed before being set.
        ///
        /// Deinitialize with `deinit`.
        ///
        /// **Input Parameters**:
        /// - `allocator`: allocator to use internally.
        /// - `shape`: shape of the array. Must have length `≤ MaxDimensions`
        /// and no dimension of size `0`.
        /// - `flags`: storage flags for the array (`ownsData` must be true).
        /// Using an empty anonymous struct `.{}` will yield the default flags.
        ///
        /// **Return Values**:
        /// - `Self`: the initialized array.
        /// - `TooManyDimensions`: `shape` is too long (`> MaxDimensions`).
        /// - `ZeroDimension`: at least one dimension is `0`.
        /// - `InvalidFlags`: flags are invalid.
        /// - `OutOfMemory`: `alloc` failed.
        pub fn init(allocator: std.mem.Allocator, shape: []const usize, flags: Flags) !Self {
            if (shape.len > MaxDimensions) {
                return Error.TooManyDimensions;
            }

            if (!flags.ownsData) {
                return Error.InvalidFlags;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return Error.ZeroDimension;
                }
            }

            const ndim = shape.len;
            var size: usize = 1;
            var shapes: [MaxDimensions]usize = [_]usize{0} ** MaxDimensions;
            var strides: [MaxDimensions]usize = [_]usize{0} ** MaxDimensions;
            if (ndim > 0) {
                for (0..ndim) |i| {
                    if (flags.order == .RowMajor) {
                        strides[shape.len - i - 1] = size;
                        size *= shape[shape.len - i - 1];
                    } else {
                        strides[i] = size;
                        size *= shape[i];
                    }
                    shapes[i] = shape[i];
                }
            }

            return Self{
                .data = try allocator.alloc(T, size),
                .ndim = ndim,
                .shape = shapes,
                .strides = strides,
                .size = size,
                .base = null,
                .flags = flags,
                .allocator = allocator,
            };
        }

        /// Initializes a view of an array.
        ///
        /// **Description**:
        ///
        /// Given a preexistsent `NDArray`, creates another array that accesses
        /// the input array's data without owning it. Do not use this array
        /// after freeing the parent. It is not needed to deinit the view.
        ///
        /// **Input Parameters**:
        /// - `parent`: The owner of the data.
        ///
        /// **Return Values**:
        /// - `Self`: the initialized array.
        pub fn view(parent: *const NDArray(T)) !Self {
            const base = if (parent.flags.ownsData) parent else parent.base;
            var flags = parent.flags;
            flags.ownsData = false;
            return Self{
                .data = parent.data,
                .ndim = parent.ndim,
                .shape = parent.shape,
                .strides = parent.strides,
                .size = parent.size,
                .base = base,
                .flags = flags,
                .allocator = null,
            };
        }

        /// Deinitializes an array
        ///
        /// **Description**:
        ///
        /// Deinitializes the array, freeing the data. The user is responsible
        /// for deinitializing the elements of `data` for custom types.
        ///
        /// If the array has `ownsData` set to false, i.e., it is a view, this
        /// only invalidates the data pointer, the shape and the strides. Note
        /// that in this case they need not be deinitialized.
        ///
        /// **Input Parameters**:
        /// - `self`: the array to be deinitialized.
        pub fn deinit(self: *Self) void {
            if (self.flags.ownsData) {
                self.allocator.?.free(self.data);
            }

            self.shape = undefined;
            self.strides = undefined;
            self.data = undefined;
        }

        /// Deinitializes an element of the array.
        ///
        /// **Description**:
        ///
        /// Deinitializes an element of the array. Only defined for custom
        /// types. For a scalar, any position will yield the only element.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to free the element.
        /// - `position`: the location of the element to be freed. Must be a
        /// valid position, unless the array is a scalar.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - `DimensionMismatch`: the length of `position` is not equalm to the
        /// number of dimentions in `self`.
        /// - `PositionOutOfBounds`: the position is out of bounds.
        pub fn deinitElement(self: *Self, position: []const usize) !void {
            const supported = core.supported.whatSupportedNumericType(T);
            switch (supported) {
                .BuiltinInt, .BuiltinFloat, .BuiltinBool, .CustomComplexFloat => @compileError("deinitElement only defined on types that need to be deinitialized."),
                .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
                    if (self.isScalar()) {
                        self.data[0].deinit();
                        return;
                    }

                    try self._checkPosition(position);

                    self.data[self._index(position)].deinit();
                },
                .Unsupported => unreachable,
            }
        }

        /// Deinitializes all elements of the array.
        ///
        /// **Description**:
        ///
        /// Deinitializes all elements of the array. Only defined for custom
        /// types.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to free the elements.
        pub fn deinitAllElements(self: *Self) void {
            const supported = core.supported.whatSupportedNumericType(T);
            switch (supported) {
                .BuiltinInt, .BuiltinFloat, .BuiltinBool, .CustomComplexFloat => @compileError("deinitAllElements only defined on types that need to be deinitialized."),
                .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
                    for (self.data) |element| {
                        element.deinit();
                    }
                },
                .Unsupported => unreachable,
            }
        }

        /// Calculates the index of the element at the given position of the
        /// array.
        ///
        /// No bounds checking is performed.
        inline fn _index(self: Self, position: []const usize) usize {
            var idx: usize = 0;
            for (0..self.ndim) |i| {
                idx += position[i] * self.strides[i];
            }

            return idx;
        }

        /// Checks if the given position is within the bounds of the array and
        /// matches its dimensions.
        ///
        /// This function should be used before accessing elements in the array
        /// to prevent `PositionOutOfBounds` and `DimensionMismatch` errors.
        inline fn _checkPosition(self: Self, position: []const usize) !void {
            if (position.len != self.ndim) {
                return Error.DimensionMismatch;
            }

            for (0..position.len) |i| {
                if (position[i] >= self.shape[i]) {
                    return Error.PositionOutOfBounds;
                }
            }
        }

        /// Sets the value at the given position, if the position is correct.
        ///
        /// **Description**:
        ///
        /// Sets the value at the given position, if the position is correct.
        /// For a scalar, any position will yield the only element.
        ///
        /// For custom types, the array will take ownership of the value and,
        /// therefore, should not be used outside of the array, although it
        /// should still be manually deinitialized before deinitializing the
        /// array. If there was already an element in that positon, it will be
        /// overwritten and a memory leak will occur, unless another reference
        /// of the value was kept elsewhere.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to set the element.
        /// - `position`: the location of the element to be set. Must be a valid
        /// position, unless the array is a scalar.
        /// - `value`: the value to set.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - `DimensionMismatch`: the length of `position` is not equalm to the
        /// number of dimentions in `self`.
        /// - `PositionOutOfBounds`: the position is out of bounds.
        pub fn set(self: *Self, position: []const usize, value: T) !void {
            if (self.isScalar()) {
                self.data[0] = value;
                return;
            }

            try self._checkPosition(position);

            self.data[self._index(position)] = value;
        }

        /// Sets all elements of the array to the input value.
        ///
        /// **Description**:
        ///
        /// Sets all elements of the array to the input value.
        ///
        /// For custom types, the array will take ownership of the value and,
        /// therefore, should not be used outside of the array, although it
        /// should still be manually deinitialized before deinitializing the
        /// array. If there was already an element at some position, it will be
        /// overwritten and a memory leak will occur, unless another copy of the
        /// value was kept elsewhere. The first element of the array is the one
        /// taking ownership of `value`, the rest are newly initialized `T`s
        /// (overwriting any existing element) which also require
        /// deinitializing.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to set the element.
        /// - `value`: the value to set.
        pub fn setAll(self: *Self, value: T) void {
            self.data[0] = value;

            const supported = core.supported.whatSupportedNumericType(T);
            switch (supported) {
                .BuiltinInt, .BuiltinFloat, .BuiltinBool, .CustomComplexFloat => {
                    for (1..self.size) |i| {
                        self.data[i] = value;
                    }
                },
                .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
                    for (1..self.size) |i| {
                        self.data[i] = try T.init(self.allocator, value);
                    }
                },
                .Unsupported => unreachable,
            }
        }

        /// Sets the value at the given position and returns the previous value.
        ///
        /// **Description**:
        ///
        /// Sets the value at the given position and returns the previous value.
        /// For a scalar, any position will yield the only element.
        ///
        /// For custom types, the array takes ownership of the new value.
        /// Therefore, the new value should not be used outside of the array
        /// after being set, although it should still be manually deinitialized
        /// before deinitializing the array. The old value is returned, so it
        /// can be used or deinitialized as needed.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to set the element.
        /// - `position`: the location of the element to be set. Must be a valid
        /// position, unless the array is a scalar.
        /// - `value`: the value to set.
        ///
        /// **Return Values**:
        /// - `T`: the previous value.
        /// - `DimensionMismatch`: the length of `position` is not equalm to the
        /// number of dimentions in `self`.
        /// - `PositionOutOfBounds`: the position is out of bounds.
        pub fn replace(self: *Self, position: []const usize, value: T) !T {
            var idx: usize = undefined;
            if (self.isScalar()) {
                idx = 0;
            } else {
                try self._checkPosition(position);
                idx = self._index(position);
            }

            const prev: T = self.data[idx];
            self.data[idx] = value;
            return prev;
        }

        /// Sets the value at the given position by updating the already
        /// prexistent element.
        ///
        /// **Description**:
        ///
        /// Sets the value at the given position by updating the already
        /// prexistent element. For a scalar, any position will yield the only
        /// element.
        ///
        /// This function is useful when you want to change the value at a
        /// position in the array without causing the array to take ownership of
        /// the new value or overwrite a prexistent element. The new value can
        /// continue to be used outside of the array after this function is
        /// called.
        ///
        /// For builting numeric types (uxx, ixx, fxx) it is equivalent to
        /// `set`.
        ///
        /// For custom types, calling this on an uninitiallized element will
        /// cause undefined behavior.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to set the element.
        /// - `position`: the location of the element to be set. Must be a valid
        /// position, unless the array is a scalar.
        /// - `value`: the value to set.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - `DimensionMismatch`: the length of `position` is not equalm to the
        /// number of dimentions in `self`.
        /// - `PositionOutOfBounds`: the position is out of bounds.
        pub fn update(self: *Self, position: []const usize, value: anytype) !void {
            var idx: usize = undefined;
            if (self.isScalar()) {
                idx = 0;
            } else {
                try self._checkPosition(position);
                idx = self._index(position);
            }
            const supported = core.supported.whatSupportedNumericType(T);
            switch (supported) {
                .BuiltinInt, .BuiltinFloat, .BuiltinBool, .CustomComplexFloat => {
                    self.data[idx] = value;
                },
                .CustomReal, .CustomComplex, .CustomExpression => {
                    try self.data[idx].set(value);
                },
                else => unreachable,
            }
        }

        /// Sets the values of all the elements of the array by updating the
        /// already prexistent element.
        ///
        /// **Description**:
        ///
        /// Sets the value at the given position by updating the already
        /// prexistent element. For a scalar, any position will yield the only
        /// element.
        ///
        /// This function is useful when you want to change all the elements in
        /// in the array without causing the array to take ownership of the new
        /// value or overwrite elements. The new value can continue to be used
        /// outside of the array after this function is called.
        ///
        /// For builting numeric types (uxx, ixx, fxx) it is equivalent to
        /// `setAll`.
        ///
        /// For custom types, calling this on an uninitiallized element will
        /// cause undefined behavior.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to set the element.
        /// - `value`: the value to set.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - ...
        pub fn updateAll(self: *Self, value: T) !void {
            const supported = core.supported.whatSupportedNumericType(T);
            switch (supported) {
                .BuiltinInt, .BuiltinFloat, .BuiltinBool, .CustomComplexFloat => {
                    for (0..self.size) |i| {
                        self.data[i] = value;
                    }
                },
                .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
                    for (0..self.size) |i| {
                        try self.data[i].set(value);
                    }
                },
                .Unsupported => unreachable,
            }
        }

        /// Gets the element at the given position, if the position is correct.
        ///
        /// **Description**:
        ///
        /// Gets the element at the given position, if the position is correct.
        /// For a scalar, any position will yield the only element.
        ///
        /// For custom types, returns a shallow copy of the element, so if it is
        /// freed from inside or outside of the array, the other reference
        /// becomes unusable.
        ///
        /// **Input Parameters**:
        /// - `self`: the array in which to set the element.
        /// - `position`: the location of the element to be set. Must be a valid
        /// position, unless the array is a scalar.
        ///
        /// **Return Values**:
        /// - `T`: the element at `position`.
        /// - `DimensionMismatch`: the length of `position` is not equalm to the
        /// number of dimentions in `self`.
        /// - `PositionOutOfBounds`: the position is out of bounds.
        pub fn get(self: Self, position: []const usize) !T {
            if (self.isScalar()) {
                return self.data[0];
            }

            try self._checkPosition(position);
            return self.data[self._index(position)];
        }

        /// Checks if the array is a scalar.
        ///
        /// **Description**:
        ///
        /// Checks if the array is a scalar, i.e., has only one element.
        ///
        /// **Input Parameters**:
        /// - `self`: the array to check.
        ///
        /// **Return Values**:
        /// - `bool`: whether the array is a scalar or not.
        pub fn isScalar(self: Self) bool {
            return self.size == 1;
        }

        /// Checks if the array is a vector.
        ///
        /// **Description**:
        ///
        /// Checks if the array is a vector, i.e., has shape `{n}`.
        ///
        /// **Input Parameters**:
        /// - `self`: the array to check.
        ///
        /// **Return Values**:
        /// - `bool`: whether the array is a vector or not.
        pub fn isVector(self: Self) bool {
            return self.shape.len == 1;
        }

        /// Checks if the array is a row vector.
        ///
        /// **Description**:
        ///
        /// Checks if the array is a vector, i.e., has shape `{ 1, n }`.
        ///
        /// **Input Parameters**:
        /// - `self`: the array to check.
        ///
        /// **Return Values**:
        /// - `bool`: whether the array is a row vector or not.
        pub fn isRowVector(self: Self) bool {
            return self.shape.len == 2 and self.shape[0] == 1 and self.shape[1] > 1;
        }

        /// Checks if the array is a column vector.
        ///
        /// **Description**:
        ///
        /// Checks if the array is a vector, i.e., has shape `{ n, 1 }`.
        ///
        /// **Input Parameters**:
        /// - `self`: the array to check.
        ///
        /// **Return Values**:
        /// - `bool`: whether the array is a column vector or not.
        pub fn isColVector(self: Self) bool {
            return self.shape.len == 2 and self.shape[0] > 1 and self.shape[1] == 1;
        }

        /// Checks if the array is a matrix.
        ///
        /// **Description**:
        ///
        /// Checks if the array is a matrix, i.e., has shape `{ m, n }`.
        ///
        /// **Input Parameters**:
        /// - `self`: the array to check.
        ///
        /// **Return Values**:
        /// - `bool`: whether the array is a matrix or not.
        pub fn isMatrix(self: Self) bool {
            return self.shape.len == 2;
        }

        // TODO: Implement the following functions (only for matrices)
        // - isSquare
        // - isUpperTriangular
        // - isLowerTriangular
        // - isSymmetric
        // - isSkewSymmetric
        // - isHermitian
        // - isSkewHermitian
        // - more property checking functions

        /// Computes elementwise addition (self = left + right).
        ///
        /// **Description**:
        ///
        /// Computes elementwise addition (self = left + right) applying
        /// broadcasting rules. Arrays must be broadcastable.
        ///
        /// **Input Parameters**:
        /// - `self`: the output array. Must have the same shape as the full
        /// broadcast.
        /// - `left`: the array on the left side.
        /// - `right`: the array on the right side.
        ///
        /// **Output Parameters**:
        /// - `self`: Overwritten by the result of the operation.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - `NotBroadcastable`: the shapes of the arrays cannot be
        /// broadcasted.
        /// - `IncompatibleDimensions`: the `out` array does not have the shape
        /// of the full broadcast.
        pub fn add(self: *Self, left: Self, right: Self) !void {
            var iter: MultiIterator(T) = try MultiIterator(T).init(&.{ self.*, left, right }, self.flags);
            if (!std.mem.eql(usize, self.shape[0..self.ndim], iter.shape[0..iter.ndim])) {
                return Error.IncompatibleDimensions;
            }

            _add(&self.data[0], left.data[0], right.data[0]);
            while (iter.next() != null) {
                _add(&self.data[iter.iterators[0].index], left.data[iter.iterators[1].index], right.data[iter.iterators[2].index]);
            }
        }

        /// Computes elementwise subtraction (self = left - right).
        ///
        /// **Description**:
        ///
        /// Computes elementwise subtraction (self = left - right) applying
        /// broadcasting rules. Arrays must be broadcastable.
        ///
        /// **Input Parameters**:
        /// - `self`: the output array. Must have the same shape as the full
        /// broadcast.
        /// - `left`: the array on the left side.
        /// - `right`: the array on the right side.
        ///
        /// **Output Parameters**:
        /// - `self`: Overwritten by the result of the operation.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - `NotBroadcastable`: the shapes of the arrays cannot be
        /// broadcasted.
        /// - `IncompatibleDimensions`: the `out` array does not have the shape
        /// of the full broadcast.
        pub fn sub(self: *Self, left: Self, right: Self) !void {
            var iter: MultiIterator(T) = try MultiIterator(T).init(&.{ self.*, left, right }, self.flags);
            if (!std.mem.eql(usize, self.shape[0..self.ndim], iter.shape[0..iter.ndim])) {
                return Error.IncompatibleDimensions;
            }

            _sub(&self.data[0], left.data[0], right.data[0]);
            while (iter.next() != null) {
                _sub(&self.data[iter.iterators[0].index], left.data[iter.iterators[1].index], right.data[iter.iterators[2].index]);
            }
        }

        /// Computes elementwise multiplication (self = left .* right).
        ///
        /// **Description**:
        ///
        /// Computes elementwise multiplication (self = left .* right) applying
        /// broadcasting rules. Arrays must be broadcastable.
        ///
        /// **Input Parameters**:
        /// - `self`: the output array. Must have the same shape as the full
        /// broadcast.
        /// - `left`: the array on the left side.
        /// - `right`: the array on the right side.
        ///
        /// **Output Parameters**:
        /// - `self`: Overwritten by the result of the operation.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - `NotBroadcastable`: the shapes of the arrays cannot be
        /// broadcasted.
        /// - `IncompatibleDimensions`: the `out` array does not have the shape
        /// of the full broadcast.
        pub fn mul(self: *Self, left: Self, right: Self) !void {
            var iter: MultiIterator(T) = try MultiIterator(T).init(&.{ self.*, left, right }, self.flags);
            if (!std.mem.eql(usize, self.shape[0..self.ndim], iter.shape[0..iter.ndim])) {
                return Error.IncompatibleDimensions;
            }

            _mul(&self.data[0], left.data[0], right.data[0]);
            while (iter.next() != null) {
                _mul(&self.data[iter.iterators[0].index], left.data[iter.iterators[1].index], right.data[iter.iterators[2].index]);
            }
        }

        /// Computes elementwise division (self = left ./ right).
        ///
        /// **Description**:
        ///
        /// Computes elementwise multiplication (self = left ./ right) applying
        /// broadcasting rules. Arrays must be broadcastable.
        ///
        /// **Input Parameters**:
        /// - `self`: the output array. Must have the same shape as the full
        /// broadcast.
        /// - `left`: the array on the left side.
        /// - `right`: the array on the right side.
        ///
        /// **Output Parameters**:
        /// - `self`: Overwritten by the result of the operation.
        ///
        /// **Return Values**:
        /// - `void`: the execution was successful.
        /// - `NotBroadcastable`: the shapes of the arrays cannot be
        /// broadcasted.
        /// - `IncompatibleDimensions`: the `out` array does not have the shape
        /// of the full broadcast.
        pub fn div(self: *Self, left: Self, right: Self) !void {
            var iter: MultiIterator(T) = try MultiIterator(T).init(&.{ self.*, left, right }, self.flags);
            if (!std.mem.eql(usize, self.shape[0..self.ndim], iter.shape[0..iter.ndim])) {
                return Error.IncompatibleDimensions;
            }

            _div(&self.data[0], left.data[0], right.data[0]);
            while (iter.next() != null) {
                _div(&self.data[iter.iterators[0].index], left.data[iter.iterators[1].index], right.data[iter.iterators[2].index]);
            }
        }

        /// Returns a view of the input array with the axes transposed.
        ///
        /// **Description**:
        ///
        /// Returns a view of the input array with the axes transposed according
        /// to the give axes, or reversed otherwise.
        ///
        /// **Input Parameters**:
        /// - `self`: the input array.
        /// - `axes`: optional axes. If provided must be a permutation of
        /// `[0, 1, ..., ndim - 1]`.
        ///
        /// **Return Values**:
        /// - `NDArray(T)`: the execution was successful.
        /// - `InvalidAxes`: the `axes` array was wrong.
        pub fn transpose(self: *const Self, axes: ?[]const usize) !NDArray(T) {
            var axess = [_]usize{0} ** MaxDimensions;
            if (axes == null) {
                for (0..self.ndim) |i| {
                    axess[i] = self.ndim - i - 1;
                }
            } else {
                for (0..self.ndim) |i| {
                    axess[i] = axes.?[i];
                }
            }

            var flags = self.flags;
            flags.ownsData = false;

            var result = Self{
                .data = self.data,
                .ndim = self.ndim,
                .shape = .{0} ** MaxDimensions,
                .strides = .{0} ** MaxDimensions,
                .size = self.size,
                .base = @constCast(self),
                .flags = flags,
                .allocator = null,
            };

            for (0..result.ndim) |i| {
                result.shape[i] = self.shape[axess[i]];
                result.strides[i] = self.strides[axess[i]];
            }

            return result;
        }

        /// Namespace for BLAS functions.
        pub const BLAS = struct {
            /// Computes the sum of magnitudes of the vector elements.
            ///
            /// **Description**:
            ///
            /// The routine computes the sum of the magnitudes of elements of a
            /// real array, or the sum of magnitudes of the real and imaginary
            /// parts of elements of a complex array:
            /// ```zig
            /// res = |x[1].Re| + |x[1].Im| + |x[2].Re| + |x[2].Im| + ... + |x[n].Re| + |x[n].Im|,
            /// ```
            /// where `x` is an `NDArray`.
            ///
            /// **Input Parameters**:
            /// - `allocator`: an optional `std.mem.Allocator`. Only needed when
            /// `T` is `BigInt`, `Fraction`, `Complex` or `Expression`.
            /// - `x`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `T`: The sum of magnitudes of real and imaginary parts of
            /// all elements of the vector.
            /// - `someError`: blabla
            pub fn asum(allocator: ?std.mem.Allocator, x: NDArray(T)) !T {
                return @import("BLAS/BLAS.zig").asum(allocator, T, x);
            }
        };

        /// Namespace for LAPACK functions.
        pub const LAPACK = struct {
            /// Computes the LU factorization of a general m-by-n matrix.
            ///
            /// **Description**:
            ///
            /// The routine computes the LU factorization of a general
            /// `m`-by-`n` matrix `A` as `A = P*L*U`, where `P` is a permutation
            /// matrix, `L` is lower triangular with unit diagonal elements
            /// (lower trapezoidal if `m > n`) and `U` is upper triangular
            /// (upper trapezoidal if `m < n`). The routine uses partial
            /// pivoting, with row interchanges.
            ///
            /// **Input Parameters**:
            /// - `a`: `NDArray` of shape `{m, n}`. Contains the matrix `A`.
            ///
            /// **Output Parameters**:
            /// - `a`: Overwritten by `L` and `U`. The unit diagonal elements of
            /// `L` are not stored.
            /// - `ipiv`: `NDArray`, of shape `{min(m, n)}`. Contains the pivot
            /// indices; for `1 ≤ i ≤ min(m, n)`, row `i` was interchanged with
            /// row `ipiv[i]`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// ...
            //pub fn getrf(a: *Self, ipiv: *NDArray(usize)) !void {
            //    return @import("ndarray/LAPACK/getrf.zig").getrf(T, a, ipiv);
            //}
            pub fn tmp() usize {
                return 3;
            }
        };
    };
}

/// Errors that can occur when working with an `NDArray`.
pub const Error = error{
    /// Invalid flags.
    InvalidFlags,
    /// Too many dimensions.
    TooManyDimensions,
    /// The dimensions of the array and the position do not match.
    DimensionMismatch,
    /// Some position is out of bounds.
    PositionOutOfBounds,
    /// A dimension is zero.
    ZeroDimension,
    /// The array is not a scalar.
    NotScalar,
    /// The array is not a vector.
    NotVector,
    /// The array is not a matrix.
    NotMatrix,
    /// The dimensions of the array are not compatible.
    IncompatibleDimensions,
    /// No allocator was provided when needed.
    NoAllocator,
};

/// Flags representing information on the storage of an array.
pub const Flags = packed struct {
    /// Row major element storage (right to left).
    order: Order = .RowMajor,
    /// The array owns the data, and it will be freed when the array is
    /// deinitialized. If it does not own it, it is assumed to be a view.
    ownsData: bool = true,
    /// The data of the array can or not be modified.
    writeable: bool = true,
};

/// Order of the elements in the array.
pub const Order = enum(u1) {
    /// Row major element storage (right to left).
    RowMajor,
    /// Column major element storage (left to right).
    ColumnMajor,
};

test "init" {
    const a: std.mem.Allocator = std.testing.allocator;

    const tooBig: [MaxDimensions + 1]usize = [_]usize{1} ** (MaxDimensions + 1);
    try std.testing.expectError(Error.TooManyDimensions, NDArray(f64).init(a, &tooBig, .{}));
    try std.testing.expectError(Error.ZeroDimension, NDArray(f64).init(a, &[_]usize{ 2, 3, 5, 6, 0, 1, 8 }, .{}));
    try std.testing.expectError(Error.InvalidFlags, NDArray(f64).init(a, &[_]usize{ 2, 2 }, .{ .order = .ColumnMajor }));
    try std.testing.expectError(Error.InvalidFlags, NDArray(f64).init(a, &[_]usize{ 2, 2 }, .{ .order = .RowMajor }));
    try std.testing.expectError(Error.InvalidFlags, NDArray(f64).init(a, &[_]usize{ 2, 2 }, .{ .ownsData = false }));

    var scalar1: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{}, .{});
    defer scalar1.deinit();
    try std.testing.expect(scalar1.data.len == 1);
    try std.testing.expect(std.mem.eql(usize, &[_]usize{}, scalar1.shape));
    try std.testing.expect(std.mem.eql(usize, &[_]usize{}, scalar1.strides));
    try std.testing.expect(scalar1.size == 1);
    try std.testing.expect(scalar1.flags.order == .RowMajor);
    try std.testing.expect(scalar1.flags.ownsData);
    try std.testing.expect(scalar1.flags.writeable);

    var scalar2: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{}, .{});
    defer scalar2.deinit();
    try std.testing.expect(scalar2.data.len == 1);
    try std.testing.expect(std.mem.eql(usize, &[_]usize{}, scalar2.shape));
    try std.testing.expect(std.mem.eql(usize, &[_]usize{}, scalar2.strides));
    try std.testing.expect(scalar2.size == 1);
    try std.testing.expect(scalar2.flags.order == .RowMajor);
    try std.testing.expect(scalar2.flags.ownsData);
    try std.testing.expect(scalar2.flags.writeable);

    var A: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 10, 5, 8, 3 }, .{});
    defer A.deinit();

    try std.testing.expect(A.data.len == 1200);
    try std.testing.expect(std.mem.eql(usize, &[_]usize{ 10, 5, 8, 3 }, A.shape));
    try std.testing.expect(std.mem.eql(usize, &[_]usize{ 120, 24, 3, 1 }, A.strides));
    try std.testing.expect(A.size == 1200);
    try std.testing.expect(A.flags.order == .RowMajor);
    try std.testing.expect(A.flags.ownsData);
    try std.testing.expect(A.flags.writeable);

    var B: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 10, 5, 8, 3 }, .{ .order = .ColumnMajor });
    defer B.deinit();

    try std.testing.expect(B.data.len == 1200);
    try std.testing.expect(std.mem.eql(usize, &[_]usize{ 10, 5, 8, 3 }, B.shape));
    try std.testing.expect(std.mem.eql(usize, &[_]usize{ 1, 10, 50, 400 }, B.strides));
    try std.testing.expect(B.size == 1200);
    try std.testing.expect(B.flags.order == .ColumnMajor);
    try std.testing.expect(B.flags.ownsData);
    try std.testing.expect(B.flags.writeable);

    var scalar: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{}, .{});
    defer scalar.deinit();
    try std.testing.expect(scalar.data.len == 1);
    try std.testing.expect(std.mem.eql(usize, &[_]usize{}, scalar.shape));
    try std.testing.expect(std.mem.eql(usize, &[_]usize{}, scalar.strides));
    try std.testing.expect(scalar.size == 1);
    try std.testing.expect(scalar.flags.order == .RowMajor);
    try std.testing.expect(scalar.flags.ownsData);
    try std.testing.expect(scalar.flags.writeable);

    // MORE FOR VIEWS WHEN THEY ARE IMPLEMENTED
}

test "_index" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 10, 5, 8, 3 }, .{});
    defer A.deinit();

    try std.testing.expect(A._index(&[_]usize{ 3, 1, 0, 2 }) == 386);

    var B: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 10, 5, 8, 3 }, .{ .order = .ColumnMajor });
    defer B.deinit();

    try std.testing.expect(B._index(&[_]usize{ 3, 1, 0, 2 }) == 813);
}

test "_checkPosition" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 10, 5, 8, 3 }, .{});
    defer A.deinit();

    try std.testing.expectError(Error.DimensionMismatch, A._checkPosition(&[_]usize{3}));
    try std.testing.expectError(Error.PositionOutOfBounds, A._checkPosition(&[_]usize{ 3, 1, 9000, 2 }));
}

test "set" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 3, 2, 4 }, .{});
    defer A.deinit();
    // This is one way of iterating through an array.
    var elem: f64 = 1;
    for (0..A.shape[0]) |i| {
        for (0..A.shape[1]) |j| {
            for (0..A.shape[2]) |k| {
                try A.set(&[_]usize{ i, j, k }, elem);
                elem += 1;
            }
        }
    }
    try std.testing.expect(std.mem.eql(f64, A.data, &[_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 }));

    var B: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 3, 2, 4 }, .{ .order = .ColumnMajor });
    defer B.deinit();
    elem = 1;
    for (0..B.shape[0]) |i| {
        for (0..B.shape[1]) |j| {
            for (0..B.shape[2]) |k| {
                try B.set(&[_]usize{ i, j, k }, elem);
                elem += 1;
            }
        }
    }
    try std.testing.expect(std.mem.eql(f64, B.data, &[_]f64{ 1, 9, 17, 5, 13, 21, 2, 10, 18, 6, 14, 22, 3, 11, 19, 7, 15, 23, 4, 12, 20, 8, 16, 24 }));
}

test "setAll" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 5, 10, 4 }, .{});
    defer A.deinit();
    // This is one way of iterating through an array.
    A.setAll(5);
    try std.testing.expect(std.mem.eql(f64, A.data, &[_]f64{5} ** 200));
}

test "replace" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &[_]usize{ 3, 2, 4 }, .{});
    defer A.deinit();
    // This is one way of iterating through an array.
    var elem: f64 = 1;
    for (0..A.shape[0]) |i| {
        for (0..A.shape[1]) |j| {
            for (0..A.shape[2]) |k| {
                try A.set(&[_]usize{ i, j, k }, elem);
                elem += 1;
            }
        }
    }
    const replaced: f64 = try A.replace(&[_]usize{ 0, 1, 1 }, 20);
    try std.testing.expect(replaced == 6);
    try std.testing.expect(try A.get(&[_]usize{ 0, 1, 1 }) == 20);
}

test {
    _ = @import("BLAS/BLAS.zig");
    std.testing.refAllDeclsRecursive(@This());
}
