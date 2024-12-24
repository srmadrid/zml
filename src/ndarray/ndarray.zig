// MAKE A `squish` FUNCTION THAT GIVEN AN ARRAY "SQUISHES" the dimensions with 1,
// i.e., (1, 1, 10, 1, 1, 9, 1) would become (10, 9). It returns a view.

// MAKE A TRANSFER OWNERSHIP, which given an array and a view of it (with the same
// shape since if the view has smaller shape then it cant happen, unless it is a
// squished shape, as explained in the above paragraph) changes the view to the
// owner and the owner to the view.

// Create a subarray view creator. Useful for i.e. LU decomposition

// FOR MOST FUNCTIONS, MAKE TWO VERSIONS (NAMES NOT FINAL)?? MAYBE JUST OPTIONAL
// PARAMETERS (?T)?
// - add, addInPlace: (self: *NDArray(T), left: NDArray(T), right: NDArray(T)) void: adds left +
//   right and stores the result in self
// - addNew, add: (allocator: std.mem.Allocator, left: NDArray(T), right: NDArray(T)) NDArray(T):
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
const scalar = core.supported.scalar;

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
        /// The data of the array.
        data: []T,
        /// The number dimensions of the array.
        ndim: usize,
        /// The shape of the array, i.e., the dimensions of the array.
        shape: [MaxDimensions]usize,
        /// The strides of the array. These are the number of elements to skip
        /// to get the next element in each dimension.
        strides: [MaxDimensions]isize,
        /// Offset of the first element of the array.
        offset: usize,
        /// Total number of elements in the array.
        size: usize,
        /// Owner of the data. Set to `null` if it is the owner.
        base: ?*const NDArray(T),
        /// Flags holding info on the storage of the array.
        flags: Flags,
        /// The allocator used for internal memory management. Set to `null` if
        /// the array is a view.
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
        /// - `NDArray(T)`: the initialized array.
        /// - `TooManyDimensions`: `shape` is too long (`> MaxDimensions`).
        /// - `ZeroDimension`: at least one dimension is `0`.
        /// - `InvalidFlags`: flags are invalid.
        /// - `OutOfMemory`: `alloc` failed.
        pub fn init(allocator: std.mem.Allocator, shape: []const usize, flags: Flags) !NDArray(T) {
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
            var shapes: [MaxDimensions]usize = .{0} ** MaxDimensions;
            var strides: [MaxDimensions]isize = .{0} ** MaxDimensions;
            if (ndim > 0) {
                for (0..ndim) |i| {
                    const idx = if (flags.order == .RowMajor) shape.len - i - 1 else i;

                    strides[idx] = @intCast(size);
                    size *= shape[idx];

                    shapes[i] = shape[i];
                }
            }

            return NDArray(T){
                .data = try allocator.alloc(T, size),
                .ndim = ndim,
                .shape = shapes,
                .strides = strides,
                .offset = 0,
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
        /// - `NDArray(T)`: the initialized array.
        pub fn view(parent: *const NDArray(T)) !NDArray(T) {
            const base = if (parent.base != null) parent.base else parent;
            var flags = parent.flags;
            flags.ownsData = false;
            return NDArray(T){
                .data = parent.data,
                .ndim = parent.ndim,
                .shape = parent.shape,
                .strides = parent.strides,
                .offset = parent.offset,
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
        pub fn deinit(self: *NDArray(T)) void {
            if (self.flags.ownsData) {
                self.allocator.?.free(self.data);
            }

            self.shape = .{0} ** MaxDimensions;
            self.strides = .{0} ** MaxDimensions;
            self.data = undefined;
        }

        /// Calculates the index of the element at the given position of the
        /// array.
        ///
        /// No bounds checking is performed.
        inline fn _index(self: NDArray(T), position: []const usize) usize {
            var idx: usize = self.offset;
            for (0..self.ndim) |i| {
                const stride = self.strides[i];
                if (stride < 0) {
                    idx -= position[i] * @abs(stride);
                } else {
                    idx += position[i] * @as(usize, @intCast(self.strides[i]));
                }
            }

            return idx;
        }

        /// Checks if the given position is within the bounds of the array and
        /// matches its dimensions.
        ///
        /// This function should be used before accessing elements in the array
        /// to prevent `PositionOutOfBounds` and `DimensionMismatch` errors.
        inline fn _checkPosition(self: NDArray(T), position: []const usize) !void {
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
        pub fn set(self: *NDArray(T), position: []const usize, value: T) !void {
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
        pub fn setAll(self: *NDArray(T), value: T) void {
            self.data[0] = value;

            const supported = core.supported.whatSupportedNumericType(T);
            switch (supported) {
                .BuiltinInt, .BuiltinFloat, .BuiltinBool, .Complex => {
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
        pub fn replace(self: *NDArray(T), position: []const usize, value: T) !T {
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
        /// - `DimensionMismatch`: the length of `position` is not equal to the
        /// number of dimentions in `self`.
        /// - `PositionOutOfBounds`: the position is out of bounds.
        pub fn update(self: *NDArray(T), position: []const usize, value: anytype) !void {
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
        pub fn updateAll(self: *NDArray(T), value: T) !void {
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
        pub fn get(self: NDArray(T), position: []const usize) !T {
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
        pub fn isScalar(self: NDArray(T)) bool {
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
        pub fn isVector(self: NDArray(T)) bool {
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
        pub fn isRowVector(self: NDArray(T)) bool {
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
        pub fn isColVector(self: NDArray(T)) bool {
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
        pub fn isMatrix(self: NDArray(T)) bool {
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
        pub fn add(self: *NDArray(T), left: NDArray(T), right: NDArray(T)) !void {
            // Fast loop if possible.
            if (self.flags.order == left.flags.order and self.flags.order == right.flags.order and
                self.ndim == left.ndim and self.ndim == right.ndim and
                std.mem.eql(usize, self.shape[0..self.ndim], left.shape[0..left.ndim]) and
                std.mem.eql(usize, self.shape[0..self.ndim], right.shape[0..right.ndim]))
            {
                for (0..self.size) |i| {
                    _add(&self.data[i], left.data[i], right.data[i]);
                }

                return;
            }

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
        pub fn sub(self: *NDArray(T), left: NDArray(T), right: NDArray(T)) !void {
            // Fast loop if possible.
            if (self.flags.order == left.flags.order and self.flags.order == right.flags.order and
                self.ndim == left.ndim and self.ndim == right.ndim and
                std.mem.eql(usize, self.shape[0..self.ndim], left.shape[0..left.ndim]) and
                std.mem.eql(usize, self.shape[0..self.ndim], right.shape[0..right.ndim]))
            {
                for (0..self.size) |i| {
                    _sub(&self.data[i], left.data[i], right.data[i]);
                }

                return;
            }

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
        pub fn mul(self: *NDArray(T), left: NDArray(T), right: NDArray(T)) !void {
            // Fast loop if possible.
            if (self.flags.order == left.flags.order and self.flags.order == right.flags.order and
                self.ndim == left.ndim and self.ndim == right.ndim and
                std.mem.eql(usize, self.shape[0..self.ndim], left.shape[0..left.ndim]) and
                std.mem.eql(usize, self.shape[0..self.ndim], right.shape[0..right.ndim]))
            {
                for (0..self.size) |i| {
                    _mul(&self.data[i], left.data[i], right.data[i]);
                }

                return;
            }

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
        pub fn div(self: *NDArray(T), left: NDArray(T), right: NDArray(T)) !void {
            // Fast loop if possible.
            if (self.flags.order == left.flags.order and self.flags.order == right.flags.order and
                self.ndim == left.ndim and self.ndim == right.ndim and
                std.mem.eql(usize, self.shape[0..self.ndim], left.shape[0..left.ndim]) and
                std.mem.eql(usize, self.shape[0..self.ndim], right.shape[0..right.ndim]))
            {
                for (0..self.size) |i| {
                    _div(&self.data[i], left.data[i], right.data[i]);
                }

                return;
            }

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
        pub fn transpose(self: *const NDArray(T), axes: ?[]const usize) !NDArray(T) {
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
            flags.order = if (self.flags.order == .RowMajor) .ColumnMajor else .RowMajor;

            var result = NDArray(T){
                .data = self.data,
                .ndim = self.ndim,
                .shape = .{0} ** MaxDimensions,
                .strides = .{0} ** MaxDimensions,
                .offset = self.offset,
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

        /// Returns a view of the input array as a 1D array.
        ///
        /// **Description**:
        ///
        /// Returns a view of the input array as a 1D array. The data is not
        /// copied, so changes in the returned array will affect the original
        /// array.
        ///
        /// **Input Parameters**:
        /// - `self`: the input array.
        ///
        /// **Return Values**:
        /// - `NDArray(T)`: the flattened array.
        pub fn flatten(self: *const NDArray(T)) NDArray(T) {
            // Once order is expanded, only work on contiguous arrays.
            var flags = self.flags;
            flags.ownsData = false;

            return NDArray(T){
                .data = self.data,
                .ndim = 1,
                .shape = .{self.size} ++ .{0} ** (MaxDimensions - 1),
                .strides = .{1} ++ .{0} ** (MaxDimensions - 1),
                .offset = self.offset,
                .size = self.size,
                .base = if (self.flags.ownsData) self else self.base,
                .flags = flags,
                .allocator = null,
            };
        }

        pub fn slice(self: *const NDArray(T), slices: []const Slice) !NDArray(T) {
            if (slices.len == 0 or slices.len > self.ndim) {
                return error.DimensionMismatch;
            }

            var ndim = self.ndim;
            var size: usize = 1;
            var shapes: [MaxDimensions]usize = .{0} ** MaxDimensions;
            var strides: [MaxDimensions]isize = .{0} ** MaxDimensions;
            var offset: usize = self.offset;

            var i: usize = 0;
            var j: usize = 0;
            while (i < self.ndim) {
                if (i >= slices.len) {
                    shapes[j] = self.shape[i];
                    strides[j] = self.strides[i];
                    size *= self.shape[i];
                    j += 1;
                    i += 1;
                    continue;
                }

                var sl = slices[i];
                if (sl.step > 0) {
                    if (sl.start >= self.shape[i]) {
                        sl.start = self.shape[i] - 1;
                        sl.stop = self.shape[i];
                    } else if (sl.stop > self.shape[i]) {
                        sl.stop = self.shape[i];
                    }
                } else {
                    if (sl.stop >= self.shape[i]) {
                        sl.stop = self.shape[i] - 1;
                        sl.start = self.shape[i];
                    } else if (sl.start > self.shape[i]) {
                        sl.start = self.shape[i];
                    }
                }

                const len = sl.len();
                if (len == 0 or len == 1) {
                    offset += sl.start * @as(usize, @intCast(self.strides[i]));
                    ndim -= 1;
                } else {
                    shapes[j] = len;
                    strides[j] = self.strides[i] * slices[i].step;
                    size *= len;
                    offset += sl.start * @as(usize, @intCast(self.strides[i]));
                    j += 1;
                }

                i += 1;
            }

            // Adjust the flags for the new array
            var flags = self.flags;
            flags.ownsData = false; // The new view does not own the data

            // Return the new sliced array
            return NDArray(T){
                .data = self.data, // Point to the original data
                .ndim = ndim, // Updated number of dimensions
                .shape = shapes, // Updated shapes
                .strides = strides, // Updated strides
                .offset = offset, // Updated offset
                .size = size, // Updated total size
                .base = if (self.flags.ownsData) self else self.base, // Keep the original base if needed
                .flags = flags, // Updated flags
                .allocator = null, // No allocator since this is a view
            };
        }

        /// Namespace for BLAS routines for simple types that do not require
        /// memory management.
        pub const BLAS = struct {
            /// Computes the sum of magnitudes of the vector elements.
            ///
            /// **Description**:
            ///
            /// The `asum` routine computes the sum of the magnitudes of
            /// elements of a 1D real array, or the sum of magnitudes of the
            /// real and imaginary parts of the elements of a complex 1D array:
            /// ```zig
            /// res = |x[1].re| + |x[1].im| + |x[2].re| + |x[2].im| + ... + |x[n].re| + |x[n].im|,
            /// ```
            /// where `x` is an `NDArray`.
            ///
            /// **Input Parameters**:
            /// - `x`: 1D `NDArray`.
            ///
            /// **Return Values**:
            /// - `T`: The sum of magnitudes of real and imaginary parts of
            /// all elements of the vector.
            /// - `IncompatibleDimensions`: the array is not 1D.
            pub fn asum(x: NDArray(T)) !scalar(T) {
                if (x.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").asum(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]));
            }

            /// Computes a vector-scalar product and adds the result to a
            /// vector.
            ///
            /// **Description**:
            ///
            /// The `axpy` routine performs a vector-vector operation defined
            /// as:
            /// ```zig
            /// y = a*x + y,
            /// ```
            /// where `a` is a scalar and `x` and `y` are 1D arrays with the
            /// same size.
            ///
            /// **Input Parameters**:
            /// - `a`: scalar.
            /// - `x`: `NDArray` of shape `{n}`.
            /// - `y`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays are not 1D or do not have the
            /// same size.
            pub fn axpy(a: T, x: NDArray(T), y: *NDArray(T)) !void {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").axpy(T, @intCast(x.shape[0]), a, x.data.ptr, @intCast(x.strides[0]), y.data.ptr, -@as(isize, @intCast(y.strides[0])));
            }

            /// Copies a vector to another vector.
            ///
            /// **Description**:
            ///
            /// The `copy` routine performs a vector-vector operation defined
            /// as:
            /// ```zig
            /// y = x,
            /// ```
            /// where `x` and `y` are 1D arrays with the same size.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            /// - `y`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            pub fn copy(x: NDArray(T), y: *NDArray(T)) !void {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").copy(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]));
            }

            /// Computes a vector-vector dot product.
            ///
            /// **Description**:
            ///
            /// The `dot` routines perform a vector-vector reduction operation
            /// defined as:
            /// ```zig
            /// res = x[1]*y[1] + x[2]*y[2] + ... + x[n]*y[n],
            /// ```
            /// where `x` and `y` are 1D arrays with the same size.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            /// - `y`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            pub fn dot(x: NDArray(T), y: NDArray(T)) !T {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").dot(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]));
            }

            /// Computes a dot product of a conjugated vector with another
            /// vector.
            ///
            /// **Description**:
            ///
            /// The `dotc` routine performs a vector-vector operation defined
            /// as:
            /// ```zig
            /// res = conj(x[1])*y[1] + conj(x[2])*y[2] + ... + conj(x[n])*y[n],
            /// ```
            /// where `x` and `y` are complex 1D arrays with the same size.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            /// - `y`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            pub fn dotc(x: NDArray(T), y: NDArray(T)) !T {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").dotc(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]));
            }

            /// Computes a dot product of a conjugated vector with another
            /// vector.
            ///
            /// **Description**:
            ///
            /// The `dotc` routine performs a vector-vector operation defined
            /// as:
            /// ```zig
            /// res = conj(x[1])*y[1] + conj(x[2])*y[2] + ... + conj(x[n])*y[n],
            /// ```
            /// where `x` and `y` are complex 1D arrays with the same size.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            /// - `y`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            pub fn dotc_sub(x: NDArray(T), y: NDArray(T), ret: *T) !void {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").dotc_sub(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]), ret);
            }

            /// Computes a complex array-array dot product.
            ///
            /// **Description**:
            ///
            /// The `dotu` routine performs an array-array operation defined as:
            /// ```zig
            /// res = x[1]*y[1] + x[2]*y[2] + ... + x[n]*y[n],
            /// ```
            /// where `x` and `y` are complex arrays each with the same size.
            /// Note that they may not have the same shape, strides or order,
            /// but keep in mind that the operation is performed following the
            /// one-dimensional data order.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of size `n`.
            /// - `y`: `NDArray` of size `n`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            pub fn dotu(x: NDArray(T), y: NDArray(T)) !T {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").dotu(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]));
            }

            /// Computes a complex array-array dot product.
            ///
            /// **Description**:
            ///
            /// The `dotu` routine performs an array-array operation defined as:
            /// ```zig
            /// res = x[1]*y[1] + x[2]*y[2] + ... + x[n]*y[n],
            /// ```
            /// where `x` and `y` are complex arrays each with the same size.
            /// Note that they may not have the same shape, strides or order,
            /// but keep in mind that the operation is performed following the
            /// one-dimensional data order.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of size `n`.
            /// - `y`: `NDArray` of size `n`.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            pub fn dotu_sub(x: NDArray(T), y: NDArray(T), ret: *T) !void {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").dotu_sub(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]), ret);
            }

            /// Computes the Euclidean norm of an array.
            ///
            /// **Description**:
            ///
            /// The `nrm2` routine performs and array reduction operation
            /// defined as:
            /// ```zig
            /// res = sqrt(|x[1]|^2 + |x[2]|^2 + ... + |x[n]|^2),
            /// ```
            /// where `x` is an `NDArray`.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of any shape.
            ///
            /// **Return Values**:
            /// - `T`: The Euclidean norm of all elements of the vector.
            pub fn nrm2(x: NDArray(T)) !scalar(T) {
                if (x.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").nrm2(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]));
            }

            /// Computes the Euclidean norm of an array.
            ///
            /// **Description**:
            ///
            /// The `nrm2` routine performs and array reduction operation
            /// defined as:
            /// ```zig
            /// res = sqrt(|x[1]|^2 + |x[2]|^2 + ... + |x[n]|^2),
            /// ```
            /// where `x` is an `NDArray`.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of any shape.
            ///
            /// **Return Values**:
            /// - `T`: The Euclidean norm of all elements of the vector.
            pub fn rot(x: *NDArray(T), y: *NDArray(T), c: scalar(T), s: scalar(T)) !void {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").rot(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]), c, s);
            }

            /// Performs modified Givens rotation of points in the plane.
            ///
            /// **Description**:
            ///
            /// Given two vectors `x` and `y`, each vector element of these
            /// vectors is replaced as follows:
            /// ```zig
            /// [ x[i] ]     [ x[i] ]
            /// [ y[i] ] = H [ y[i] ]
            /// ```
            /// for `i=1` to `n`, where `H` is a modified Givens transformation
            /// matrix whose values are stored in the `param[1]` through
            /// `param[4]` array. See discussion on the param argument.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            /// - `y`: `NDArray` of shape `{n}`.
            /// - `param`: `NDArray` of shape `{5}`. The last four elements
            /// contain the modified Givens transformation matrix. The first
            /// element contains a flag that determines the form of the
            /// Givens rotation.
            ///
            /// **Output Parameters**:
            /// - `x`: Overwritten by the result of the operation.
            /// - `y`: Overwritten by the result of the operation.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            /// - `IncompatibleDimensions`: the arrays are not 1D.
            pub fn rotm(x: *NDArray(T), y: *NDArray(T), param: NDArray(T)) !void {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                if (param.ndim != 1 or param.shape[0] != 5) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").rotm(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]), param.data.ptr);
            }

            /// Computes the product of a vector by a scalar.
            ///
            /// **Description**:
            ///
            /// The scal routine performs a vector operation defined as:
            /// ```zig
            /// x = a*x,
            /// ```
            /// where `a` is a scalar and `x` is a 1D array.
            ///
            /// **Input Parameters**:
            /// - `a`: scalar.
            /// - `x`: `NDArray` of shape `{n}`.
            ///
            /// **Output Parameters**:
            /// - `x`: Overwritten by the result of the operation.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleDimensions`: the array is not 1D.
            pub fn scal(a: T, x: *NDArray(T)) !void {
                if (x.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").scal(T, @intCast(x.shape[0]), a, x.data.ptr, @intCast(x.strides[0]));
            }

            /// Swaps a vector with another vector.
            ///
            /// **Description**:
            ///
            /// Given two vectors `x` and `y`, the swap routine returns vectors
            /// `y` and `x` swapped, each replacing the other.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            /// - `y`: `NDArray` of shape `{n}`.
            ///
            /// **Output Parameters**:
            /// - `x`: Overwritten by the result of the operation.
            /// - `y`: Overwritten by the result of the operation.
            ///
            /// **Return Values**:
            /// - `void`: the execution was successful.
            /// - `IncompatibleDimensions`: the arrays are not 1D.
            /// - `IncompatibleSize`: the arrays do not have the same size.
            pub fn swap(x: *NDArray(T), y: *NDArray(T)) !void {
                if (x.ndim != 1 or y.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                if (x.shape[0] != y.shape[0]) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").swap(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]), y.data.ptr, @intCast(y.strides[0]));
            }

            /// Finds the index of the element with maximum absolute value.
            ///
            /// Given a vector `x`, the iamax function returns the position of
            /// the vector element `x[i]` that has the largest absolute value
            /// for real flavors, or the largest sum `|x[i].re|+|x[i].im| for
            /// complex flavors. If more than one vector element is found with
            /// the same largest absolute value, the index of the first one
            /// encountered is returned.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `usize`: the index of the element with maximum absolute value.
            /// - `IncompatibleDimensions`: the array is not 1D.
            pub fn iamax(x: NDArray(T)) !usize {
                if (x.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").iamax(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]));
            }

            /// Finds the index of the element with minimum absolute value.
            ///
            /// Given a vector `x`, the iamin function returns the position of
            /// the vector element `x[i]` that has the smallest absolute value
            /// for real flavors, or the smallest sum `|x[i].re|+|x[i].im| for
            /// complex flavors. If more than one vector element is found with
            /// the same smallest absolute value, the index of the first one
            /// encountered is returned.
            ///
            /// **Input Parameters**:
            /// - `x`: `NDArray` of shape `{n}`.
            ///
            /// **Return Values**:
            /// - `usize`: the index of the element with minimum absolute value.
            /// - `IncompatibleDimensions`: the array is not 1D.
            pub fn iamin(x: NDArray(T)) !usize {
                if (x.ndim != 1) {
                    return Error.IncompatibleDimensions;
                }

                return @import("BLAS/BLAS.zig").iamin(T, @intCast(x.shape[0]), x.data.ptr, @intCast(x.strides[0]));
            }
        };

        /// Namespace for LAPACK functions for simple types that do not require
        /// memory management.
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
            //pub fn getrf(a: *NDArray(T), ipiv: *NDArray(usize)) !void {
            //    return @import("ndarray/LAPACK/getrf.zig").getrf(T, a, ipiv);
            //}
            pub fn tmp() usize {
                return 3;
            }
        };

        /// Namespace for BLAS functions for types that require memory
        /// management.
        pub const BLASA = struct {
            pub fn tmp() usize {
                return 3;
            }
        };

        /// Namespace for LAPACK functions for types that require memory
        /// management.
        pub const LAPACKA = struct {
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
    /// The arrays do not have the same size.
    IncompatibleSize,
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
pub const Order = enum(u7) {
    /// Row major element storage (right to left).
    RowMajor = 101,
    /// Column major element storage (left to right).
    ColumnMajor = 102,
};

/// Transposition of matrices.
pub const Transpose = enum(u7) {
    NoTrans = 111,
    Trans = 112,
    ConjTrans = 113,
    ConjNoTrans = 114,
};

/// Uplo parameter for symmetric matrices.
pub const Uplo = enum(u7) {
    Upper = 121,
    Lower = 122,
};

/// Diag parameter for triangular matrices.
pub const Diag = enum(u8) {
    NonUnit = 131,
    Unit = 132,
};

/// Side parameter for triangular matrices.
pub const Side = enum(u8) {
    Left = 141,
    Right = 142,
};

///
pub const Slice = struct {
    /// Start of the slice.
    start: usize,
    /// End of the slice, non-inclusive.
    stop: usize,
    /// Step of the slice.
    step: isize,

    /// Initializes a slice with the given parameters.
    pub fn init(start: usize, stop: usize, step: isize) !Slice {
        if (step == 0 and start != stop) {
            return Error.ZeroDimension;
        }

        if (step > 0) {
            if (start >= stop) {
                return Error.DimensionMismatch;
            }
        } else if (step < 0) {
            if (start <= stop) {
                return Error.DimensionMismatch;
            }
        }

        return Slice{ .start = start, .stop = stop, .step = step };
    }

    /// Returns the length of the slice.
    pub fn len(self: Slice) usize {
        if (self.start == self.stop) {
            return 0;
        }

        if (self.step > 0) {
            return (self.stop - self.start + @as(usize, @intCast(self.step)) - 1) / @as(usize, @intCast(self.step));
        }

        return (self.start - self.stop + @abs(self.step) - 1) / @abs(self.step);
    }
};

test "init" {
    const a: std.mem.Allocator = std.testing.allocator;

    const tooBig: [MaxDimensions + 1]usize = [_]usize{1} ** (MaxDimensions + 1);
    try std.testing.expectError(Error.TooManyDimensions, NDArray(f64).init(a, &tooBig, .{}));
    try std.testing.expectError(Error.ZeroDimension, NDArray(f64).init(a, &.{ 2, 3, 5, 6, 0, 1, 8 }, .{}));
    try std.testing.expectError(Error.InvalidFlags, NDArray(f64).init(a, &.{ 2, 2 }, .{ .ownsData = false }));

    var scalar1: NDArray(f64) = try NDArray(f64).init(a, &.{}, .{});
    defer scalar1.deinit();
    try std.testing.expect(scalar1.data.len == 1);
    try std.testing.expect(std.mem.eql(usize, &.{}, scalar1.shape[0..scalar1.ndim]));
    try std.testing.expect(std.mem.eql(isize, &.{}, scalar1.strides[0..scalar1.ndim]));
    try std.testing.expect(scalar1.size == 1);
    try std.testing.expect(scalar1.flags.order == .RowMajor);
    try std.testing.expect(scalar1.flags.ownsData);
    try std.testing.expect(scalar1.flags.writeable);

    var scalar2: NDArray(f64) = try NDArray(f64).init(a, &.{}, .{});
    defer scalar2.deinit();
    try std.testing.expect(scalar2.data.len == 1);
    try std.testing.expect(std.mem.eql(usize, &.{}, scalar2.shape[0..scalar2.ndim]));
    try std.testing.expect(std.mem.eql(isize, &.{}, scalar2.strides[0..scalar1.ndim]));
    try std.testing.expect(scalar2.size == 1);
    try std.testing.expect(scalar2.flags.order == .RowMajor);
    try std.testing.expect(scalar2.flags.ownsData);
    try std.testing.expect(scalar2.flags.writeable);

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 10, 5, 8, 3 }, .{});
    defer A.deinit();

    try std.testing.expect(A.data.len == 1200);
    try std.testing.expect(std.mem.eql(usize, &.{ 10, 5, 8, 3 }, A.shape[0..A.ndim]));
    try std.testing.expect(std.mem.eql(isize, &.{ 120, 24, 3, 1 }, A.strides[0..A.ndim]));
    try std.testing.expect(A.size == 1200);
    try std.testing.expect(A.flags.order == .RowMajor);
    try std.testing.expect(A.flags.ownsData);
    try std.testing.expect(A.flags.writeable);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 10, 5, 8, 3 }, .{ .order = .ColumnMajor });
    defer B.deinit();

    try std.testing.expect(B.data.len == 1200);
    try std.testing.expect(std.mem.eql(usize, &.{ 10, 5, 8, 3 }, B.shape[0..B.ndim]));
    try std.testing.expect(std.mem.eql(isize, &.{ 1, 10, 50, 400 }, B.strides[0..B.ndim]));
    try std.testing.expect(B.size == 1200);
    try std.testing.expect(B.flags.order == .ColumnMajor);
    try std.testing.expect(B.flags.ownsData);
    try std.testing.expect(B.flags.writeable);

    var scalar3: NDArray(f64) = try NDArray(f64).init(a, &.{}, .{});
    defer scalar3.deinit();
    try std.testing.expect(scalar3.data.len == 1);
    try std.testing.expect(std.mem.eql(usize, &.{}, scalar3.shape[0..scalar3.ndim]));
    try std.testing.expect(std.mem.eql(isize, &.{}, scalar3.strides[0..scalar3.ndim]));
    try std.testing.expect(scalar3.size == 1);
    try std.testing.expect(scalar3.flags.order == .RowMajor);
    try std.testing.expect(scalar3.flags.ownsData);
    try std.testing.expect(scalar3.flags.writeable);

    // MORE FOR VIEWS WHEN THEY ARE IMPLEMENTED
}

test "_index" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 10, 5, 8, 3 }, .{});
    defer A.deinit();

    try std.testing.expect(A._index(&.{ 3, 1, 0, 2 }) == 386);

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 10, 5, 8, 3 }, .{ .order = .ColumnMajor });
    defer B.deinit();

    try std.testing.expect(B._index(&.{ 3, 1, 0, 2 }) == 813);
}

test "_checkPosition" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 10, 5, 8, 3 }, .{});
    defer A.deinit();

    try std.testing.expectError(Error.DimensionMismatch, A._checkPosition(&.{3}));
    try std.testing.expectError(Error.PositionOutOfBounds, A._checkPosition(&.{ 3, 1, 9000, 2 }));
}

test "set" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 3, 2, 4 }, .{});
    defer A.deinit();
    // This is one way of iterating through an array.
    var elem: f64 = 1;
    for (0..A.shape[0]) |i| {
        for (0..A.shape[1]) |j| {
            for (0..A.shape[2]) |k| {
                try A.set(&.{ i, j, k }, elem);
                elem += 1;
            }
        }
    }
    try std.testing.expect(std.mem.eql(f64, A.data, &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 }));

    var B: NDArray(f64) = try NDArray(f64).init(a, &.{ 3, 2, 4 }, .{ .order = .ColumnMajor });
    defer B.deinit();
    elem = 1;
    for (0..B.shape[0]) |i| {
        for (0..B.shape[1]) |j| {
            for (0..B.shape[2]) |k| {
                try B.set(&.{ i, j, k }, elem);
                elem += 1;
            }
        }
    }
    try std.testing.expect(std.mem.eql(f64, B.data, &.{ 1, 9, 17, 5, 13, 21, 2, 10, 18, 6, 14, 22, 3, 11, 19, 7, 15, 23, 4, 12, 20, 8, 16, 24 }));
}

test "setAll" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 5, 10, 4 }, .{});
    defer A.deinit();
    // This is one way of iterating through an array.
    A.setAll(5);
    try std.testing.expect(std.mem.eql(f64, A.data, &[_]f64{5} ** 200));
}

test "replace" {
    const a: std.mem.Allocator = std.testing.allocator;

    var A: NDArray(f64) = try NDArray(f64).init(a, &.{ 3, 2, 4 }, .{});
    defer A.deinit();
    // This is one way of iterating through an array.
    var elem: f64 = 1;
    for (0..A.shape[0]) |i| {
        for (0..A.shape[1]) |j| {
            for (0..A.shape[2]) |k| {
                try A.set(&.{ i, j, k }, elem);
                elem += 1;
            }
        }
    }
    const replaced: f64 = try A.replace(&.{ 0, 1, 1 }, 20);
    try std.testing.expect(replaced == 6);
    try std.testing.expect(try A.get(&.{ 0, 1, 1 }) == 20);
}

test {
    _ = @import("BLAS/BLAS.zig");
    std.testing.refAllDeclsRecursive(@This());
}
