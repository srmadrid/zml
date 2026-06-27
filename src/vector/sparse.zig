const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const vector = @import("../vector.zig");

/// Sparse vector type, represented as a contiguous array of non-zero elements
/// of type `N` along with their corresponding indices, in ascending order.
pub fn Sparse(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.vector.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        idx: [*]usize,
        nnz: usize,
        len: usize,
        _dlen: usize, // allocated length of data
        _ilen: usize, // allocated length of idx
        flags: vector.Flags,

        // Type signatures
        pub const is_vector = true;
        pub const is_sparse = true;

        // Numeric type
        pub const Numeric = N;

        pub const empty = vector.Sparse(N){
            .data = &.{},
            .idx = &.{},
            .nnz = 0,
            .len = 0,
            ._dlen = 0,
            ._ilen = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new sparse vector with the specified length and an
        /// initial capacity for non-zero elements.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `len` (`usize`): The length of the vector.
        /// * `nnz` (`usize`): The initial capacity for non-zero elements.
        ///
        /// ## Returns
        /// `vector.Sparse(N)`: The newly initialized vector.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `vector.Error.ZeroLength`: If `len` is zero.
        /// * `vector.Error.DimensionMismatch`: If `nnz` is greater than `len`.
        pub fn init(allocator: std.mem.Allocator, len: usize, nnz: usize) !vector.Sparse(N) {
            if (len == 0)
                return vector.Error.ZeroLength;

            if (nnz > len)
                return vector.Error.DimensionMismatch;

            const data: []N = try allocator.alloc(N, nnz);
            errdefer allocator.free(data);

            return .{
                .data = data.ptr,
                .idx = (try allocator.alloc(usize, nnz)).ptr,
                .nnz = 0,
                .len = len,
                ._dlen = nnz,
                ._ilen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new `vector.Sparse(N)` with the given buffers.
        ///
        /// ## Arguments
        /// * `data_buffer` (`[]N`): The buffer for `data`.
        /// * `idx_buffer` (`[]usize`): The buffer for `idx`.
        /// * `len` (`usize`): The length of the vector.
        ///
        /// ## Returns
        /// `vector.Sparse(N)`: The newly initialized vector.
        ///
        /// ## Errors
        /// * `vector.Error.ZeroLength`: If if the length of any buffer is zero.
        pub fn initBuffer(data_buffer: []N, idx_buffer: []usize, len: usize) !vector.Sparse(N) {
            if (data_buffer.len == 0 or idx_buffer.len == 0)
                return vector.Error.ZeroLength;

            return .{
                .data = data_buffer.ptr,
                .idx = idx_buffer.ptr,
                .nnz = 0,
                .len = len,
                ._dlen = data_buffer.len,
                ._ilen = idx_buffer.len,
                .flags = .{ .owns_data = false },
            };
        }

        /// Deinitializes the vector, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Sparse(N)`): A pointer to the vector to
        ///   deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///   `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *Sparse(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.idx[0..self._ilen]);
            }

            self.* = undefined;
        }

        /// Reserves space for at least `new_nnz` non-zero elements. If `self`
        /// does not own its data or if `new_nnz` is less than or equal to the
        /// current capacity, this function is a no-op.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Sparse(N)`): A pointer to the vector to reserve
        ///   space for.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `new_nnz` (`usize`): The new capacity for non-zero elements.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `vector.Error.DimensionMismatch`: If `new_nnz` is greater than the
        ///   length of the vector.
        pub fn reserve(self: *vector.Sparse(N), allocator: std.mem.Allocator, new_nnz: usize) !void {
            if (!self.flags.owns_data)
                return;

            if (new_nnz <= self._dlen and new_nnz <= self._ilen)
                return;

            if (new_nnz > self.len)
                return vector.Error.DimensionMismatch;

            if (new_nnz > self._dlen) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnz)).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._ilen) {
                self.idx = (try allocator.realloc(self.idx[0..self._ilen], new_nnz)).ptr;
                self._ilen = new_nnz;
            }
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`vector.Sparse(N)`): A pointer to the vector to get the
        ///   element from.
        /// * `index` (`usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `vector.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn get(self: vector.Sparse(N), index: usize) !N {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            return self.getAssumeInBounds(index);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`vector.Sparse(N)`): A pointer to the vector to get
        ///   the element from.
        /// * `index` (`usize`): The index of the element to get. Assumed to be
        ///   within bounds.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: vector.Sparse(N), index: usize) N {
            var i: usize = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index)
                    return if (self.flags.noconj) self.data[i] else numeric.conj(self.data[i])
                else if (self.idx[i] > index)
                    break;
            }

            return numeric.zero(N);
        }

        /// Sets the element at the specified index, inserting it if it does not
        /// already exist and shifting elements as necessary to maintain index
        /// order.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Sparse(N)`): A pointer to the vector to set the
        ///   element in.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `index` (`usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails
        ///   when inserting a new element.
        /// * `vector.Error.PositionOutOfBounds`: If `index` is out of bounds.
        /// * `vector.Error.DataNotOwned`: If the vector does not own its data
        ///   and a resize is required.
        pub fn set(self: *vector.Sparse(N), allocator: std.mem.Allocator, index: usize, value: N) !void {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            var i: usize = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index) {
                    self.data[i] = if (self.flags.noconj) value else numeric.conj(value);

                    return;
                } else if (self.idx[i] > index) {
                    break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._ilen) {
                if (!self.flags.owns_data)
                    return vector.Error.DataNotOwned;

                // Need more space
                var new_nnz = if (self.nnz * 2 > self.len) self.len else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to the right to make space for the new element
            var j: usize = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.idx[j] = self.idx[j - 1];
            }

            self.data[i] = if (self.flags.noconj) value else numeric.conj(value);
            self.idx[i] = index;
            self.nnz += 1;
        }

        /// Sets the element at the specified index without bounds or space
        /// checking, inserting it if it does not already exist and shifting
        /// elements as necessary to maintain index order.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Sparse(N)`): A pointer to the vector to set the
        ///   element in. If needed, assumed to have enough space to insert a
        ///   new element.
        /// * `index` (`usize`): The index of the element to set. Assumed to be
        ///   within bounds.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// Returns
        /// `void`
        pub fn setAssumeInBounds(self: *vector.Sparse(N), index: usize, value: N) void {
            var i: usize = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index) {
                    self.data[i] = if (self.flags.noconj) value else numeric.conj(value);
                    return;
                } else if (self.idx[i] > index) {
                    break;
                }
            }

            // Shift elements to the right to make space for the new element
            var j: usize = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.idx[j] = self.idx[j - 1];
            }

            self.data[i] = if (self.flags.noconj) value else numeric.conj(value);
            self.idx[i] = index;
            self.nnz += 1;
        }

        /// Accumulates the specified value at the given index, adding to the
        /// existing value if it exists, or inserting it if it does not.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Sparse(N)`): A pointer to the vector to set the
        ///   element in.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `index` (`usize`): The index of the element to accumulate at.
        /// * `value` (`N`): The value to accumulate at the specified index.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails
        ///   when inserting a new element.
        /// * `vector.Error.PositionOutOfBounds`: If `index` is out of bounds.
        /// * `vector.Error.DataNotOwned`: If the vector does not own its data
        ///   and a resize is required.
        pub fn accumulate(self: *vector.Sparse(N), allocator: std.mem.Allocator, index: usize, value: N) !void {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            var i: usize = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index) {
                    numeric.addInto(
                        &self.data[i],
                        self.data[i],
                        if (self.flags.noconj)
                            value
                        else
                            numeric.conj(value),
                    );

                    return;
                } else if (self.idx[i] > index) {
                    break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._ilen) {
                if (!self.flags.owns_data)
                    return vector.Error.DataNotOwned;

                // Need more space
                var new_nnz = if (self.nnz * 2 > self.len) self.len else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to the right to make space for the new element
            var j: usize = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.idx[j] = self.idx[j - 1];
            }

            self.data[i] = if (self.flags.noconj) value else numeric.conj(value);
            self.idx[i] = index;
            self.nnz += 1;
        }
    };
}
