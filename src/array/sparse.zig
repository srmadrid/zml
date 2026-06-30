const std = @import("std");

const meta = @import("../meta.zig");

const numeric = @import("../numeric.zig");

const array = @import("../array.zig");

const arrutils = @import("utils.zig");

/// Sparse `n`-dimensional array type in CSF format, the natural generalization
/// of CSC (`order == .f`) and CSR (`order == .c`) to arbitrary dimensions.
pub fn Sparse(N: type, comptime order: array.Order) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.array.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        idx: [array.max_dimensions][*]usize,
        ptr: [array.max_dimensions][*]usize,
        nnz: usize,
        ndim: usize,
        shape: [array.max_dimensions]usize,
        _dlen: usize,
        _ilen: [array.max_dimensions]usize,
        _plen: [array.max_dimensions]usize,
        flags: array.Flags,

        // Type signatures
        pub const is_array = true;
        pub const is_sparse = true;

        // Numeric type
        pub const Numeric = N;
        pub const storage_order: array.Order = order;

        pub const empty: array.Sparse(N, order) = .{
            .data = &.{},
            .idx = .{&.{}} ** array.max_dimensions,
            .ptr = .{&.{}} ** array.max_dimensions,
            .nnz = 0,
            .ndim = 0,
            .shape = .{0} ** array.max_dimensions,
            ._dlen = 0,
            ._ilen = .{0} ** array.max_dimensions,
            ._plen = .{0} ** array.max_dimensions,
            .flags = .{ .owns_data = false },
        };

        /// Deinitializes the array, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*array.Sparse(N, order)`): A pointer to the array to
        ///   deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to compile `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *array.Sparse(N, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);

                for (0..self.ndim) |k| {
                    if (self._ilen[k] > 0)
                        allocator.free(self.idx[k][0..self._ilen[k]]);
                }

                if (self.ndim >= 2) {
                    for (0..self.ndim - 1) |k| {
                        if (self._plen[k] > 0)
                            allocator.free(self.ptr[k][0..self._plen[k]]);
                    }
                }
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`array.Sparse(N, order)`): The array to get the element
        ///   from.
        /// * `index` (`[]const usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn get(self: array.Sparse(N, order), index: []const usize) !N {
            try arrutils.checkIndex(self.shape[0..self.ndim], index);

            return self.getAssumeInBounds(index);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`array.Sparse(N, order)`): The array to get the element
        ///   from.
        /// * `index` (`[]const usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: array.Sparse(N, order), index: []const usize) N {
            if (self.nnz == 0)
                return numeric.zero(N);

            var range_start: usize = 0;
            var range_end: usize = self._ilen[0];

            var level: usize = 0;
            while (level < self.ndim) : (level += 1) {
                const mode: usize = if (comptime order == .c) level else self.ndim - 1 - level;
                const target = index[mode];

                var lo = range_start;
                var hi = range_end;
                var found = false;
                var found_pos: usize = 0;

                while (lo < hi) {
                    const mid = lo + (hi - lo) / 2;
                    const mv = self.idx[level][mid];
                    if (mv == target) {
                        found = true;
                        found_pos = mid;
                        break;
                    } else if (mv < target) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }

                if (!found)
                    return numeric.zero(N);

                if (level == self.ndim - 1)
                    return if (self.flags.noconj) self.data[found_pos] else numeric.conj(self.data[found_pos]);

                range_start = self.ptr[level][found_pos];
                range_end = self.ptr[level][found_pos + 1];
            }

            return numeric.zero(N);
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*array.Sparse(N, order)`): A pointer to the array to set
        ///   the element in.
        /// * `index` (`[]const usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        /// * `array.Error.BreaksStructure`: If the position is not already
        ///   stored as a non-zero.
        pub fn set(self: *array.Sparse(N, order), index: []const usize, value: N) !void {
            try arrutils.checkIndex(self.shape[0..self.ndim], index);

            return self.setAssumeInBounds(index, value);
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*array.Sparse(N, order)`): A pointer to the array to set
        ///   the element in.
        /// * `index` (`[]const usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *array.Sparse(N, order), index: []const usize, value: N) !void {
            if (self.nnz == 0)
                return array.Error.BreaksStructure;

            var range_start: usize = 0;
            var range_end: usize = self._ilen[0];

            var level: usize = 0;
            while (level < self.ndim) : (level += 1) {
                const mode: usize = if (comptime order == .c) level else self.ndim - 1 - level;
                const target = index[mode];

                var lo = range_start;
                var hi = range_end;
                var found = false;
                var found_pos: usize = 0;

                while (lo < hi) {
                    const mid = lo + (hi - lo) / 2;
                    const mv = self.idx[level][mid];
                    if (mv == target) {
                        found = true;
                        found_pos = mid;
                        break;
                    } else if (mv < target) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }

                if (!found)
                    return array.Error.BreaksStructure;

                if (level == self.ndim - 1) {
                    self.data[found_pos] = if (self.flags.noconj) value else numeric.conj(value);
                    return;
                }

                range_start = self.ptr[level][found_pos];
                range_end = self.ptr[level][found_pos + 1];
            }

            return array.Error.BreaksStructure;
        }
    };
}
