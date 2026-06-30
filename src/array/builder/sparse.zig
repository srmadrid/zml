const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const array = @import("../../array.zig");

const arrutils = @import("../utils.zig");

/// Sparse builder array type, represented in COO format. Two arrays are used
/// to store the indices and values of the non-zero entries. Indices are not
/// sorted and duplicate entries are allowed, getting summed at compilation.
/// This type cannot be used for array computations directly; it must first be
/// compiled into a standard sparse array.
pub fn Sparse(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.array.builder.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        idx: [array.max_dimensions][*]usize,
        nnz: usize,
        ndim: usize,
        shape: [array.max_dimensions]usize,
        _dlen: usize,
        _ilen: usize,

        // Type signatures
        pub const is_array = true;
        pub const is_builder = true;
        pub const is_sparse = true;

        // Numeric type
        pub const Numeric = N;

        pub const empty: array.builder.Sparse(N) = .{
            .data = &.{},
            .idx = .{&.{}} ** array.max_dimensions,
            .nnz = 0,
            .ndim = 0,
            .shape = .{0} ** array.max_dimensions,
            ._dlen = 0,
            ._ilen = 0,
        };

        /// Initializes a new `array.builder.Sparse(N)` with the specified shape
        /// and an initial capacity for non-zero entries.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for
        ///   memory allocations.
        /// * `shape` (`[]const usize`): The shape of the array.
        /// * `nnz` (`usize`): Initial capacity for non-zero entries.
        ///
        /// ## Returns
        /// `array.builder.Sparse(N)`: The newly initialized builder array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `array.Error.TooManyDimensions`: If `shape.len` is larger than
        ///   `array.max_dimensions`.
        /// * `array.Error.ZeroDimension`: If `shape.len` or any dimension is
        ///   zero.
        pub fn init(allocator: std.mem.Allocator, shape: []const usize, nnz: usize) !array.builder.Sparse(N) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len == 0)
                return array.Error.ZeroDimension;

            for (shape) |dim| {
                if (dim == 0)
                    return array.Error.ZeroDimension;
            }

            const ndim = shape.len;
            var builder_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            for (0..ndim) |d| {
                builder_shape[d] = shape[d];
            }

            const data: []N = try allocator.alloc(N, nnz);
            errdefer allocator.free(data);

            var idx: [array.max_dimensions][*]usize = undefined;
            var n_alloc: usize = 0;
            errdefer for (0..n_alloc) |d| allocator.free(idx[d][0..nnz]);

            for (0..ndim) |d| {
                idx[d] = (try allocator.alloc(usize, nnz)).ptr;
                n_alloc += 1;
            }

            return .{
                .data = data.ptr,
                .idx = idx,
                .nnz = 0,
                .ndim = ndim,
                .shape = builder_shape,
                ._dlen = nnz,
                ._ilen = nnz,
            };
        }

        // initBuffer

        /// Deinitializes the builder array, freeing all allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*array.builder.Sparse(N)`): A pointer to the builder
        ///   array to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///   `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *array.builder.Sparse(N), allocator: std.mem.Allocator) void {
            allocator.free(self.data[0..self._dlen]);
            for (0..self.ndim) |d| {
                allocator.free(self.idx[d][0..self._ilen]);
            }

            self.* = undefined;
        }

        /// Reserves capacity for at least `new_nnz` non-zero entries.
        ///
        /// ## Arguments
        /// * `self` (`*array.builder.Sparse(N)`): A pointer to the builder
        ///   array to reserve space for.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `new_nnz` (`usize`): The new capacity for non-zero entries.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn reserve(self: *array.builder.Sparse(N), allocator: std.mem.Allocator, new_nnz: usize) !void {
            if (new_nnz <= self._dlen and new_nnz <= self._ilen)
                return;

            if (new_nnz > self._dlen) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnz)).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._ilen) {
                for (0..self.ndim) |d| {
                    self.idx[d] = (try allocator.realloc(self.idx[d][0..self._ilen], new_nnz)).ptr;
                }

                self._ilen = new_nnz;
            }
        }

        /// Appends a new non-zero entry to the builder array.
        ///
        /// ## Arguments
        /// * `self` (`*array.builder.Sparse(N)`): A pointer to the builder
        ///   array.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `index` (`[]const usize`): The index of the entry.
        /// * `value` (`N`): The value to insert.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails
        ///   when resizing the internal buffers.
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn append(self: *array.builder.Sparse(N), allocator: std.mem.Allocator, index: []const usize, value: N) !void {
            try arrutils.checkIndex(self.shape[0..self.ndim], index);

            if (self.nnz == self._dlen or self.nnz == self._ilen) {
                var new_nnz = self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            return self.appendAssumeCapacity(index, value);
        }

        /// Appends a new non-zero element without performing bounds checks
        /// or verifying capacity.
        ///
        /// ## Arguments
        /// * `self` (`*array.builder.Sparse(N)`): A pointer to the builder
        ///   array.
        /// * `index` (`[]const usize`): The index of the entry. Assumed to be
        ///   within bounds.
        /// * `value` (`N`): The value to insert.
        ///
        /// ## Returns
        /// `void`
        pub fn appendAssumeCapacity(self: *array.builder.Sparse(N), index: []const usize, value: N) void {
            self.data[self.nnz] = value;
            for (0..self.ndim) |d| {
                self.idx[d][self.nnz] = index[d];
            }
            self.nnz += 1;
        }

        /// Creates a copy of the builder.
        ///
        /// ## Arguments
        /// * `self` (`array.builder.Sparse(N)`): The builder array to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `array.builder.Sparse(N)`: The copied builder array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn clone(self: array.builder.Sparse(N), allocator: std.mem.Allocator) !array.builder.Sparse(N) {
            const data: []N = try allocator.alloc(N, self.nnz);
            errdefer allocator.free(data);
            var idx: [array.max_dimensions][*]usize = undefined;
            var n_alloc: usize = 0;
            errdefer for (0..n_alloc) |d| allocator.free(idx[d][0..self.nnz]);

            for (0..self.ndim) |d| {
                idx[d] = (try allocator.alloc(usize, self.nnz)).ptr;
                n_alloc += 1;
            }

            for (0..self.nnz) |i| {
                data[i] = self.data[i];
                for (0..self.ndim) |d| {
                    idx[d][i] = self.idx[d][i];
                }
            }

            return .{
                .data = data.ptr,
                .idx = idx,
                .nnz = self.nnz,
                .ndim = self.ndim,
                .shape = self.shape,
                ._dlen = self.nnz,
                ._ilen = self.nnz,
            };
        }

        /// Compiles the builder array into an `array.Sparse(N, order)`,
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// ## Arguments
        /// * `self` (`*array.builder.Sparse(N)`): A pointer to the builder
        ///   array to compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `order` (`comptime array.Order`): Storage order of the compiled
        ///   array. `.f` makes mode `ndim-1` the outermost fiber level
        ///   (CSC-like); `.c` makes mode 0 the outermost (CSR-like).
        ///
        /// ## Returns
        /// `array.Sparse(N, order)`: The compiled sparse array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn compile(self: *array.builder.Sparse(N), allocator: std.mem.Allocator, comptime order: array.Order) !array.Sparse(N, order) {
            return self.compileInternal(allocator, order);
        }

        /// Compiles the builder array into an `array.Sparse(N, order)` by
        /// copying the data, leaving the builder intact.
        ///
        /// ## Arguments
        /// * `self` (`array.builder.Sparse(N)`): The builder array to compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `order` (`comptime array.Order`): Storage order of the compiled
        ///   array. `.f` makes mode `ndim-1` the outermost fiber level
        ///   (CSC-like); `.c` makes mode 0 the outermost (CSR-like).
        ///
        /// ## Returns
        /// `array.Sparse(N, order)`: The compiled sparse array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn compileClone(self: array.builder.Sparse(N), allocator: std.mem.Allocator, comptime order: array.Order) !array.Sparse(N, order) {
            var copy = try self.clone(allocator);
            errdefer copy.deinit(allocator);
            return copy.compileInternal(allocator, order);
        }

        fn compileInternal(self: *array.builder.Sparse(N), allocator: std.mem.Allocator, comptime order: array.Order) !array.Sparse(N, order) {
            if (self.nnz == 0) {
                allocator.free(self.data[0..self._dlen]);
                for (0..self.ndim) |d| allocator.free(self.idx[d][0..self._ilen]);

                self.* = undefined;

                return .{
                    .data = &.{},
                    .idx = .{&[_]usize{}} ** array.max_dimensions,
                    .ptr = .{&[_]usize{}} ** array.max_dimensions,
                    .nnz = 0,
                    .ndim = self.ndim,
                    .shape = self.shape,
                    ._dlen = 0,
                    ._ilen = .{0} ** array.max_dimensions,
                    ._plen = .{0} ** array.max_dimensions,
                    .flags = .{ .owns_data = false },
                };
            }

            const perm: []usize = try allocator.alloc(usize, self.nnz);
            defer allocator.free(perm);

            var i: usize = 0;
            while (i < perm.len) : (i += 1) {
                perm[i] = i;
            }

            const Context = struct {
                idx: *const [array.max_dimensions][*]usize,
                ndim: usize,

                pub fn lessThan(ctx: @This(), a: usize, b: usize) bool {
                    var level: usize = 0;
                    while (level < ctx.ndim) : (level += 1) {
                        const mode: usize = if (comptime order == .c) level else ctx.ndim - 1 - level;
                        const ia = ctx.idx[mode][a];
                        const ib = ctx.idx[mode][b];

                        if (ia < ib)
                            return true;

                        if (ia > ib)
                            return false;
                    }

                    return false;
                }
            };

            std.mem.sortUnstable(usize, perm, Context{ .idx = &self.idx, .ndim = self.ndim }, Context.lessThan);

            var unique_nnz: usize = 1;
            var fiber_count: [array.max_dimensions]usize = .{0} ** array.max_dimensions;

            if (self.ndim >= 2) {
                for (0..self.ndim - 1) |k| {
                    fiber_count[k] = 1;
                }
            }

            var last_idx: [array.max_dimensions]usize = undefined;
            for (0..self.ndim) |d| {
                last_idx[d] = self.idx[d][perm[0]];
            }

            i = 1;
            while (i < self.nnz) : (i += 1) {
                var fdl: usize = self.ndim;
                var level: usize = 0;
                while (level < self.ndim) : (level += 1) {
                    const mode: usize = if (comptime order == .c) level else self.ndim - 1 - level;
                    if (self.idx[mode][perm[i]] != last_idx[mode]) {
                        fdl = level;
                        break;
                    }
                }

                if (fdl == self.ndim)
                    continue;

                unique_nnz += 1;

                if (self.ndim >= 2) {
                    for (fdl..self.ndim - 1) |k| {
                        fiber_count[k] += 1;
                    }
                }

                for (0..self.ndim) |d| {
                    last_idx[d] = self.idx[d][perm[i]];
                }
            }

            var idx: [array.max_dimensions][*]usize = undefined;
            var ptr: [array.max_dimensions][*]usize = undefined;
            var _ilen: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var _plen: [array.max_dimensions]usize = .{0} ** array.max_dimensions;

            var n_idx_alloc: usize = 0;
            var n_ptr_alloc: usize = 0;

            errdefer {
                for (0..n_idx_alloc) |k| allocator.free(idx[k][0.._ilen[k]]);
                for (0..n_ptr_alloc) |k| allocator.free(ptr[k][0.._plen[k]]);
            }

            if (self.ndim >= 2) {
                for (0..self.ndim - 1) |k| {
                    const fc = fiber_count[k];

                    const idx_s = try allocator.alloc(usize, fc);
                    idx[k] = idx_s.ptr;
                    _ilen[k] = fc;
                    n_idx_alloc += 1;

                    const ptr_s = try allocator.alloc(usize, fc + 1);
                    @memset(ptr_s, 0);
                    ptr[k] = ptr_s.ptr;
                    _plen[k] = fc + 1;
                    n_ptr_alloc += 1;
                }
            }

            const leaf_s = try allocator.alloc(usize, unique_nnz);
            idx[self.ndim - 1] = leaf_s.ptr;
            _ilen[self.ndim - 1] = unique_nnz;
            n_idx_alloc += 1;

            const data = try allocator.alloc(N, unique_nnz);
            errdefer allocator.free(data);

            const leaf_mode: usize = if (comptime order == .c) self.ndim - 1 else 0;

            var fiber_pos: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var write_pos: usize = 0;

            var acc_val: N = undefined;
            var acc_idx: [array.max_dimensions]usize = undefined;
            var has_acc = false;

            i = 0;
            while (i < self.nnz) : (i += 1) {
                const p = perm[i];

                var entry_idx: [array.max_dimensions]usize = undefined;
                for (0..self.ndim) |d| {
                    entry_idx[d] = self.idx[d][p];
                }

                var fdl: usize = if (!has_acc) 0 else self.ndim;
                if (has_acc) {
                    var level: usize = 0;
                    while (level < self.ndim) : (level += 1) {
                        const mode: usize = if (comptime order == .c) level else self.ndim - 1 - level;
                        if (entry_idx[mode] != acc_idx[mode]) {
                            fdl = level;
                            break;
                        }
                    }
                }

                if (fdl == self.ndim) {
                    numeric.addInto(&acc_val, acc_val, self.data[p]);
                    continue;
                }

                if (has_acc) {
                    idx[self.ndim - 1][write_pos] = acc_idx[leaf_mode];
                    data[write_pos] = acc_val;
                    write_pos += 1;

                    if (self.ndim >= 2)
                        ptr[self.ndim - 2][fiber_pos[self.ndim - 2]] += 1;
                }

                if (self.ndim >= 2) {
                    for (fdl..self.ndim - 1) |k| {
                        const mode: usize = if (comptime order == .c) k else self.ndim - 1 - k;
                        idx[k][fiber_pos[k]] = entry_idx[mode];
                        fiber_pos[k] += 1;
                        if (k > 0)
                            ptr[k - 1][fiber_pos[k - 1]] += 1;
                    }
                }

                acc_val = self.data[p];
                acc_idx = entry_idx;
                has_acc = true;
            }

            if (has_acc) {
                idx[self.ndim - 1][write_pos] = acc_idx[leaf_mode];
                data[write_pos] = acc_val;
                write_pos += 1;
                if (self.ndim >= 2)
                    ptr[self.ndim - 2][fiber_pos[self.ndim - 2]] += 1;
            }

            std.debug.assert(write_pos == unique_nnz);

            if (self.ndim >= 2) {
                for (0..self.ndim - 1) |k| {
                    const fc = fiber_count[k];
                    var j: usize = 0;
                    while (j < fc) : (j += 1) {
                        ptr[k][j + 1] += ptr[k][j];
                    }
                }
            }

            allocator.free(self.data[0..self._dlen]);
            for (0..self.ndim) |d| allocator.free(self.idx[d][0..self._ilen]);

            var arr: array.Sparse(N, order) = undefined;
            arr.data = data.ptr;
            arr._dlen = unique_nnz;
            arr.idx = idx;
            arr._ilen = _ilen;
            arr.ptr = ptr;
            arr._plen = _plen;
            arr.nnz = unique_nnz;
            arr.flags = .{};

            arr.shape = self.shape;
            arr.ndim = self.ndim;

            self.* = undefined;

            return arr;
        }
    };
}
