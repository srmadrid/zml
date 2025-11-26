//! Storage scheme:
//!
//! COO (Coordinate List) format. Order chooses ordering of (row, col, data)
//! arrays:
//! - row_major: (row, col, data) sorted by row, then by col -> compiles to CSR
//! - col_major: (col, row, data) sorted by col, then by row -> compiles to CSC

const std = @import("std");

const types = @import("../../types.zig");
const Order = types.Order;
const Uplo = types.Uplo;
const Diag = types.Diag;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

/// Sparse builder matrix type, represented in COO format. Three arrays are
/// used to store the row indices, column indices, and values of the non-zero
/// elements. The arrays are sorted according to the specified order (row major
/// means sorted by row, then column; column major means sorted by column, then
/// row). This type cannot be used for matrix computations directly; it must be
/// first compiled into a standard sparse matrix.
pub fn Sparse(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("matrix.builder.Sparse requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        row: [*]u32,
        col: [*]u32,
        nnz: u32,
        rows: u32,
        cols: u32,
        _dlen: u32, // allocated length of data
        _rlen: u32, // allocated length of rows
        _clen: u32, // allocated length of cols
        flags: Flags = .{},

        pub const empty = Sparse(T, order){
            .data = &.{},
            .row = &.{},
            .col = &.{},
            .nnz = 0,
            .rows = 0,
            .cols = 0,
            ._dlen = 0,
            ._rlen = 0,
            ._clen = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new builder matrix with the specified rows and
        /// columns, and an initial capacity for non-zero elements.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `rows` (`u32`):
        /// The rows of the builder matrix.
        ///
        /// `cols` (`u32`):
        /// The columns of the builder matrix.
        ///
        /// `nnz` (`u32`):
        /// The initial capacity for non-zero elements.
        ///
        /// Returns
        /// -------
        /// `matrix.builder.Sparse(T, order)`:
        /// The newly initialized builder matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.ZeroDimension`:
        /// If either `rows` or `cols` is zero.
        ///
        /// `matrix.Error.DimensionMismatch`:
        /// If `nnz` is zero or greater than `rows * cols`.
        pub fn init(allocator: std.mem.Allocator, rows: u32, cols: u32, nnz: u32) !Sparse(T, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            if (nnz == 0 or nnz > rows * cols)
                return matrix.Error.DimensionMismatch;

            const data: []T = try allocator.alloc(T, nnz);
            errdefer allocator.free(data);

            const row: []u32 = try allocator.alloc(u32, nnz);
            errdefer allocator.free(row);

            return .{
                .data = data.ptr,
                .row = row.ptr,
                .col = (try allocator.alloc(u32, nnz)).ptr,
                .nnz = 0,
                .rows = rows,
                .cols = cols,
                ._dlen = nnz,
                ._rlen = nnz,
                ._clen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        /// Deinitializes the builder matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to deinitialize.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory deallocation. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, `cleanup` must be
        /// called before `deinit` to properly deinitialize the elements.
        pub fn deinit(self: *Sparse(T, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.row[0..self._rlen]);
                allocator.free(self.col[0..self._clen]);
            }

            self.* = undefined;
        }

        /// Reserves space for at least `new_nnz` non-zero elements.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `new_nnz` (`u32`):
        /// The new capacity for non-zero elements.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.DimensionMismatch`:
        /// If `new_nnz` is greater than `rows * cols`.
        ///
        /// Notes
        /// -----
        /// If `self` does not own its data, this function does nothing, and
        /// each array is only reallocated if `new_nnz` is greater than its
        /// current capacity.
        pub fn reserve(self: *Sparse(T, order), allocator: std.mem.Allocator, new_nnz: u32) !void {
            if (!self.flags.owns_data)
                return;

            if (new_nnz <= self._dlen and new_nnz <= self._rlen and new_nnz <= self._clen)
                return;

            if (new_nnz > self.rows * self.cols)
                return matrix.Error.DimensionMismatch;

            if (new_nnz > self._dlen) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnz)).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._rlen) {
                self.row = (try allocator.realloc(self.row[0..self._rlen], new_nnz)).ptr;
                self._rlen = new_nnz;
            }

            if (new_nnz > self._clen) {
                self.col = (try allocator.realloc(self.col[0..self._clen], new_nnz)).ptr;
                self._clen = new_nnz;
            }
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `r` or `c` is out of bounds.
        pub fn get(self: *const Sparse(T, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c)
                    return self.data[i];

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be within bounds.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be within bounds.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub fn at(self: *Sparse(T, order), r: u32, c: u32) T {
            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c)
                    return self.data[i];

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        /// Sets the element at the specified location, inserting it if it does
        /// not already exist and shifting elements as necessary to maintain
        /// index order.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T)`):
        /// A pointer to the builder matrix to set the element in.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set.
        ///
        /// `value` (`T`):
        /// The value to set the element to.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails when inserting a new element.
        ///
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `r` or `c` is out of bounds.
        ///
        /// `matrix.Error.DataNotOwned`:
        /// If the builder matrix does not own its data and a resize is
        /// required.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type and an existing
        /// element is being overwritten at `(r, c)`, the existing is not
        /// deinitialized. The user must ensure that no memory leaks occur.
        ///
        /// If the elements are of arbitrary precision type, the builder matrix
        /// takes ownership of `value`.
        pub fn set(self: *Sparse(T, order), allocator: std.mem.Allocator, r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c) {
                    self.data[i] = value;
                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._rlen or self.nnz == self._clen) {
                if (!self.flags.owns_data)
                    return matrix.Error.DataNotOwned;
                // Need more space
                var new_nnz = if (self.nnz * 2 > self.rows * self.cols) self.rows * self.cols else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i] = value;
            self.row[i] = r;
            self.col[i] = c;
            self.nnz += 1;
        }

        /// Sets the element at the specified location without bounds or space
        /// checking, inserting it if it does not already exist and shifting
        /// elements as necessary to maintain index order.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrxi.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to set the element in. If needed,
        /// assumed to have enough space to insert a new element.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be within bounds.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be within bounds.
        ///
        /// `value` (`T`):
        /// The value to set the element to.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Notes
        /// -----
        /// Attempting to set a location that does not exist or when there is no
        /// space will result in undefined behavior.
        ///
        /// If the elements are of arbitrary precision type and an existing
        /// element is being overwritten at `(r, c)`, the existing is not
        /// deinitialized. The user must ensure that no memory leaks occur.
        ///
        /// If the elements are of arbitrary precision type, the builder matrix
        /// takes ownership of `value`.
        pub fn put(self: *Sparse(T, order), r: u32, c: u32, value: T) void {
            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c) {
                    self.data[i] = value;
                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i] = value;
            self.row[i] = r;
            self.col[i] = c;
            self.nnz += 1;
        }

        /// Accumulates the specified value at the given location, adding to the
        /// existing value if it exists, or inserting it if it does not.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T)`):
        /// A pointer to the builder matrix to set the element in.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `r` (`u32`):
        /// The row index of the element to accumulate at.
        ///
        /// `c` (`u32`):
        /// The column index of the element to accumulate at.
        ///
        /// `value` (`T`):
        /// The value to accumulate at the specified index.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails when inserting a new element.
        ///
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `r` or `c` is out of bounds.
        ///
        /// `matrix.Error.DataNotOwned`:
        /// If the builder matrix does not own its data and a resize is required.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type and the element is
        /// being inserted, the builder matrix takes ownership of `value`.
        ///
        /// When `T` is of arbitrary precision, the context may provide an
        /// optional pre-allocated buffer to store intermediate results of the
        /// addition, avoiding repeated allocations in scenarios where
        /// `accumulate` is called multiple times. If no buffer is provided, the
        /// operation will allocate a temporary buffer internally, using the
        /// allocator specified in the context.
        pub fn accumulate(self: *Sparse(T, order), allocator: std.mem.Allocator, r: u32, c: u32, value: anytype, ctx: anytype) !void {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                            .buffer = .{ .type = ?*T, .required = false, .default = null },
                        },
                    );
                },
            };

            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c) {
                    try ops.add_(
                        &self.data[i],
                        self.data[i],
                        value,
                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                    );

                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._rlen or self.nnz == self._clen) {
                if (!self.flags.owns_data)
                    return matrix.Error.DataNotOwned;

                // Need more space
                var new_nnz = if (self.nnz * 2 > self.rows * self.cols) self.rows * self.cols else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i] = value;
            self.row[i] = r;
            self.col[i] = c;
            self.nnz += 1;
        }

        /// Creates a copy of the builder matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to copy.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.builder.Sparse(T, order)`:
        /// The copied builder matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, they are deep
        /// copied.
        pub fn copy(self: *const Sparse(T, order), allocator: std.mem.Allocator, ctx: anytype) !Sparse(T, order) {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            var data: []T = try allocator.alloc(T, self.nnz);
            errdefer allocator.free(data);
            var row: []u32 = try allocator.alloc(u32, self.nnz);
            errdefer allocator.free(row);
            var col: []u32 = try allocator.alloc(u32, self.nnz);
            errdefer allocator.free(col);

            var i: u32 = 0;

            errdefer _cleanup(data.ptr, i, ctx);

            while (i < self.nnz) : (i += 1) {
                data[i] = try ops.copy(
                    self.data[i],
                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                );
                row[i] = self.row[i];
                col[i] = self.col[i];
            }

            return .{
                .data = data.ptr,
                .row = row.ptr,
                .col = col.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                ._dlen = self.nnz,
                ._rlen = self.nnz,
                ._clen = self.nnz,
                .flags = .{ .owns_data = true },
            };
        }

        /// Compiles the builder matrix into a general sparse matrix,
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// Returns
        /// -------
        /// `matrix.general.Sparse(T, order)`:
        /// The compiled general sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        pub fn compile(self: *Sparse(T, order), allocator: std.mem.Allocator) !matrix.general.Sparse(T, order) {
            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.general.Sparse(T, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        /// Compiles the builder matrix into a general sparse matrix by copying
        /// the data, leaving the builder intact.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.general.Sparse(T, order)`:
        /// The compiled general sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        pub fn compileCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.general.Sparse(T, order) {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            var data: []T = try allocator.alloc(T, self.nnz);
            errdefer allocator.free(data);
            var idx: []u32 = try allocator.alloc(u32, self.nnz);
            errdefer allocator.free(idx);

            i = 0;

            errdefer _cleanup(data.ptr, i, ctx);

            while (i < self.nnz) : (i += 1) {
                data[i] = try ops.copy(
                    self.data[i],
                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                );
                idx[i] = if (comptime order == .col_major) self.row[i] else self.col[i];
            }

            const result = matrix.general.Sparse(T, order){
                .data = data.ptr,
                .idx = idx.ptr,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        fn removeTriangle(self: *Sparse(T, order), comptime uplo: Uplo, comptime diagonal: bool, ctx: anytype) void {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            var i: u32 = 0;
            var j: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                const r = self.row[i];
                const c = self.col[i];

                const keep = switch (comptime uplo) {
                    .upper => if (diagonal) r > c else r >= c,
                    .lower => if (diagonal) r < c else r <= c,
                };

                if (comptime types.isArbitraryPrecision(T)) {
                    if (!keep) {
                        // Deinitialize element
                        ops.deinit(
                            &self.data[i],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }

                if (keep) {
                    if (i != j) {
                        self.data[j] = self.data[i];
                        self.row[j] = self.row[i];
                        self.col[j] = self.col[i];
                    }

                    j += 1;
                }
            }

            self.nnz = j;
        }

        /// Compiles the builder matrix into a symmetric sparse matrix, keeping
        /// only the specified triangle part and discarding the other, and
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies which triangle part to keep.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.symmetric.Sparse(T, uplo, order)`:
        /// The compiled symmetric sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.NotSquare`:
        /// If the builder matrix is not square.
        pub fn compileSymmetric(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, ctx: anytype) !matrix.symmetric.Sparse(T, uplo, order) {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            self.removeTriangle(comptime uplo.invert(), false, ctx);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.symmetric.Sparse(T, uplo, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .size = self.rows,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        /// Compiles the builder matrix into a symmetric sparse matrix by
        /// copying the data, leaving the builder intact, keeping only the
        /// specified triangle part and discarding the other.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies which triangle part to keep.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.symmetric.Sparse(T, uplo, order)`:
        /// The compiled symmetric sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.NotSquare`:
        /// If the builder matrix is not square.
        pub fn compileSymmetricCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, ctx: anytype) !matrix.symmetric.Sparse(T, uplo, order) {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            // May copy unneeded data. Eventually optimize this
            var cpy = try self.copy(allocator, ctx);
            errdefer cpy.deinit(allocator);
            errdefer cpy.cleanup(ctx);

            cpy.removeTriangle(comptime uplo.invert(), false, ctx);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) cpy.cols + 1 else cpy.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < cpy.nnz and cpy.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < cpy.nnz and cpy.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(cpy.col[0..cpy._clen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._rlen > cpy.nnz)
                    cpy.row = (try allocator.realloc(cpy.row[0..cpy._rlen], cpy.nnz)).ptr;
            } else {
                allocator.free(cpy.row[0..cpy._rlen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._clen > cpy.nnz)
                    cpy.col = (try allocator.realloc(cpy.col[0..cpy._clen], cpy.nnz)).ptr;
            }

            const result = matrix.symmetric.Sparse(T, uplo, order){
                .data = cpy.data,
                .idx = if (comptime order == .col_major) cpy.row else cpy.col,
                .ptr = ptr.ptr,
                .nnz = cpy.nnz,
                .size = cpy.rows,
                .flags = .{ .owns_data = true },
            };

            cpy.* = undefined;

            return result;
        }

        /// Compiles the builder matrix into a hermitian sparse matrix, keeping
        /// only the specified triangle part and discarding the other, and
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies which triangle part to keep.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.hermitian.Sparse(T, uplo, order)`:
        /// The compiled hermitian sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.NotSquare`:
        /// If the builder matrix is not square.
        pub fn compileHermitian(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, ctx: anytype) !matrix.hermitian.Sparse(T, uplo, order) {
            if (comptime !types.isComplex(T))
                @compileError("T must be a complex type");

            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            self.removeTriangle(comptime uplo.invert(), false, ctx);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.hermitian.Sparse(T, uplo, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .size = self.rows,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        /// Compiles the builder matrix into a hermitian sparse matrix by
        /// copying the data, leaving the builder intact, keeping only the
        /// specified triangle part and discarding the other.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies which triangle part to keep.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.hermitian.Sparse(T, uplo, order)`:
        /// The compiled hermitian sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.NotSquare`:
        /// If the builder matrix is not square.
        pub fn compileHermitianCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, ctx: anytype) !matrix.hermitian.Sparse(T, uplo, order) {
            if (comptime !types.isComplex(T))
                @compileError("T must be a complex type");

            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            var cpy = try self.copy(allocator, ctx);
            errdefer cpy.deinit(allocator);
            errdefer cpy.cleanup(ctx);

            cpy.removeTriangle(comptime uplo.invert(), false);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) cpy.cols + 1 else cpy.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < cpy.nnz and cpy.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < cpy.nnz and cpy.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(cpy.col[0..cpy._clen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._rlen > cpy.nnz)
                    cpy.row = (try allocator.realloc(cpy.row[0..cpy._rlen], cpy.nnz)).ptr;
            } else {
                allocator.free(cpy.row[0..cpy._rlen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._clen > cpy.nnz)
                    cpy.col = (try allocator.realloc(cpy.col[0..cpy._clen], cpy.nnz)).ptr;
            }

            const result = matrix.hermitian.Sparse(T, uplo, order){
                .data = cpy.data,
                .idx = if (comptime order == .col_major) cpy.row else cpy.col,
                .ptr = ptr.ptr,
                .nnz = cpy.nnz,
                .size = cpy.rows,
                .flags = .{ .owns_data = true },
            };

            cpy.* = undefined;

            return result;
        }

        /// Compiles the builder matrix into a triangular sparse matrix, keeping
        /// only the specified triangle part and discarding the other, and
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies which triangle part to keep.
        ///
        /// `diag` (`Diag`):
        /// Specifies whether the diagonal is unit or non-unit.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.triangular.Sparse(T, uplo, diag, order)`:
        /// The compiled triangular sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        pub fn compileTriangular(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime diag: Diag, ctx: anytype) !matrix.triangular.Sparse(T, uplo, diag, order) {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            self.removeTriangle(comptime uplo.invert(), diag == .unit, ctx);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.triangular.Sparse(T, uplo, diag, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        /// Compiles the builder matrix into a triangular sparse matrix by
        /// copying the data, leaving the builder intact, keeping only the
        /// specified triangle part and discarding the other.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to compile.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies which triangle part to keep.
        ///
        /// `diag` (`Diag`):
        /// Specifies whether the diagonal is unit or non-unit.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.triangular.Sparse(T, uplo, diag, order)`:
        /// The compiled triangular sparse matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        pub fn compileTriangularCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime diag: Diag, ctx: anytype) !matrix.triangular.Sparse(T, uplo, diag, order) {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            var cpy = try self.copy(allocator, ctx);
            errdefer cpy.deinit(allocator);
            errdefer cpy.cleanup(ctx);

            cpy.removeTriangle(comptime uplo.invert(), diag == .unit, ctx);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) cpy.cols + 1 else cpy.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < cpy.nnz and cpy.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < cpy.nnz and cpy.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(cpy.col[0..cpy._clen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._rlen > cpy.nnz)
                    cpy.row = (try allocator.realloc(cpy.row[0..cpy._rlen], cpy.nnz)).ptr;
            } else {
                allocator.free(cpy.row[0..cpy._rlen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._clen > cpy.nnz)
                    cpy.col = (try allocator.realloc(cpy.col[0..cpy._clen], cpy.nnz)).ptr;
            }

            const result = matrix.triangular.Sparse(T, uplo, diag, order){
                .data = cpy.data,
                .idx = if (comptime order == .col_major) cpy.row else cpy.col,
                .ptr = ptr.ptr,
                .nnz = cpy.nnz,
                .rows = cpy.rows,
                .cols = cpy.cols,
                .flags = .{ .owns_data = true },
            };

            cpy.* = undefined;

            return result;
        }

        /// Cleans up the elements of the builder matrix, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.builder.Sparse(T, order)`):
        /// A pointer to the builder matrix to clean up.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Notes
        /// -----
        /// This function must be called before `deinit` if the elements are of
        /// arbitrary precision type to properly deinitialize them.
        pub fn cleanup(self: *Sparse(T, order), ctx: anytype) void {
            return _cleanup(self.data, self.nnz, ctx);
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary, in the specified order up to position (i, i), exclusive.
        fn _cleanup(data: [*]T, i: u32, ctx: anytype) void {
            switch (comptime types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    // No cleanup needed for fixed precision types.
                },
                .integer, .rational, .real, .complex => {
                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );

                    var _i: u32 = 0;
                    while (_i < i) : (_i += 1) {
                        ops.deinit(
                            &data[_i],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                },
            }
        }
    };
}
