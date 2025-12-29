const std = @import("std");

const types = @import("../../types.zig");
const Order = types.Order;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

/// Tridiagonal matrix type, represented as a contiguous array of `3 × size - 2`
/// elements of type `T`. This array contains the subdiagonal, diagonal, and
/// superdiagonal elements of the matrix, in that order.
pub fn Tridiagonal(T: type) type {
    if (!types.isNumeric(T))
        @compileError("matrix.Tridiagonal requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        osize: u32,
        offset: u32,
        sdoffset: u32,
        flags: Flags = .{},

        /// Type signatures
        pub const is_matrix = {};
        pub const is_tridiagonal = {};
        pub const is_general = {};

        /// Numeric type
        pub const Numeric = T;

        pub const empty: Tridiagonal(T) = .{
            .data = &.{},
            .size = 0,
            .osize = 0,
            .offset = 0,
            .sdoffset = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new matrix with the specified size.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `size` (`u32`):
        /// The size of the (square) matrix.
        ///
        /// Returns
        /// -------
        /// `matrix.Tridiagonal(T)`:
        /// The newly initialized matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.ZeroDimension`:
        /// If `size` is zero.
        ///
        /// Notes
        /// -----
        /// The elements are not initialized.
        pub fn init(
            allocator: std.mem.Allocator,
            size: u32,
        ) !Tridiagonal(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(T, 3 * size - 2)).ptr,
                .size = size,
                .osize = size,
                .offset = 0,
                .sdoffset = (size - 1) + size,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new matrix with the specified rows and columns, with
        /// the tridiagonal part filled with the specified value.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `size` (`u32`):
        /// The size of the (square) matrix.
        ///
        /// `value` (`anytype`):
        /// The value to fill the matrix with.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T` and the
        /// type of `value`. If  the context is missing required fields or
        /// contains unnecessary or wrongly typed fields, the compiler will emit
        /// a detailed error message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `matrix.Tridiagonal(T)`:
        /// The newly initialized matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.ZeroDimension`:
        /// If `size` is zero.
        ///
        /// Notes
        /// -----
        /// The matrix does not take ownership of `value` if it is an arbitrary
        /// precision type.
        pub fn full(
            allocator: std.mem.Allocator,
            size: u32,
            value: anytype,
            ctx: anytype,
        ) !Tridiagonal(T) {
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

            var mat: Tridiagonal(T) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;

            errdefer mat._cleanup(i, ctx);

            while (i < 3 * size - 2) : (i += 1) {
                mat.data[i] = try ops.init(
                    T,
                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                );

                try ops.set(
                    &mat.data[i],
                    value,
                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                );
            }

            return mat;
        }

        /// Initializes a new identity matrix of the specified size.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `size` (`u32`):
        /// The size of the (square) matrix.
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
        /// `matrix.Tridiagonal(T)`:
        /// The newly initialized identity matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.ZeroDimension`:
        /// If `size` is zero.
        pub fn eye(
            allocator: std.mem.Allocator,
            size: u32,
            ctx: anytype,
        ) !Tridiagonal(T) {
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

            if (size == 0)
                return matrix.Error.ZeroDimension;

            var mat: Tridiagonal(T) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;

            errdefer mat._cleanup(i, ctx);

            while (i < size - 1) : (i += 1) {
                mat.data[i] = try constants.zero(
                    T,
                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                );
            }

            while (i < 2 * size - 1) : (i += 1) {
                mat.data[i] = try constants.one(
                    T,
                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                );
            }

            while (i < 3 * size - 2) : (i += 1) {
                mat.data[i] = try constants.zero(
                    T,
                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                );
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Tridiagonal(T)`):
        /// A pointer to the matrix to deinitialize.
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
        pub fn deinit(self: *Tridiagonal(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..(3 * self.size - 2)]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Tridiagonal(T)`):
        /// A pointer to the matrix to get the element from.
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
        pub fn get(self: *const Tridiagonal(T), r: u32, c: u32) !T {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            return if (self._index(r, c)) |idx|
                self.data[idx]
            else
                constants.zero(T, .{}) catch unreachable;
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Tridiagonal(T)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be within bounds and
        /// in banded range.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be within bounds
        /// and in banded range.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub inline fn at(self: *const Tridiagonal(T), r: u32, c: u32) T {
            return self.data[self._index(r, c).?];
        }

        /// Sets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Tridiagonal(T)`):
        /// A pointer to the matrix to set the element in.
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
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `r` or `c` is out of bounds.
        ///
        /// `matrix.Error.BreaksStructure`:
        /// If the position `(r, c)` is outside the tridiagonal structure.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn set(self: *Tridiagonal(T), r: u32, c: u32, value: T) !void {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            if (self._index(r, c)) |idx| {
                self.data[idx] = value;
            } else {
                return matrix.Error.BreaksStructure;
            }
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Tridiagonal(T)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be within bounds and
        /// in banded range.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be within bounds
        /// and in banded range.
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
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub inline fn put(self: *Tridiagonal(T), r: u32, c: u32, value: T) void {
            self.data[self._index(r, c).?] = value;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Tridiagonal(T)`):
        /// The matrix to transpose.
        ///
        /// Returns
        /// -------
        /// `matrix.Tridiagonal(T)`:
        /// The transposed matrix.
        pub fn transpose(self: *const Tridiagonal(T)) Tridiagonal(T) {
            return .{
                .data = self.data,
                .size = self.size,
                .osize = self.osize,
                .offset = self.offset,
                .sdoffset = self.size + (self.size - 1) - self.sdoffset,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Tridiafonal(T)`):
        /// The matrix to get the submatrix from.
        ///
        /// `start` (`u32`):
        /// The starting diagonal index of the submatrix (inclusive).
        ///
        /// `end` (`u32`):
        /// The ending diagonal index of the submatrix (exclusive). Must be
        /// greater than `start`.
        ///
        /// Returns
        /// -------
        /// `matrix.Tridiagonal(T)`:
        /// The submatrix.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.InvalidRange`:
        /// If the specified range is invalid.
        pub fn submatrix(
            self: *const Tridiagonal(T),
            start: u32,
            end: u32,
        ) !Tridiagonal(T) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size: u32 = end - start;

            return .{
                .data = self.data,
                .size = sub_size,
                .osize = self.osize,
                .offset = self.offset + start,
                .sdoffset = self.sdoffset,
                .flags = .{ .owns_data = false },
            };
        }

        /// Copies the tridiagonal matrix to a general dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Tridiagonal(T)`):
        /// A pointer to the matrix to copy.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `order` (`Order`):
        /// The storage order of the resulting matrix.
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
        /// `matrix.general.Dense(T, order)`:
        /// The copied matrix.
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
        pub fn copyToGeneralDenseMatrix(
            self: Tridiagonal(T),
            allocator: std.mem.Allocator,
            comptime order: Order,
            ctx: anytype,
        ) !matrix.general.Dense(T, order) {
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

            var mat: matrix.general.Dense(T, order) = try .init(allocator, self.size, self.size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                while (j < self.size) : (j += 1) {
                    i = 0;
                    while (i < self.size) : (i += 1) {
                        if (i == j) { // Diagonal
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self.offset + j + (self.osize - 1)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else if (i == j + 1) { // Subdiagonal
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else if (i + 1 == j) { // Superdiagonal
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self.offset + j + self.sdoffset - 1],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            } else {
                while (i < self.size) : (i += 1) {
                    j = 0;
                    while (j < self.size) : (j += 1) {
                        if (i == j) { // Diagonal
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self.offset + j + (self.osize - 1)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else if (i == j + 1) { // Subdiagonal
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else if (i + 1 == j) { // Superdiagonal
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self.offset + j + self.sdoffset - 1],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            }

            return mat;
        }

        pub fn copyToDenseArray(
            self: *const Tridiagonal(T),
            allocator: std.mem.Allocator,
            comptime order: Order,
            ctx: anytype,
        ) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            if (i == j) { // Diagonal
                                result.data[j + j * result.strides[1]] = self.data[self.offset + j + (self.osize - 1)];
                            } else if (i == j + 1) { // Subdiagonal
                                result.data[i + j * result.strides[1]] = self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1];
                            } else if (i + 1 == j) { // Superdiagonal
                                result.data[i + j * result.strides[1]] = self.data[self.offset + j + self.sdoffset - 1];
                            } else {
                                result.data[i + j * result.strides[1]] = constants.zero(T, .{}) catch unreachable;
                            }
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            if (i == j) { // Diagonal
                                result.data[i * result.strides[0] + i] = self.data[self.offset + i + (self.osize - 1)];
                            } else if (i == j + 1) { // Subdiagonal
                                result.data[i * result.strides[0] + j] = self.data[self.offset + i + self.osize + (self.osize - 1) - self.sdoffset - 1];
                            } else if (i + 1 == j) { // Superdiagonal
                                result.data[i * result.strides[0] + j] = self.data[self.offset + j + self.sdoffset - 1];
                            } else {
                                result.data[i * result.strides[0] + j] = constants.zero(T, .{}) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub inline fn _index(self: *const Tridiagonal(T), r: u32, c: u32) ?u32 {
            const diff: i32 = types.scast(i32, c) - types.scast(i32, r);

            return switch (diff) {
                -1 => self.offset + r + self.osize + (self.osize - 1) - self.sdoffset - 1, // Subdiagonal
                0 => self.offset + r + (self.osize - 1), // Diagonal
                1 => self.offset + c + self.sdoffset - 1, // Superdiagonal
                else => null, // Trash value, should not happen
            };
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Tridiagonal(T)`):
        /// A pointer to the matrix to clean up.
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
        pub fn cleanup(self: *Tridiagonal(T), ctx: anytype) void {
            return self._cleanup(3 * self.size - 2, ctx);
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary, sequentially, up to the specified index.
        pub fn _cleanup(self: *Tridiagonal(T), i: u32, ctx: anytype) void {
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
                            &self.data[_i],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                },
            }
        }
    };
}
