//! Storage scheme:
//!
//! A 1-d array of size `n` stores the permutation of 0..n-1. If `direction` is
//! `.forward`, the element at index `i` indicates the column index of the 1 in
//! row `i`, i.e., if `data[i] = j`, then the element at row `i` and column `j`
//! is 1, and all other elements in row `i` are 0. If `direction` is
//! `.backward`, the same applies but for columns, i.e., if `data[j] = i`,
//! then the element at row `i` and column `j` is 1, and all other elements in
//! column `j` are 0.

const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const Uplo = types.Uplo;
const Diag = types.Diag;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../array.zig");

pub const Direction = enum {
    forward,
    backward,
};

/// Permutation matrix type, represented as a contiguous array of `size`
/// elements of type `u32` holding a permutation of `0 .. size - 1`. If
/// `direction` is forward, the element at index `i` indicates the column
/// index of the 1 in row `i`, i.e., if `data[i] = j`, then the element at
/// row `i` and column `j` is 1, and all other elements in row `i` are 0. If
/// `direction` is backward, the same applies but for columns, i.e., if
/// `data[j] = i`, then the element at row `i` and column `j` is 1, and all
/// other elements in column `j` are 0.
pub fn Permutation(T: type) type {
    if (!types.isNumeric(T))
        @compileError("matrix.Permutation requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]u32,
        size: u32,
        direction: Direction = .forward,
        flags: Flags = .{},

        pub const empty: Permutation(T) = .{
            .data = &.{},
            .size = 0,
            .direction = .forward,
            .flags = .{ .owns_data = false },
        };

        pub fn tp() type {
            return T;
        }

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
        /// `matrix.Permutation(T)`:
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
        ) !Permutation(T) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(u32, size)).ptr,
                .size = size,
                .direction = .forward,
                .flags = .{ .owns_data = true },
            };
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
        /// Returns
        /// -------
        /// `matrix.Permutation(T)`:
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
        ) !Permutation(T) {
            var mat: Permutation(T) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            while (i < size) : (i += 1) {
                mat.data[i] = i;
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Permutation(T)`):
        /// A pointer to the matrix to deinitialize.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory deallocation. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// Returns
        /// -------
        /// `void`
        pub fn deinit(self: *Permutation(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.size]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Permutation(T)`):
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
        pub fn get(self: *const Permutation(T), r: u32, c: u32) !T {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            if (self.direction == .forward) {
                if (self.data[r] == c) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            } else {
                if (self.data[c] == r) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            }
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Permutation(T)`):
        /// A pointer to the matrix to get the element from.
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
        pub inline fn at(self: *const Permutation(T), row: u32, col: u32) T {
            if (self.direction == .forward) {
                if (self.data[row] == col) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            } else {
                if (self.data[col] == row) {
                    return constants.one(T, .{}) catch unreachable;
                } else {
                    return constants.zero(T, .{}) catch unreachable;
                }
            }
        }

        // pub fn set(self: *Permutation(T), row: u32, col: u32, value: u32) !void {
        //     if (row >= self.size or col >= self.size)
        //         return matrix.Error.PositionOutOfBounds;

        //     if (value != 0 and value != 1)
        //         return matrix.Error.BreaksStructure;
        // }

        // pub inline fn put(self: *Permutation(T), row: u32, col: u32, value: u32) void {
        //     // Unchecked version of set. Assumes row and col are valid and
        //     // in banded range.
        //     if (value == 1) {
        //         self.data[row] = col;
        //     }
        // }

        /// Returns a transposed view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Permutation(T)`):
        /// The matrix to transpose.
        ///
        /// Returns
        /// -------
        /// `matrix.Permutation(T)`:
        /// The transposed matrix.
        pub fn transpose(self: Permutation(T)) Permutation(T) {
            return .{
                .data = self.data,
                .size = self.size,
                .direction = if (self.direction == .forward) .backward else .forward,
                .flags = self.flags,
            };
        }

        /// Copies the permutation matrix to a general dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Permutation(T)`):
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
            self: Permutation(T),
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
                if (self.direction == .forward) {
                    while (j < self.size) : (j += 1) {
                        i = 0;
                        while (i < self.size) : (i += 1) {
                            mat.data[mat._index(i, j)] = if (self.data[i] == j)
                                try constants.one(
                                    T,
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                )
                            else
                                try constants.zero(
                                    T,
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                );
                        }
                    }
                } else {
                    while (i < self.size) : (i += 1) {
                        j = 0;
                        while (j < self.size) : (j += 1) {
                            mat.data[mat._index(i, j)] = if (self.data[j] == i)
                                try constants.one(
                                    T,
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                )
                            else
                                try constants.zero(
                                    T,
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                );
                        }
                    }
                }
            } else {
                if (self.direction == .forward) {
                    while (i < self.size) : (i += 1) {
                        j = 0;
                        while (j < self.size) : (j += 1) {
                            mat.data[mat._index(i, j)] = if (self.data[i] == j)
                                try constants.one(
                                    T,
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                )
                            else
                                try constants.zero(
                                    T,
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                );
                        }
                    }
                } else {
                    while (j < self.size) : (j += 1) {
                        i = 0;
                        while (i < self.size) : (i += 1) {
                            mat.data[mat._index(i, j)] = if (self.data[j] == i)
                                try constants.one(
                                    T,
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                )
                            else
                                try constants.zero(
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
            self: *const Permutation(T),
            allocator: std.mem.Allocator,
            comptime order: Order,
            ctx: anytype,
        ) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
            errdefer result.deinit(allocator);

            if (comptime order == .col_major) {
                if (self.direction == .forward) {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i + j * result.strides[0]] = if (self.data[i] == j)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[i + j * result.strides[0]] = if (self.data[j] == i)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                if (self.direction == .forward) {
                    var i: u32 = 0;
                    while (i < self.size) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[i * result.strides[0] + j] = if (self.data[i] == j)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var j: u32 = 0;
                    while (j < self.size) : (j += 1) {
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i * result.strides[0] + j] = if (self.data[j] == i)
                                constants.one(T, ctx) catch unreachable
                            else
                                constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            }

            return result;
        }

        // pub fn submatrix(
        //     self: *const Permutation(T),
        //     start: u32,
        //     end: u32,
        // ) !? {
        //     if (start >= self.size or end > self.size or start >= end)
        //         return matrix.Error.InvalidRange;

        //     const sub_size = end - start;

        //     return .{
        //         .data = self.data,
        //         .size = sub_size,
        //         .osize = self.osize,
        //         .offset = self.offset + start,
        //         .sdoffset = self.sdoffset,
        //         .flags = .{
        //             .owns_data = false,
        //         },
        //     };
        // }
    };
}
