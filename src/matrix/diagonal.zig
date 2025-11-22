const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../array.zig");

/// Diagonal matrix type, represented as a contiguous array of `min(rows, cols)`
/// elements of type `T`.
pub fn Diagonal(T: type) type {
    if (!types.isNumeric(T))
        @compileError("matrix.Diagonal requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        flags: Flags = .{},

        pub const empty: Diagonal(T) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new matrix with the specified rows and columns.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `rows` (`u32`):
        /// The rows of the matrix.
        ///
        /// `cols` (`u32`):
        /// The columns of the matrix.
        ///
        /// Returns
        /// -------
        /// `matrix.Diagonal(T)`:
        /// The newly initialized matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.ZeroDimension`:
        /// If either `rows` or `cols` is zero.
        ///
        /// Notes
        /// -----
        /// The elements are not initialized.
        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
        ) !Diagonal(T) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return Diagonal(T){
                .data = (try allocator.alloc(T, int.min(rows, cols))).ptr,
                .rows = rows,
                .cols = cols,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new matrix with the specified rows and columns, with
        /// the diagonal part of the matrix filled with the specified value.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `rows` (`u32`):
        /// The rows of the matrix.
        ///
        /// `cols` (`u32`):
        /// The columns of the matrix.
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
        /// `matrix.Diagonal(T)`:
        /// The newly initialized matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.ZeroDimension`:
        /// If either `rows` or `cols` is zero.
        ///
        /// Notes
        /// -----
        /// The matrix does not take ownership of `value` if it is an arbitrary
        /// precision type.
        pub fn full(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            value: anytype,
            ctx: anytype,
        ) !Diagonal(T) {
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

            var mat: Diagonal(T) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;

            errdefer mat._cleanup(i, ctx);

            while (i < int.min(rows, cols)) : (i += 1) {
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
        /// `matrix.Diagonal(T)`:
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
        ) !Diagonal(T) {
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

            var mat: Diagonal(T) = try .init(allocator, size, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;

            errdefer mat._cleanup(i, ctx);

            while (i < size) : (i += 1) {
                mat.data[i] = try constants.one(
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
        /// `self` (`*matrix.Diagonal(T)`):
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
        pub fn deinit(self: *Diagonal(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..(int.min(self.rows, self.cols))]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Diagonal(T)`):
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
        pub fn get(self: *const Diagonal(T), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (r != c)
                return constants.zero(T, .{}) catch unreachable;

            return self.data[r];
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Diagonal(T)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be within bounds and
        /// equal to `c`.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be within bounds
        /// and equal to `r`.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub inline fn at(self: *const Diagonal(T), r: u32, c: u32) T {
            _ = c;
            return self.data[r];
        }

        /// Sets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Diagonal(T)`):
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
        /// If `r` is not equal to `c`.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn set(self: *Diagonal(T), r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (r != c)
                return matrix.Error.BreaksStructure;

            self.data[r] = value;
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Diagonal(T)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be within bounds and
        /// equal to `c`.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be within bounds
        /// and equal to `r`.
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
        pub inline fn put(self: *Diagonal(T), r: u32, c: u32, value: T) void {
            _ = c;
            self.data[r] = value;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Diagonal(T)`):
        /// The matrix to transpose.
        ///
        /// Returns
        /// -------
        /// `matrix.Diagonal(T)`:
        /// The transposed matrix.
        pub fn transpose(self: Diagonal(T)) Diagonal(T) {
            return Diagonal(T){
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Diagonal(T)`):
        /// The matrix to get the submatrix from.
        ///
        /// `start` (`u32`):
        /// The starting diagonal index of the submatrix (inclusive).
        ///
        /// `row_end` (`u32`):
        /// The ending row index of the submatrix (exclusive). Must be greater
        /// than `start`.
        ///
        /// `col_end` (`u32`):
        /// The ending column index of the submatrix (exclusive). Must be
        /// greater than `start`.
        ///
        /// Returns
        /// -------
        /// `matrix.Diagonal(T)`:
        /// The submatrix.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.InvalidRange`:
        /// If the specified range is invalid.
        pub fn submatrix(
            self: *const Diagonal(T),
            start: u32,
            row_end: u32,
            col_end: u32,
        ) !Diagonal(T) {
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - start;
            const sub_cols = col_end - start;

            return Diagonal(T){
                .data = self.data + start,
                .rows = sub_rows,
                .cols = sub_cols,
                .flags = .{ .owns_data = false },
            };
        }

        /// Copies the symmetric matrix to a general dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Diagonal(T)`):
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
            self: Diagonal(T),
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
                while (j < mat.cols) : (j += 1) {
                    i = 0;
                    while (i < int.min(j, mat.rows)) : (i += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    mat.data[mat._index(j, j)] = try ops.copy(
                        self.data[j],
                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                    );

                    i = j + 1;
                    while (i < mat.rows) : (i += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            } else {
                while (i < mat.rows) : (i += 1) {
                    j = 0;
                    while (j < int.min(i, mat.cols)) : (j += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    mat.data[mat._index(i, i)] = try ops.copy(
                        self.data[i],
                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                    );

                    j = i + 1;
                    while (j < mat.cols) : (j += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            }

            return mat;
        }

        pub fn copyToDenseArray(
            self: *const Diagonal(T),
            allocator: std.mem.Allocator,
            comptime order: Order,
            ctx: anytype,
        ) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.rows, self.cols });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    var j: u32 = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: u32 = 0;
                        while (i < int.min(j, self.rows)) : (i += 1) {
                            result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (j < int.min(self.rows, self.cols)) {
                            result.data[j * result.ld + j] = self.data[j];
                        }

                        i = j + 1;
                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.ld] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < int.min(i, self.cols)) : (j += 1) {
                            result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        if (i < int.min(self.rows, self.cols)) {
                            result.data[i * result.ld + i] = self.data[i];
                        }

                        j = i + 1;
                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.ld + j] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Diagonal(T, order)`):
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
        pub fn cleanup(self: *Diagonal(T), ctx: anytype) void {
            return self._cleanup(int.min(self.rows, self.cols), ctx);
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary, in the specified order up to position (i, i), exclusive.
        pub fn _cleanup(self: *Diagonal(T), i: u32, ctx: anytype) void {
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
                    while (_i < int.min(i, int.min(self.rows, self.cols))) : (_i += 1) {
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
