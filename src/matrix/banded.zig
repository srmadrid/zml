const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../array.zig");

/// Banded matrix type, represented as a contiguous `n × (lower + upper + 1)`
/// array of elements of type `T`, stored in either column-major (`n = cols`) or
/// row-major (`n = rows`) order. For column-major storage, each element
/// `(i, j)` is stored at `(upper + i - j, j)`. For row-major storage, each
/// element `(i, j)` is stored at `(i, lower + j - i)`.
pub fn Banded(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("matrix.Banded requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        ld: u32, // leading dimension
        lower: u32,
        upper: u32,
        flags: Flags = .{},

        pub const empty: Banded(T, order) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
            .lower = 0,
            .upper = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new matrix with the specified rows, columns, and lower
        /// and upper bandwidths.
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
        /// `lower` (`u32`):
        /// The number of sub-diagonals.
        ///
        /// `upper` (`u32`):
        /// The number of super-diagonals.
        ///
        /// Returns
        /// -------
        /// `matrix.Banded(T, order)`:
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
        /// `matrix.Error.InvalidBandwidth`:
        /// If `lower` is greater than or equal to `rows`, or if `upper` is
        /// greater than or equal to `cols`.
        ///
        /// Notes
        /// -----
        /// The elements are not initialized.
        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
            lower: u32,
            upper: u32,
        ) !Banded(T, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            if (lower >= rows or upper >= cols)
                return matrix.Error.InvalidBandwidth;

            return .{
                .data = (try allocator.alloc(T, (lower + upper + 1) * (if (comptime order == .col_major) cols else rows))).ptr,
                .rows = rows,
                .cols = cols,
                .ld = lower + upper + 1,
                .lower = lower,
                .upper = upper,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new matrix with the specified rows, columns, and lower
        /// and upper bandwidths, filled with the specified value.
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
        /// `lower` (`u32`):
        /// The number of sub-diagonals.
        ///
        /// `upper` (`u32`):
        /// The number of super-diagonals.
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
        /// `matrix.Banded(T, order)`:
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
            lower: u32,
            upper: u32,
            value: anytype,
            ctx: anytype,
        ) !Banded(T, order) {
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

            var mat: Banded(T, order) = try .init(allocator, rows, cols, lower, upper);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                while (j < cols) : (j += 1) {
                    i = if (j < upper) 0 else j - upper;
                    while (i < int.min(rows, j + lower + 1)) : (i += 1) {
                        mat.data[mat._index(i, j)] = try ops.init(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );

                        try ops.set(
                            &mat.data[mat._index(i, j)],
                            value,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            } else {
                while (i < rows) : (i += 1) {
                    j = if (i < lower) 0 else i - lower;
                    while (j < int.min(cols, i + upper + 1)) : (j += 1) {
                        mat.data[mat._index(i, j)] = try ops.init(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );

                        try ops.set(
                            &mat.data[mat._index(i, j)],
                            value,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            }

            return mat;
        }

        /// Initializes a new identity matrix of the specified size, and lower
        /// and upper bandwidths.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `size` (`u32`):
        /// The size of the (square) matrix.
        ///
        /// `lower` (`u32`):
        /// The number of sub-diagonals.
        ///
        /// `upper` (`u32`):
        /// The number of super-diagonals.
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
        /// `matrix.Banded(T, order)`:
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
            lower: u32,
            upper: u32,
            ctx: anytype,
        ) !Banded(T, order) {
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

            var mat: Banded(T, order) = try .init(allocator, size, size, lower, upper);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                while (j < size) : (j += 1) {
                    i = if (j < upper) 0 else j - upper;
                    while (i < j) : (i += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    mat.data[mat._index(j, j)] = try constants.one(
                        T,
                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                    );

                    i += 1;

                    while (i < int.min(size, j + lower + 1)) : (i += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            } else {
                while (i < size) : (i += 1) {
                    j = if (i < lower) 0 else i - lower;
                    while (j < i) : (j += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    mat.data[mat._index(i, i)] = try constants.one(
                        T,
                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                    );

                    j += 1;

                    while (j < int.min(size, i + upper + 1)) : (j += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Banded(T, order)`):
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
        pub fn deinit(self: *Banded(T, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. (self.lower + self.upper + 1) * (if (comptime order == .col_major) self.cols else self.rows)]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Banded(T, order)`):
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
        pub fn get(self: *const Banded(T, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (c + self.lower < r or c > r + self.upper) {
                return constants.zero(T, .{}) catch unreachable;
            } else {
                return self.data[self._index(r, c)];
            }
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Banded(T, order)`):
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
        pub inline fn at(self: *const Banded(T, order), r: u32, c: u32) T {
            return self.data[self._index(r, c)];
        }

        /// Sets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Banded(T, order)`):
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
        /// If the position `(r, c)` is outside the banded structure.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn set(self: *Banded(T, order), r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (c + self.lower < r or c > r + self.upper) {
                return matrix.Error.BreaksStructure;
            } else {
                self.data[self._index(r, c)] = value;
            }
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Banded(T, order)`):
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
        pub inline fn put(self: *Banded(T, order), r: u32, c: u32, value: T) void {
            self.data[self._index(r, c)] = value;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Banded(T, order)`):
        /// The matrix to transpose.
        ///
        /// Returns
        /// -------
        /// `matrix.Banded(T, order.invert())`:
        /// The transposed matrix.
        pub fn transpose(self: Banded(T, order)) Banded(T, order.invert()) {
            return .{
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .ld = self.ld,
                .lower = self.upper,
                .upper = self.lower,
                .flags = .{ .owns_data = false },
            };
        }

        fn submatrix(
            self: *const Banded(T, order),
            start: u32,
            row_end: u32,
            col_end: u32,
        ) !Banded(T, order) {
            _ = self;
            _ = start;
            _ = row_end;
            _ = col_end;
            @compileError("Banded submatrix does not work!\n");
            // if (start >= int.min(self.rows, self.cols) or
            //     row_end > self.rows or col_end > self.cols or
            //     row_end < start or col_end < start)
            //     return matrix.Error.InvalidRange;

            // const sub_rows = row_end - start;
            // const sub_cols = col_end - start;

            // return .{
            //     .data = self.data + start * self.strides[0] + start * self.strides[1] -
            //         (if (self.flags.order == .col_major) self.upper else self.lower),
            //     .rows = sub_rows,
            //     .cols = sub_cols,
            //     .strides = self.strides,
            //     .lower = int.min(self.lower, sub_rows - 1),
            //     .upper = int.min(self.upper, sub_cols - 1),
            //     .flags = .{
            //         .order = self.flags.order,
            //         .owns_data = false,
            //     },
            // };
        }

        /// Copies the banded matrix to a general dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.Banded(T, order)`):
        /// A pointer to the matrix to copy.
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
            self: *const Banded(T, order),
            allocator: std.mem.Allocator,
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

            var mat: matrix.general.Dense(T, order) = try .init(allocator, self.rows, self.cols);
            errdefer mat.deinit(allocator);

            if (comptime order == .col_major) {
                var j: u32 = 0;
                while (j < self.cols) : (j += 1) {
                    var i: u32 = 0;
                    while (i < if (j < self.upper) 0 else j - self.upper) : (i += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    while (i < int.min(self.rows, j + self.lower + 1)) : (i += 1) {
                        mat.data[mat._index(i, j)] = try ops.copy(
                            self.data[self._index(i, j)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    while (i < self.rows) : (i += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            } else {
                var i: u32 = 0;
                while (i < self.rows) : (i += 1) {
                    var j: u32 = 0;
                    while (j < if (i < self.lower) 0 else i - self.lower) : (j += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    while (j < int.min(self.cols, i + self.upper + 1)) : (j += 1) {
                        mat.data[mat._index(i, j)] = try ops.copy(
                            self.data[self._index(i, j)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }

                    while (j < self.cols) : (j += 1) {
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
            self: *const Banded(T, order),
            allocator: std.mem.Allocator,
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
                        while (i < if (j < self.upper) 0 else j - self.upper) : (i += 1) {
                            result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                        }

                        while (i <= int.min(self.rows - 1, j + self.lower)) : (i += 1) {
                            result.data[i + j * result.strides[1]] = self.data[(self.upper + i - j) + j * self.ld];
                        }

                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < if (i < self.lower) 0 else i - self.lower) : (j += 1) {
                            result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                        }

                        while (j <= int.min(self.cols - 1, i + self.upper)) : (j += 1) {
                            result.data[i * result.strides[0] + j] = self.data[i * self.ld + (self.lower + j - i)];
                        }

                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub inline fn _index(self: *const Banded(T, order), r: u32, c: u32) u32 {
            return if (comptime order == .col_major)
                (self.upper + r - c) + c * self.ld
            else
                r * self.ld + (self.lower + c - r);
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.Banded(T, order)`):
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
        pub fn cleanup(self: *Banded(T, order), ctx: anytype) void {
            return self._cleanup(self.rows, self.cols, order, ctx);
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary, in the specified order up to position (i, j), exclusive.
        /// In other words, if iter_order is column-major, all elements from
        /// (0, 0) to (i - 1, j) are cleaned up, and if iter_order is row-major,
        /// all elements from (0, 0) to (i, j - 1) are cleaned up.
        pub fn _cleanup(self: *Banded(T, order), i: u32, j: u32, iter_order: Order, ctx: anytype) void {
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

                    if (iter_order == .col_major) {
                        var _j: u32 = 0;
                        while (_j <= int.min(j, self.cols - 1)) : (_j += 1) {
                            var _i: u32 = if (_j < self.upper) 0 else _j - self.upper;
                            if (_j == j) {
                                while (_i < int.min(int.min(i, self.rows), _j + self.lower + 1)) : (_i += 1) {
                                    ops.deinit(
                                        &self.data[self._index(_i, _j)],
                                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                    );
                                }
                            } else {
                                while (_i < int.min(self.rows, _j + self.lower + 1)) : (_i += 1) {
                                    ops.deinit(
                                        &self.data[self._index(_i, _j)],
                                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                    );
                                }
                            }
                        }
                    } else {
                        var _i: u32 = 0;
                        while (_i <= int.min(i, self.rows - 1)) : (_i += 1) {
                            var _j: u32 = if (_i < self.lower) 0 else _i - self.lower;
                            if (_i == i) {
                                while (_j < int.min(int.min(j, self.cols), _i + self.upper + 1)) : (_j += 1) {
                                    ops.deinit(
                                        &self.data[self._index(_i, _j)],
                                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                    );
                                }
                            } else {
                                while (_j < int.min(self.cols, _i + self.upper + 1)) : (_j += 1) {
                                    ops.deinit(
                                        &self.data[self._index(_i, _j)],
                                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                    );
                                }
                            }
                        }
                    }
                },
            }
        }
    };
}
