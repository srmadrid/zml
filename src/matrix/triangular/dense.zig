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

const array = @import("../../array.zig");

/// Dense triangular matrix type, represented as a contiguous array of
/// `min(rows, cols) × cols` or `rows × min(rows, cols)` elements of type `T`,
/// depending on `uplo`, stored in either column-major or row-major order with a
/// specified leading dimension. Only the upper or lower triangular part of the
/// matrix is accessed, depending on the `uplo` parameter, and the diagonal can
/// be either unit, meaning all diagonal elements are assumed to be 1 and not
/// accessed, or non-unit, meaning the diagonal elements are accessed normally.
pub fn Dense(T: type, uplo: Uplo, diag: Diag, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("matrix.triangular.Dense requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        ld: u32, // leading dimension
        flags: Flags = .{},

        pub const empty: Dense(T, uplo, diag, order) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
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
        /// `matrix.triangular.Dense(T, uplo, diag, order)`:
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
        ) !Dense(T, uplo, diag, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            const m: u32 = if (comptime uplo == .upper)
                int.min(rows, cols)
            else
                rows;
            const n: u32 = if (comptime uplo == .upper)
                cols
            else
                int.min(rows, cols);

            return .{
                .data = (try allocator.alloc(T, m * n)).ptr,
                .rows = rows,
                .cols = cols,
                .ld = if (comptime order == .col_major) m else n,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new matrix with the specified rows and columns, with
        /// the triangular part of the matrix filled with the specified value.
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
        /// `matrix.triangular.Dense(T, uplo, diag, order)`:
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
        ) !Dense(T, uplo, diag, order) {
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

            var mat: Dense(T, uplo, diag, order) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        while (j < cols) : (j += 1) {
                            i = 0;
                            while (i < int.min(j, rows)) : (i += 1) {
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
                    } else { // cun
                        while (j < cols) : (j += 1) {
                            i = 0;
                            while (i < int.min(j + 1, rows)) : (i += 1) {
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
                } else {
                    if (comptime diag == .unit) { // clu
                        while (j < int.min(rows, cols)) : (j += 1) {
                            i = j + 1;
                            while (i < rows) : (i += 1) {
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
                    } else { // cln
                        while (j < int.min(rows, cols)) : (j += 1) {
                            i = j;
                            while (i < rows) : (i += 1) {
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
                }
            } else {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // ruu
                        while (i < int.min(rows, cols)) : (i += 1) {
                            j = i + 1;
                            while (j < cols) : (j += 1) {
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
                    } else { // run
                        while (i < int.min(rows, cols)) : (i += 1) {
                            j = i;
                            while (j < cols) : (j += 1) {
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
                } else {
                    if (comptime diag == .unit) { // rlu
                        while (i < rows) : (i += 1) {
                            j = 0;
                            while (j < int.min(i, cols)) : (j += 1) {
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
                    } else { // rln
                        while (i < rows) : (i += 1) {
                            j = 0;
                            while (j < int.min(i + 1, cols)) : (j += 1) {
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
                }
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
        /// `matrix.triangular.Dense(T, uplo, diag, order)`:
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
        ) !Dense(T, uplo, diag, order) {
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

            var mat: Dense(T, uplo, diag, order) = try .init(allocator, size, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                if (comptime uplo == .upper) { // cu
                    while (j < size) : (j += 1) {
                        i = 0;
                        while (i < j) : (i += 1) {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        if (comptime diag == .non_unit) {
                            mat.data[mat._index(j, j)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // cl
                    while (j < size) : (j += 1) {
                        if (comptime diag == .non_unit) {
                            mat.data[mat._index(j, j)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        i = j + 1;

                        while (i < size) : (i += 1) {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    while (i < size) : (i += 1) {
                        if (comptime diag == .non_unit) {
                            mat.data[mat._index(i, i)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        j = i + 1;

                        while (j < size) : (j += 1) {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // rl
                    while (i < size) : (i += 1) {
                        j = 0;
                        while (j < i) : (j += 1) {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        if (comptime diag == .non_unit) {
                            mat.data[mat._index(i, i)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
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
        /// `self` (`*matrix.triangular.Dense(T, uplo, diag, order)`):
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
        pub fn deinit(self: *Dense(T, uplo, diag, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..if (comptime uplo == .upper)
                    int.min(self.rows, self.cols) * self.cols
                else
                    self.rows * int.min(self.rows, self.cols)]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.triangular.Dense(T, uplo, diag, order)`):
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
        pub fn get(self: *const Dense(T, uplo, diag, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (r > c)
                    return constants.zero(T, .{}) catch unreachable;
            } else {
                if (r < c)
                    return constants.zero(T, .{}) catch unreachable;
            }

            if (comptime diag == .unit) {
                if (r == c)
                    return constants.one(T, .{}) catch unreachable;
            }

            return self.data[self._index(r, c)];
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.triangular.Dense(T, uplo, diag, order)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be within bounds, on
        /// the correct triangular part, and outside the diagonal if `diag` is
        /// unit.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be within bounds,
        /// on the correct triangular part, and outside the diagonal if `diag`
        /// is unit.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub inline fn at(self: *const Dense(T, uplo, diag, order), r: u32, c: u32) T {
            return self.data[self._index(r, c)];
        }

        /// Sets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.triangular.Dense(T, uplo, diag, order)`):
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
        /// If `r == c` and `diag` is unit, or if the position `(r, c)` is
        /// outside the correct triangular part of the matrix.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn set(self: *Dense(T, uplo, diag, order), row: u32, col: u32, value: T) !void {
            if (row >= self.rows or col >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (row > col)
                    return matrix.Error.BreaksStructure;
            } else {
                if (row < col)
                    return matrix.Error.BreaksStructure;
            }

            if (comptime diag == .unit) {
                if (row == col)
                    return matrix.Error.BreaksStructure;
            }

            self.data[self._index(row, col)] = value;
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.triangular.Dense(T, uplo, diag, order)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be within bounds, on
        /// the correct triangular part, and outside the diagonal if `diag` is
        /// unit.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be within bounds,
        /// on the correct triangular part, and outside the diagonal if `diag`
        /// is unit.
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
        pub inline fn put(self: *Dense(T, uplo, diag, order), row: u32, col: u32, value: T) void {
            self.data[self._index(row, col)] = value;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.triangular.Dense(T, uplo, diag, order)`):
        /// The matrix to transpose.
        ///
        /// Returns
        /// -------
        /// `matrix.triangular.Dense(T, uplo.invert(), diag, order.invert())`:
        /// The transposed matrix.
        pub fn transpose(self: Dense(T, uplo, diag, order)) Dense(T, uplo.invert(), diag, order.invert()) {
            return .{
                .data = self.data,
                .rows = self.cols,
                .cols = self.rows,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.triangular.Dense(T, uplo, diag, order)`):
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
        /// `matrix.triangular.Dense(T, uplo, diag, order)`:
        /// The submatrix.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.InvalidRange`:
        /// If the specified range is invalid.
        pub fn submatrix(
            self: *const Dense(T, uplo, diag, order),
            start: u32,
            row_end: u32,
            col_end: u32,
        ) !Dense(T, uplo, diag, order) {
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - start;
            const sub_cols = col_end - start;

            return .{
                .data = self.data + (start + start * self.ld),
                .rows = sub_rows,
                .cols = sub_cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Copies the triangular matrix to a general dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.triangular.Dense(T, uplo, diag, order)`):
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
            self: Dense(T, uplo, diag, order),
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

            var mat: matrix.general.Dense(T, order) = try .init(allocator, self.size, self.size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                if (comptime uplo == .upper) { // cu
                    while (j < mat.cols) : (j += 1) {
                        i = 0;
                        while (i < j) : (i += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        if (comptime diag == .unit) {
                            mat.data[mat._index(j, j)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else {
                            mat.data[mat._index(j, j)] = try ops.copy(
                                self.data[self._index(j, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        i = j + 1;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // cl
                    while (j < mat.cols) : (j += 1) {
                        i = 0;
                        while (i < j) : (i += 1) {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        if (comptime diag == .unit) {
                            mat.data[mat._index(j, j)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else {
                            mat.data[mat._index(j, j)] = try ops.copy(
                                self.data[self._index(j, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        i = j + 1;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    while (i < mat.rows) : (i += 1) {
                        j = 0;
                        while (j < i) : (j += 1) {
                            mat.data[mat._index(i, j)] = try constants.zero(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        if (comptime diag == .unit) {
                            mat.data[mat._index(i, i)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else {
                            mat.data[mat._index(i, i)] = try ops.copy(
                                self.data[self._index(i, i)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        j = i + 1;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // rl
                    while (i < mat.rows) : (i += 1) {
                        j = 0;
                        while (j < i) : (j += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        if (comptime diag == .unit) {
                            mat.data[mat._index(i, i)] = try constants.one(
                                T,
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        } else {
                            mat.data[mat._index(i, i)] = try ops.copy(
                                self.data[self._index(i, i)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        j = i + 1;
                        while (j < mat.cols) : (j += 1) {
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
            self: *const Dense(T, uplo, diag, order),
            allocator: std.mem.Allocator,
            ctx: anytype,
        ) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.rows, self.cols });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, self.rows)) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                            }

                            if (j < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
                                }
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: u32 = 0;
                            while (i < int.min(j, self.rows)) : (i += 1) {
                                result.data[i + j * result.strides[1]] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (j < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[j + j * result.strides[1]] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
                                }
                            }

                            i = j + 1;
                            while (i < self.rows) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, self.cols)) : (j += 1) {
                                result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                            }

                            if (i < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];
                                }
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: u32 = 0;
                            while (j < int.min(i, self.cols)) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                            }

                            if (i < int.min(self.rows, self.cols)) {
                                if (comptime diag == .unit) {
                                    result.data[i * result.strides[0] + i] = constants.one(T, ctx) catch unreachable;
                                } else {
                                    result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];
                                }
                            }

                            j = i + 1;
                            while (j < self.cols) : (j += 1) {
                                result.data[i * result.strides[0] + j] = constants.zero(T, ctx) catch unreachable;
                            }
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub inline fn _index(self: *const Dense(T, uplo, diag, order), r: u32, c: u32) u32 {
            return if (comptime order == .col_major)
                r + c * self.ld
            else
                r * self.ld + c;
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.triangular.Dense(T, uplo, diag, order)`):
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
        pub fn cleanup(self: *Dense(T, uplo, diag, order), ctx: anytype) void {
            return self._cleanup(self.rows, self.cols, order, ctx);
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary, in the specified order up to position (i, j), exclusive.
        /// In other words, if iter_order is column-major, all elements from
        /// (0, 0) to (i - 1, j) are cleaned up, and if iter_order is row-major,
        /// all elements from (0, 0) to (i, j - 1) are cleaned up.
        pub fn _cleanup(self: *Dense(T, uplo, diag, order), i: u32, j: u32, iter_order: Order, ctx: anytype) void {
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
                        if (comptime uplo == .upper) {
                            var _j: u32 = 0;
                            while (_j <= int.min(j, self.cols - 1)) : (_j += 1) {
                                var _i: u32 = 0;
                                if (_j == j) {
                                    const l: u32 = if (comptime diag == .unit) _j else _j + 1;
                                    while (_i < int.min(i, l)) : (_i += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    const l: u32 = if (comptime diag == .unit) _j else _j + 1;
                                    while (_i < l) : (_i += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                }
                            }
                        } else {
                            var _j: u32 = 0;
                            while (_j <= int.min(j, self.cols - 1)) : (_j += 1) {
                                var _i: u32 = if (comptime diag == .unit) _j + 1 else _j;
                                if (_j == j) {
                                    while (_i < int.min(i, self.rows)) : (_i += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    while (_i < self.rows) : (_i += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime uplo == .upper) {
                            var _i: u32 = 0;
                            while (_i <= int.min(i, self.rows - 1)) : (_i += 1) {
                                var _j: u32 = if (comptime diag == .unit) _i + 1 else _i;
                                if (_i == i) {
                                    while (_j < int.min(j, self.cols)) : (_j += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    while (_j < self.cols) : (_j += 1) {
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
                                var _j: u32 = 0;
                                if (_i == i) {
                                    const l: u32 = if (comptime diag == .unit) _i else _i + 1;
                                    while (_j < int.min(j, l)) : (_j += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    const l: u32 = if (comptime diag == .unit) _i else _i + 1;
                                    while (_j < l) : (_j += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                }
                            }
                        }
                    }
                },
            }
        }
    };
}
