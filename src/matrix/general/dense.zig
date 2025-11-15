const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const Diag = types.Diag;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const vector = @import("../../vector.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

const linalg = @import("../../linalg.zig");

/// Dense general matrix type, represented as a contiguous array of elements of
/// type `T`, stored in either column-major or row-major order.
pub fn Dense(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("Dense requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        ld: u32, // leading dimension
        flags: Flags = .{},

        pub const empty: Dense(T, order) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.general.Dense(T, order)` with the
        /// specified rows and columns.
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
        /// `matrix.general.Dense(T, order)`:
        /// The newly initialized matrix.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `matrix.Error.ZeroDimension`:
        /// If either `rows` or `cols` is zero.
        pub fn init(
            allocator: std.mem.Allocator,
            rows: u32,
            cols: u32,
        ) !Dense(T, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(T, rows * cols)).ptr,
                .rows = rows,
                .cols = cols,
                .ld = if (comptime order == .col_major) rows else cols,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new `matrix.general.Dense(T, order)` with the
        /// specified rows and columns, filled with the specified value.
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
        /// `matrix.general.Dense(T, order)`:
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
        ) !Dense(T, order) {
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

            var mat: Dense(T, order) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(mat._index(i, j), order, ctx);

            if (comptime order == .col_major) {
                while (j < cols) : (j += 1) {
                    i = 0;
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
            } else {
                while (i < rows) : (i += 1) {
                    j = 0;
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

            return mat;
        }

        /// Initializes a new identity `matrix.general.Dense(T, order)` of the
        /// specified size.
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
        /// `matrix.general.Dense(T, order)`:
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
        ) !Dense(T, order) {
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

            var mat: Dense(T, order) = try .init(allocator, size, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(mat._index(i, j), order, ctx);

            if (comptime order == .col_major) {
                while (j < size) : (j += 1) {
                    i = 0;
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

                    while (i < size) : (i += 1) {
                        mat.data[mat._index(i, j)] = try constants.zero(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            } else {
                while (i < size) : (i += 1) {
                    j = 0;
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

                    while (j < size) : (j += 1) {
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
        /// `self` (`*matrix.general.Dense(T, order)`):
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
        pub fn deinit(self: *Dense(T, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
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
        pub fn get(self: *const Dense(T, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            return self.data[self._index(r, c)];
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be valid.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be valid.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub inline fn at(self: *const Dense(T, order), r: u32, c: u32) T {
            return self.data[self._index(r, c)];
        }

        /// Sets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.general.Dense(T, order)`):
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
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn set(self: *Dense(T, order), r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            self.data[self._index(r, c)] = value;
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.general.Dense(T, order)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be valid.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be valid.
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
        pub inline fn put(self: *Dense(T, order), r: u32, c: u32, value: T) void {
            self.data[self._index(r, c)] = value;
        }

        /// Creates a copy of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
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
        pub fn copy(self: *const Dense(T, order), allocator: std.mem.Allocator, ctx: anytype) !Dense(T, order) {
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

            var mat: Dense(T, order) = try .init(allocator, self.rows, self.cols);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(mat._index(i, j), order, ctx);

            if (comptime order == .col_major) {
                while (j < mat.cols) : (j += 1) {
                    i = 0;
                    while (i < mat.rows) : (i += 1) {
                        mat.data[mat._index(i, j)] = try ops.copy(
                            self.data[self._index(i, j)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            } else {
                while (i < mat.rows) : (i += 1) {
                    j = 0;
                    while (j < mat.cols) : (j += 1) {
                        mat.data[mat._index(i, j)] = try ops.copy(
                            self.data[self._index(i, j)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            }

            return mat;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// The matrix to transpose.
        ///
        /// Returns
        /// -------
        /// `matrix.general.Dense(T, order.invert())`:
        /// The transposed matrix.
        pub fn transpose(self: *const Dense(T, order)) Dense(T, order.invert()) {
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
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// The matrix to get the submatrix from.
        ///
        /// `row_start` (`u32`):
        /// The starting row index of the submatrix (inclusive).
        ///
        /// `row_end` (`u32`):
        /// The ending row index of the submatrix (exclusive). Must be greater
        /// than `row_start`.
        ///
        /// `col_start` (`u32`):
        /// The starting column index of the submatrix (inclusive).
        ///
        /// `col_end` (`u32`):
        /// The ending column index of the submatrix (exclusive). Must be
        /// greater than `col_start`.
        ///
        /// Returns
        /// -------
        /// `matrix.general.Dense(T, order)`:
        /// The submatrix.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.InvalidRange`:
        /// If the specified range is invalid.
        pub fn submatrix(
            self: *const Dense(T, order),
            row_start: u32,
            row_end: u32,
            col_start: u32,
            col_end: u32,
        ) !Dense(T, order) {
            if (row_start >= self.rows or col_start >= self.cols or
                row_end > self.rows or col_end > self.cols or
                row_start >= row_end or col_start >= col_end)
                return matrix.Error.InvalidRange;

            const sub_rows = row_end - row_start;
            const sub_cols = col_end - col_start;

            return .{
                .data = self.data + self._index(row_start, col_start),
                .rows = sub_rows,
                .cols = sub_cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the specified row as a dense vector.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// The matrix to get the row from.
        ///
        /// `r` (`u32`):
        /// The row index to get.
        ///
        /// Returns
        /// -------
        /// `vector.Dense(T)`:
        /// The specified row as a dense vector.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `r` is out of bounds.
        pub fn row(self: *const Dense(T, order), r: u32) !vector.Dense(T) {
            if (r >= self.rows)
                return matrix.Error.PositionOutOfBounds;

            return .{
                .data = self.data + self._index(r, 0),
                .len = self.cols,
                .inc = if (comptime order == .col_major)
                    types.scast(i32, self.ld)
                else
                    1,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the specified column as a dense vector.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// The matrix to get the column from.
        ///
        /// `c` (`u32`):
        /// The column index to get.
        ///
        /// Returns
        /// -------
        /// `vector.Dense(T)`:
        /// The specified column as a dense vector.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `c` is out of bounds.
        pub fn col(self: *const Dense(T, order), c: u32) !vector.Dense(T) {
            if (c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            return .{
                .data = self.data + self._index(0, c),
                .len = self.rows,
                .inc = if (comptime order == .col_major)
                    1
                else
                    types.scast(i32, self.ld),
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the matrix as a symmetric dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// The matrix to convert.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies whether the upper or lower triangle of the matrix is used,
        /// the other triangle is ignored.
        ///
        /// Returns
        /// -------
        /// `matrix.symmetric.Dense(T, uplo, order)`:
        /// The symmetric dense matrix view.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.NotSquare`:
        /// If the matrix is not square.
        pub fn asSymmetricDenseMatrix(self: *const Dense(T, order), comptime uplo: Uplo) !matrix.symmetric.Dense(T, uplo, order) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            return .{
                .data = self.data,
                .size = self.rows,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the matrix as a Hermitian dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// The matrix to convert.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies whether the upper or lower triangle of the matrix is used,
        /// the other triangle is ignored.
        ///
        /// Returns
        /// -------
        /// `matrix.hermitian.Dense(T, uplo, order)`:
        /// The Hermitian dense matrix view.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.NotSquare`:
        /// If the matrix is not square.
        pub fn asHermitianDenseMatrix(self: *const Dense(T, order), comptime uplo: Uplo) !matrix.hermitian.Dense(T, uplo, order) {
            comptime if (!types.isComplex(T))
                @compileError("Hermitian matrices require a complex type, got " ++ @typeName(T));

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            return .{
                .data = self.data,
                .size = self.rows,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the matrix as a triangular dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.general.Dense(T, order)`):
        /// The matrix to convert.
        ///
        /// `uplo` (`Uplo`):
        /// Specifies whether the upper or lower triangle of the matrix is used,
        /// the other triangle is ignored.
        ///
        /// `diag` (`Diag`):
        /// Specifies whether the matrix is unit triangular (diagonal elements
        /// are assumed to be 1 and are ignored) or non-unit triangular.
        ///
        /// Returns
        /// -------
        /// `matrix.triangular.Dense(T, uplo, diag, order)`:
        /// The triangular dense matrix view.
        pub fn asTriangularDenseMatrix(self: *const Dense(T, order), comptime uplo: Uplo, comptime diag: Diag) matrix.triangular.Dense(T, uplo, diag, order) {
            return .{
                .data = self.data,
                .rows = self.rows,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn asDenseArray(self: *const Dense(T, order)) array.Dense(T, order) {
            return .{
                .data = self.data,
                .ndim = 2,
                .shape = .{ self.rows, self.cols } ++ .{0} ** (array.max_dim - 2),
                .strides = if (comptime order == .col_major)
                    .{ 1, self.ld } ++ .{0} ** (array.max_dim - 2)
                else
                    .{ self.ld, 1 } ++ .{0} ** (array.max_dim - 2),
                .flags = .{
                    .order = self.flags.order,
                    .owns_data = false,
                },
            };
        }

        pub fn copyToDenseArray(
            self: *const Dense(T, order),
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
                        while (i < self.rows) : (i += 1) {
                            result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                        }
                    }
                } else {
                    var i: u32 = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: u32 = 0;
                        while (j < self.cols) : (j += 1) {
                            result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub inline fn _index(self: *const Dense(T, order), r: u32, c: u32) u32 {
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
        /// `self` (`*matrix.general.Dense(T, order)`):
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
        pub fn cleanup(self: *Dense(T, order), ctx: anytype) void {
            return _cleanup(self, self.rows * self.cols, order, ctx);
        }

        pub fn _cleanup(self: *Dense(T, order), num_elems: u32, iter_order: Order, ctx: anytype) void {
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

                    var n: u32 = 0;
                    if (iter_order == .col_major) {
                        var j: u32 = 0;
                        while (j < self.cols and n < num_elems) : (j += 1) {
                            var i: u32 = 0;
                            while (i < self.rows and n < num_elems) : (i += 1) {
                                ops.deinit(
                                    &self.data[self._index(i, j)],
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                );

                                n += 1;
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < self.rows and n < num_elems) : (i += 1) {
                            var j: u32 = 0;
                            while (j < self.cols and n < num_elems) : (j += 1) {
                                ops.deinit(
                                    &self.data[self._index(i, j)],
                                    types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                );

                                n += 1;
                            }
                        }
                    }
                },
            }
        }
    };
}
