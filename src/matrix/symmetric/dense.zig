const std = @import("std");

const types = @import("../../types.zig");
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

/// Dense symmetric matrix type, represented as a contiguous array of
/// `size × size` elements of type `T`, stored in either column-major or
/// row-major order with a specified leading dimension. Only the upper or lower
/// triangular part of the matrix is accessed, depending on the `uplo`
/// parameter.
pub fn Dense(T: type, uplo: Uplo, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("matrix.symmetric.Dense requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        ld: u32, // leading dimension
        flags: Flags = .{},

        pub const empty: Dense(T, uplo, order) = .{
            .data = &.{},
            .size = 0,
            .ld = 0,
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
        /// `matrix.symmetric.Dense(T, uplo, order)`:
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
        ) !Dense(T, uplo, order) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(T, size * size)).ptr,
                .size = size,
                .ld = size,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new matrix with the specified size, filled  with the
        /// specified value.
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
        /// `matrix.symmetric.Dense(T, uplo, order)`:
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
        ) !Dense(T, uplo, order) {
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

            var mat: Dense(T, uplo, order) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                if (comptime uplo == .upper) { // cu
                    while (j < size) : (j += 1) {
                        i = 0;
                        while (i <= j) : (i += 1) {
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
                } else { // cl
                    while (j < size) : (j += 1) {
                        i = j;
                        while (i < size) : (i += 1) {
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
                if (comptime uplo == .upper) { // ru
                    while (i < size) : (i += 1) {
                        j = i;
                        while (j < size) : (j += 1) {
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
                } else { // rl
                    while (i < size) : (i += 1) {
                        j = 0;
                        while (j <= i) : (j += 1) {
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
        /// `matrix.symmetric.Dense(T, uplo, order)`:
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
        ) !Dense(T, uplo, order) {
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

            var mat: Dense(T, uplo, order) = try .init(allocator, size);
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

                        mat.data[mat._index(j, j)] = try constants.one(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                } else { // cl
                    while (j < size) : (j += 1) {
                        mat.data[mat._index(j, j)] = try constants.one(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );

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
                        mat.data[mat._index(i, i)] = try constants.one(
                            T,
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );

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

                        mat.data[mat._index(i, i)] = try constants.one(
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
        /// `self` (`*matrix.symmetric.Dense(T, uplo, order)`):
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
        pub fn deinit(self: *Dense(T, uplo, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.size * self.size]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.symmetric.Dense(T, uplo, order)`):
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
        pub fn get(self: *const Dense(T, uplo, order), r: u32, c: u32) !T {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = r;
            var j: u32 = c;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            }

            return self.data[self._index(i, j)];
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.symmetric.Dense(T, uplo, order)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be within bounds and
        /// on the correct triangular part.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be within bounds
        /// and on the correct triangular part.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub inline fn at(self: *const Dense(T, uplo, order), r: u32, c: u32) T {
            return self.data[self._index(r, c)];
        }

        /// Sets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.symmetric.Dense(T, uplo, order)`):
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
        pub fn set(self: *Dense(T, uplo, order), row: u32, col: u32, value: T) !void {
            if (row >= self.size or col >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = row;
            var j: u32 = col;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: u32 = i;
                    i = j;
                    j = temp;
                }
            }

            self.data[self._index(i, j)] = value;
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.symmetric.Dense(T, uplo, order)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be within bounds and
        /// on the correct triangular part.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be within bounds
        /// and on the correct triangular part.
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
        pub inline fn put(self: *Dense(T, uplo, order), r: u32, c: u32, value: T) void {
            self.data[self._index(r, c)] = value;
        }

        /// Creates a copy of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.symmetric.Dense(T, uplo, order)`):
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
        /// `matrix.symmetric.Dense(T, uplo, order)`:
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
        pub fn copy(self: *const Dense(T, uplo, order), allocator: std.mem.Allocator, ctx: anytype) !Dense(T, uplo, order) {
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

            var mat: Dense(T, uplo, order) = try .init(allocator, self.size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                if (comptime uplo == .upper) { // cu
                    while (j < mat.size) : (j += 1) {
                        i = 0;
                        while (i <= j) : (i += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // cl
                    while (j < mat.size) : (j += 1) {
                        i = j;
                        while (i < mat.size) : (i += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    while (i < mat.size) : (i += 1) {
                        j = i;
                        while (j < mat.size) : (j += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // rl
                    while (i < mat.size) : (i += 1) {
                        j = 0;
                        while (j <= i) : (j += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            }

            return mat;
        }

        /// Creates a copy of the matrix with inverted `uplo`.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.symmetric.Dense(T, uplo, order)`):
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
        /// `matrix.symmetric.Dense(T, uplo.invert(), order)`:
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
        pub fn copyInverseUplo(
            self: *const Dense(T, uplo, order),
            allocator: std.mem.Allocator,
            ctx: anytype,
        ) !Dense(T, uplo.invert(), order) {
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

            var mat: Dense(T, uplo.invert(), order) = try .init(allocator, self.size);
            errdefer mat.deinit(allocator);

            var i: u32 = 0;
            var j: u32 = 0;

            errdefer mat._cleanup(i, j, order, ctx);

            if (comptime order == .col_major) {
                if (comptime uplo.invert() == .upper) { // cl -> cu
                    while (j < mat.size) : (j += 1) {
                        i = 0;
                        while (i <= j) : (i += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(j, i)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // cu -> cl
                    while (j < mat.size) : (j += 1) {
                        i = j;
                        while (i < mat.size) : (i += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(j, i)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            } else {
                if (comptime uplo.invert() == .upper) { // rl -> ru
                    while (i < mat.size) : (i += 1) {
                        j = i;
                        while (j < mat.size) : (j += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(j, i)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                } else { // ru -> rl
                    while (i < mat.size) : (i += 1) {
                        j = 0;
                        while (j <= i) : (j += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(j, i)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            }

            return mat;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.symmetric.Dense(T, uplo, order)`):
        /// The matrix to transpose.
        ///
        /// Returns
        /// -------
        /// `matrix.general.Dense(T, uplo.invert(), order.invert())`:
        /// The transposed matrix.
        pub fn transpose(self: Dense(T, uplo, order)) Dense(T, uplo.invert(), order.invert()) {
            return .{
                .data = self.data,
                .size = self.size,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.symmetric.Dense(T, uplo, order)`):
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
        /// `matrix.symmetric.Dense(T, uplo, order)`:
        /// The submatrix.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.InvalidRange`:
        /// If the specified range is invalid.
        pub fn submatrix(
            self: *const Dense(T, uplo, order),
            start: u32,
            end: u32,
        ) !Dense(T, uplo, order) {
            if (start >= self.size or end > self.size or start >= end)
                return matrix.Error.InvalidRange;

            const sub_size: u32 = end - start;

            return .{
                .data = self.data + self._index(start, start),
                .size = sub_size,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Copies the symmetric matrix to a general dense matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.symmetric.Dense(T, uplo, order)`):
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
            self: Dense(T, uplo, order),
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

                            mat.data[mat._index(j, i)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        mat.data[mat._index(j, j)] = try ops.copy(
                            self.data[self._index(j, j)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                } else { // cl
                    while (j < mat.cols) : (j += 1) {
                        mat.data[mat._index(j, j)] = try ops.copy(
                            self.data[self._index(j, j)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );

                        i = j + 1;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );

                            mat.data[mat._index(j, i)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    while (i < mat.rows) : (i += 1) {
                        mat.data[mat._index(i, i)] = try ops.copy(
                            self.data[self._index(i, i)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );

                        j = i + 1;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[mat._index(i, j)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );

                            mat.data[mat._index(j, i)] = try ops.copy(
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

                            mat.data[mat._index(j, i)] = try ops.copy(
                                self.data[self._index(i, j)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }

                        mat.data[mat._index(i, i)] = try ops.copy(
                            self.data[self._index(i, i)],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                }
            }

            return mat;
        }

        pub fn copyToDenseArray(
            self: *const Dense(T, uplo, order),
            allocator: std.mem.Allocator,
            ctx: anytype,
        ) !array.Dense(T, order) {
            var result: array.Dense(T, order) = try .init(allocator, &.{ self.size, self.size });
            errdefer result.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime order == .col_major) {
                    if (comptime uplo == .upper) { // cu
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                                result.data[j + i * result.strides[1]] = self.data[i + j * self.ld];
                            }

                            result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
                        }
                    } else { // cl
                        var j: u32 = 0;
                        while (j < self.size) : (j += 1) {
                            result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];

                            var i: u32 = j + 1;
                            while (i < self.size) : (i += 1) {
                                result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
                                result.data[j + i * result.strides[1]] = self.data[i + j * self.ld];
                            }
                        }
                    }
                } else {
                    if (comptime uplo == .upper) { // ru
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];

                            var j: u32 = i + 1;
                            while (j < self.size) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                                result.data[j * result.strides[0] + i] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // rl
                        var i: u32 = 0;
                        while (i < self.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
                                result.data[j * result.strides[0] + i] = self.data[i * self.ld + j];
                            }

                            result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];
                        }
                    }
                }
            } else {
                @compileError("Arbitrary precision types not implemented yet");
            }

            return result;
        }

        pub inline fn _index(self: *const Dense(T, uplo, order), r: u32, c: u32) u32 {
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
        /// `self` (`*matrix.symmetric.Dense(T, order)`):
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
        pub fn cleanup(self: *Dense(T, uplo, order), ctx: anytype) void {
            return self._cleanup(self.size, self.size, order, ctx);
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary, in the specified order up to position (i, j), exclusive.
        /// In other words, if iter_order is column-major, all elements from
        /// (0, 0) to (i - 1, j) are cleaned up, and if iter_order is row-major,
        /// all elements from (0, 0) to (i, j - 1) are cleaned up.
        pub fn _cleanup(self: *Dense(T, uplo, order), i: u32, j: u32, iter_order: Order, ctx: anytype) void {
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
                            while (_j <= int.min(j, self.size - 1)) : (_j += 1) {
                                var _i: u32 = 0;
                                if (_j == j) {
                                    while (_i < int.min(i, _j + 1)) : (_i += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    while (_i <= _j) : (_i += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                }
                            }
                        } else {
                            var _j: u32 = 0;
                            while (_j <= int.min(j, self.size - 1)) : (_j += 1) {
                                var _i: u32 = _j;
                                if (_j == j) {
                                    while (_i < int.min(i, self.size)) : (_i += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    while (_i < self.size) : (_i += 1) {
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
                            while (_i <= int.min(i, self.size - 1)) : (_i += 1) {
                                var _j: u32 = _i;
                                if (_i == i) {
                                    while (_j < int.min(j, self.size)) : (_j += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    while (_j < self.size) : (_j += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                }
                            }
                        } else {
                            var _i: u32 = 0;
                            while (_i <= int.min(i, self.size - 1)) : (_i += 1) {
                                var _j: u32 = 0;
                                if (_i == i) {
                                    while (_j < int.min(j, _i + 1)) : (_j += 1) {
                                        ops.deinit(
                                            &self.data[self._index(_i, _j)],
                                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                                        );
                                    }
                                } else {
                                    while (_j <= _i) : (_j += 1) {
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
