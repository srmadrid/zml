const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");
const array = @import("../../array.zig");

/// Dense general matrix type, represented as a contiguous array of
/// `rows × cols` elements of type `N`, stored in either column-major or
/// row-major order with a specified leading dimension.
pub fn Dense(N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.general.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        rows: usize,
        cols: usize,
        ld: usize,
        flags: matrix.Flags,

        // Type signatures
        pub const is_matrix = true;
        pub const is_dense = true;
        pub const is_general = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = null;
        pub const storage_diag: ?matrix.Diag = null;

        // Numeric type
        pub const Numeric = N;

        pub const empty: Dense(N, layout) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.general.Dense(N, layout)` with the
        /// specified rows and columns.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the matrix.
        /// * `cols` (`usize`): The columns of the matrix.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The newly initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !matrix.general.Dense(N, layout) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(N, rows * cols)).ptr,
                .rows = rows,
                .cols = cols,
                .ld = if (comptime layout == .col_major) rows else cols,
                .flags = .{ .owns_data = true },
            };
        }

        // pub fn initBuffer

        /// Initializes a new `matrix.general.Dense(N, layout)` with the
        /// specified rows and columns, with all elements set to the specified
        /// value.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the matrix.
        /// * `cols` (`usize`): The columns of the matrix.
        /// * `value` (`N`): The value to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The newly initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn initValue(allocator: std.mem.Allocator, rows: usize, cols: usize, value: N) !matrix.general.Dense(N, layout) {
            const mat: matrix.general.Dense(N, layout) = try .init(allocator, rows, cols);

            if (comptime layout == .col_major) {
                var j: usize = 0;
                while (j < cols) : (j += 1) {
                    var i: usize = 0;
                    while (i < rows) : (i += 1) {
                        mat.data[i + j * mat.ld] = value;
                    }
                }
            } else {
                var i: usize = 0;
                while (i < rows) : (i += 1) {
                    var j: usize = 0;
                    while (j < cols) : (j += 1) {
                        mat.data[i * mat.ld + j] = value;
                    }
                }
            }

            return mat;
        }

        /// Initializes a new `matrix.general.Dense(N, layout)` with the
        /// specified rows and columns, with all elements set by calling the
        /// specified function with the given arguments.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the matrix.
        /// * `cols` (`usize`): The columns of the matrix.
        /// * `@"fn"` (`anytype`): The function to call to fill the matrix.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The newly initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn initFn(allocator: std.mem.Allocator, rows: usize, cols: usize, comptime @"fn": anytype, args: anytype) !matrix.general.Dense(N, layout) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.general.Dense(N, layout).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.general.Dense(N, layout) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            if (comptime layout == .col_major) {
                var j: usize = 0;
                while (j < cols) : (j += 1) {
                    var i: usize = 0;
                    while (i < rows) : (i += 1) {
                        mat.data[i + j * mat.ld] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                            try @call(.auto, @"fn", args)
                        else
                            @call(.auto, @"fn", args);
                    }
                }
            } else {
                var i: usize = 0;
                while (i < rows) : (i += 1) {
                    var j: usize = 0;
                    while (j < cols) : (j += 1) {
                        mat.data[i * mat.ld + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                            try @call(.auto, @"fn", args)
                        else
                            @call(.auto, @"fn", args);
                    }
                }
            }

            return mat;
        }

        /// Initializes a new identity `matrix.general.Dense(N, layout)` of the
        /// specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The newly initialized identity
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initIdentity(allocator: std.mem.Allocator, size: usize) !matrix.general.Dense(N, layout) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            const mat: Dense(N, layout) = try .init(allocator, size, size);

            if (comptime layout == .col_major) {
                var j: usize = 0;
                while (j < size) : (j += 1) {
                    var i: usize = 0;
                    while (i < j) : (i += 1) {
                        mat.data[i + j * mat.ld] = numeric.zero(N);
                    }

                    mat.data[j + j * mat.ld] = numeric.one(N);

                    i += 1;

                    while (i < size) : (i += 1) {
                        mat.data[i + j * mat.ld] = numeric.zero(N);
                    }
                }
            } else {
                var i: usize = 0;
                while (i < size) : (i += 1) {
                    var j: usize = 0;
                    while (j < i) : (j += 1) {
                        mat.data[i * mat.ld + j] = numeric.zero(N);
                    }

                    mat.data[i * mat.ld + j] = numeric.one(N);

                    j += 1;

                    while (j < size) : (j += 1) {
                        mat.data[i * mat.ld + j] = numeric.zero(N);
                    }
                }
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Dense(N, layout)`):  A pointer to the
        ///   matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.general.Dense(N, layout), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.general.Dense(N, layout), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            return self.data[self._index(r, c)];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: matrix.general.Dense(N, layout), r: usize, c: usize) N {
            return self.data[self._index(r, c)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Dense(N, layout)`): A pointer to the
        ///   matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set.
        /// * `c` (`usize`): The column index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn set(self: *matrix.general.Dense(N, layout), r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            self.data[self._index(r, c)] = value;
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Dense(N, layout)`): A pointer to the
        ///   matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.general.Dense(N, layout), r: usize, c: usize, value: N) void {
            self.data[self._index(r, c)] = value;
        }

        /// Sets all elements of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Dense(N, layout)`): A pointer to the
        ///   matrix to set the elements in.
        /// * `value` (`N`): The value to set the elements to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAll(self: *matrix.general.Dense(N, layout), value: N) void {
            if (comptime layout == .col_major) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    var i: usize = 0;
                    while (i < self.rows) : (i += 1) {
                        self.data[i + j * self.ld] = value;
                    }
                }
            } else {
                var i: usize = 0;
                while (i < self.rows) : (i += 1) {
                    var j: usize = 0;
                    while (j < self.cols) : (j += 1) {
                        self.data[i * self.ld + j] = value;
                    }
                }
            }
        }

        /// Creates a copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The copied matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn clone(self: matrix.general.Dense(N, layout), allocator: std.mem.Allocator) !matrix.general.Dense(N, layout) {
            const mat: matrix.general.Dense(N, layout) = try .init(allocator, self.rows, self.cols);

            if (comptime layout == .col_major) {
                var j: usize = 0;
                while (j < mat.cols) : (j += 1) {
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                    }
                }
            } else {
                var i: usize = 0;
                while (i < mat.rows) : (i += 1) {
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                    }
                }
            }

            return mat;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to
        ///   transpose.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout.invert())`: The transposed matrix.
        pub fn transposeView(self: matrix.general.Dense(N, layout)) matrix.general.Dense(N, layout.invert()) {
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
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   submatrix from.
        /// * `row_start` (`usize`): The starting row index of the submatrix
        ///   (inclusive).
        /// * `row_end` (`usize`): The ending row index of the submatrix
        ///   (exclusive). Must be greater than `row_start`.
        /// * `col_start` (`usize`): The starting column index of the submatrix
        ///   (inclusive).
        /// * `col_end` (`usize`): The ending column index of the submatrix
        ///   (exclusive). Must be greater than `col_start`.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The submatrix.
        ///
        /// ## Errors
        /// * `matrix.Error.InvalidRange`: If the specified range is invalid.
        pub fn submatrixView(self: matrix.general.Dense(N, layout), row_start: usize, row_end: usize, col_start: usize, col_end: usize) !matrix.general.Dense(N, layout) {
            if (row_start >= self.rows or col_start >= self.cols or
                row_end > self.rows or col_end > self.cols or
                row_start >= row_end or col_start >= col_end)
                return matrix.Error.InvalidRange;

            return .{
                .data = self.data + self._index(row_start, col_start),
                .rows = row_end - row_start,
                .cols = col_end - col_start,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the specified row as a dense vector.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   row from.
        /// * `r` (`usize`): The row index to get.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The specified row as a dense vector.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` is out of bounds.
        pub fn rowView(self: matrix.general.Dense(N, layout), r: usize) !vector.Dense(N) {
            if (r >= self.rows)
                return matrix.Error.PositionOutOfBounds;

            return .{
                .data = self.data + self._index(r, 0),
                .len = self.cols,
                .inc = if (comptime layout == .col_major)
                    numeric.cast(isize, self.ld)
                else
                    1,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the specified column as a dense vector.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   column from.
        /// * `c` (`usize`): The column index to get.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The specified column as a dense vector.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `c` is out of bounds.
        pub fn colView(self: matrix.general.Dense(N, layout), c: usize) !vector.Dense(N) {
            if (c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            return .{
                .data = self.data + self._index(0, c),
                .len = self.rows,
                .inc = if (comptime layout == .col_major)
                    1
                else
                    numeric.cast(isize, self.ld),
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the matrix as a symmetric dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   view of.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The symmetric dense
        /// matrix view.
        ///
        /// ## Errors
        /// * `matrix.Error.NotSquare`: If the matrix is not square.
        pub fn symmetricView(self: matrix.general.Dense(N, layout), comptime uplo: matrix.Uplo) !matrix.symmetric.Dense(N, uplo, layout) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            return .{
                .data = self.data,
                .rows = self.rows,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the matrix as a Hermitian dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   view of.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The Hermitian dense
        /// matrix view.
        ///
        /// ## Errors
        /// * `matrix.Error.NotSquare`: If the matrix is not square.
        pub fn hermitianView(self: matrix.general.Dense(N, layout), comptime uplo: matrix.Uplo) !matrix.hermitian.Dense(N, uplo, layout) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            return .{
                .data = self.data,
                .rows = self.rows,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the matrix as a triangular dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Dense(N, layout)`): The matrix to get the
        ///   view of.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        /// * `diag` (`comptime matrix.Diag`): Specifies whether the matrix is
        ///   unit triangular (diagonal elements are assumed to be 1 and are
        ///   ignored) or non-unit triangular.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The triangular
        /// dense matrix view.
        pub fn triangularView(self: matrix.general.Dense(N, layout), comptime uplo: matrix.Uplo, comptime diag: matrix.Diag) matrix.triangular.Dense(N, uplo, diag, layout) {
            return .{
                .data = self.data,
                .rows = self.rows,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn _index(self: matrix.general.Dense(N, layout), r: usize, c: usize) usize {
            return if (comptime layout == .col_major)
                r + c * self.ld
            else
                r * self.ld + c;
        }
    };
}
