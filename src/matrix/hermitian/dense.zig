const std = @import("std");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

const array = @import("../../array.zig");

/// Dense Hermitian matrix type, represented as a contiguous array of
/// `size × size` elements of type `N`, stored in either column-major or
/// row-major order with a specified leading dimension. Only the upper or lower
/// triangular part of the matrix is accessed, depending on the `uplo`
/// parameter.
pub fn Dense(N: type, uplo: Uplo, layout: Layout) type {
    if (!meta.isNumeric(N) or !meta.isComplex(N))
        @compileError("zsl.matrix.hermitian.Dense: N must be a complex numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        rows: usize,
        cols: usize,
        ld: usize,
        flags: matrix.Flags = .{},

        // Type signatures
        pub const is_matrix = true;
        pub const is_dense = true;
        pub const is_hermitian = true;
        pub const storage_layout = layout;
        pub const storage_uplo = uplo;
        pub const storage_diag = meta.default_diag;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.hermitian.Dense(N, uplo, layout) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.hermitian.Dense(N, uplo, layout)` with the
        /// specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn init(allocator: std.mem.Allocator, size: usize) !matrix.hermitian.Dense(N, uplo, layout) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(N, size * size)).ptr,
                .rows = size,
                .cols = size,
                .ld = size,
                .flags = .{ .owns_data = true },
            };
        }

        // pub fn initBuffer

        /// Initializes a new `matrix.hermitian.Dense(N, uplo, layout)` with the
        /// specified size, with all elements in the stored triangular part set
        /// to the specified value.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        /// * `value` (`N`): The value to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initValue(allocator: std.mem.Allocator, size: usize, value: N) !matrix.hermitian.Dense(N, uplo, layout) {
            const mat: matrix.hermitian.Dense(N, uplo, layout) = try .init(allocator, size);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            mat.data[i + j * mat.ld] = value;
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        var i: usize = j;
                        while (i < size) : (i += 1) {
                            mat.data[i + j * mat.ld] = value;
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < size) : (i += 1) {
                        var j: usize = i;
                        while (j < size) : (j += 1) {
                            mat.data[i * mat.ld + j] = value;
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < size) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            mat.data[i * mat.ld + j] = value;
                        }
                    }
                }
            }

            return mat;
        }

        /// Initializes a new `matrix.hermitian.Dense(N, uplo, layout)` with the
        /// specified size, with all elements in the stored triangular part set
        /// by calling the specifie function with the given arguments.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        /// * `@"fn"` (`anytype`): The function to call to fill the matrix.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initFn(allocator: std.mem.Allocator, size: usize, comptime @"fn": anytype, args: anytype) !matrix.hermitian.Dense(N, uplo, layout) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.hermitian.Dense(N, uplo, layout).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.hermitian.Dense(N, uplo, layout) = try .init(allocator, size);
            errdefer mat.deinit(allocator);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            mat.data[i + j * mat.ld] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        var i: usize = j;
                        while (i < size) : (i += 1) {
                            mat.data[i + j * mat.ld] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < size) : (i += 1) {
                        var j: usize = i;
                        while (j < size) : (j += 1) {
                            mat.data[i * mat.ld + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < size) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            mat.data[i * mat.ld + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                }
            }

            return mat;
        }

        /// Initializes a new identity `matrix.hermitian.Dense(N, uplo, layout)`
        /// of the specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The newly initialized
        /// identity matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initIdentity(allocator: std.mem.Allocator, size: usize) !matrix.hermitian.Dense(N, uplo, layout) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            const mat: matrix.hermitian.Dense(N, uplo, layout) = try .init(allocator, size);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            mat.data[i + j * mat.ld] = numeric.zero(N);
                        }

                        mat.data[j + j * mat.ld] = numeric.one(N);
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        mat.data[j + j * mat.ld] = numeric.one(N);

                        var i: usize = j + 1;

                        while (i < size) : (i += 1) {
                            mat.data[i + j * mat.ld] = numeric.zero(N);
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < size) : (i += 1) {
                        mat.data[i * mat.ld + i] = numeric.one(N);

                        var j: usize = i + 1;

                        while (j < size) : (j += 1) {
                            mat.data[i * mat.ld + j] = numeric.zero(N);
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < size) : (i += 1) {
                        var j: usize = 0;
                        while (j < i) : (j += 1) {
                            mat.data[i * mat.ld + j] = numeric.zero(N);
                        }

                        mat.data[i * mat.ld + i] = numeric.one(N);
                    }
                }
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.hermitian.Dense(N, uplo, layout)`): A pointer to
        ///   the matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///   `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.hermitian.Dense(N, uplo, layout), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.hermitian.Dense(n, uplo, layout)`): The matrix to
        ///   get the element from.
        /// *  `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.hermitian.Dense(N, uplo, layout), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: usize = r;
            var j: usize = c;
            var noconj: bool = true;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                    noconj = false;
                }
            } else {
                if (i < j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                    noconj = false;
                }
            }

            return if (noconj)
                self.data[self._index(i, j)]
            else
                numeric.conj(self.data[self._index(i, j)]);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.hermitian.Dense(N, uplo, layout)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds and on the correct triangular part.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds and on the correct triangular part.
        ///
        /// Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: matrix.hermitian.Dense(N, uplo, layout), r: usize, c: usize) N {
            return self.data[self._index(r, c)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.hermitian.Dense(N, uplo, layout)`): The matrix to
        ///   set the element in.
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
        /// * `matrix.Error.BreaksStructure`: If `r == c` and the imaginary part
        ///   of `value` is not zero.
        pub fn set(self: *matrix.hermitian.Dense(N, uplo, layout), r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (r == c and numeric.ne(value.im, 0))
                return matrix.Error.BreaksStructure;

            var i: usize = r;
            var j: usize = c;
            var conj: bool = false;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                    conj = true;
                }
            } else {
                if (i < j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                    conj = true;
                }
            }

            self.data[self._index(i, j)] = value;

            if (conj) {
                numeric.conjInto(&self.data[self._index(i, j)], self.data[self._index(i, j)]);
            }
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.hermitian.Dense(N, uplo, layout)`): A pointer to
        ///   the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds and on the correct triangular part.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds and on the correct triangular part.
        /// * `value` (`N`): The value to set the element to. Assumed to have
        ///   zero imaginary part if `r == c`.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.hermitian.Dense(N, uplo, layout), r: usize, c: usize, value: N) void {
            self.data[self._index(r, c)] = value;
        }

        /// Sets all elements of the stored triangle of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Dense(N, layout)`): A pointer to the
        ///   matrix to set the elements in.
        /// * `value` (`N`): The value to set the elements to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAll(self: *matrix.hermitian.Dense(N, uplo, layout), value: N) void {
            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            self.data[i + j * self.ld] = value;
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: usize = j;
                        while (i < self.rows) : (i += 1) {
                            self.data[i + j * self.ld] = value;
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: usize = i;
                        while (j < self.cols) : (j += 1) {
                            self.data[i * self.ld + j] = value;
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            self.data[i * self.ld + j] = value;
                        }
                    }
                }
            }
        }

        /// Creates a copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.hermitian.Dense(N, uplo, layout)`): The matrix to
        ///   copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The copied matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn copy(self: matrix.hermitian.Dense(N, uplo, layout), allocator: std.mem.Allocator) !matrix.hermitian.Dense(N, uplo, layout) {
            const mat: matrix.hermitian.Dense(N, uplo, layout) = try .init(allocator, self.rows);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = j;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = i;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                        }
                    }
                }
            }

            return mat;
        }

        /// Creates a copy of the matrix with inverted `uplo`.
        ///
        /// ## Arguments
        /// * `self` (`matrix.hermitian.Dense(N, uplo, layout)`): The matrix to
        ///   copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo.invert(), layout)`: The copied
        /// matrix.
        ///
        /// ## Errors
        /// `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn copyInverseUplo(self: matrix.hermitian.Dense(N, uplo, layout), allocator: std.mem.Allocator) !matrix.hermitian.Dense(N, uplo.invert(), layout) {
            const mat: matrix.hermitian.Dense(N, uplo.invert(), layout) = try .init(allocator, self.rows);

            if (comptime layout == .col_major) {
                if (comptime uplo.invert() == .upper) { // cl -> cu
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            numeric.conjInto(&mat.data[i + j * mat.ld], self.data[j + i * self.ld]);
                        }
                    }
                } else { // cu -> cl
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = j;
                        while (i < mat.rows) : (i += 1) {
                            numeric.conjInto(&mat.data[i + j * mat.ld], self.data[j + i * self.ld]);
                        }
                    }
                }
            } else {
                if (comptime uplo.invert() == .upper) { // rl -> ru
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = i;
                        while (j < mat.cols) : (j += 1) {
                            numeric.conjInto(&mat.data[i * mat.ld + j], self.data[j * self.ld + i]);
                        }
                    }
                } else { // ru -> rl
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            numeric.conjInto(&mat.data[i * mat.ld + j], self.data[j * self.ld + i]);
                        }
                    }
                }
            }

            return mat;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.hermitian.Dense(N, uplo, layout)`): The matrix to
        ///   transpose.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, uplo.invert(), layout.invert())`: The
        /// transposed matrix.
        pub fn transpose(self: matrix.hermitian.Dense(N, uplo, layout)) matrix.hermitian.Dense(N, uplo.invert(), layout.invert()) {
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
        /// * `self` (`matrix.hermitian.Dense(N, uplo, layout)`): The matrix to
        ///   get the submatrix from.
        /// * `start` (`usize`): The starting diagonal index of the submatrix
        ///   (inclusive).
        /// * `end` (`usize`): The ending diagonal index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The submatrix.
        ///
        /// Errors
        /// * `matrix.Error.InvalidRange`: If the specified range is invalid.
        pub fn submatrix(self: matrix.hermitian.Dense(N, uplo, layout), start: usize, end: usize) !matrix.hermitian.Dense(N, uplo, layout) {
            if (start >= self.rows or end > self.rows or start >= end)
                return matrix.Error.InvalidRange;

            return .{
                .data = self.data + self._index(start, start),
                .rows = end - start,
                .cols = end - start,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Copies the hermitian matrix to a general dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.hermitian.Dense(N, uplo, layout)`): The matrix to
        ///   copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The copied matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn copyToGeneralDenseMatrix(self: matrix.hermitian.Dense(N, uplo, layout), allocator: std.mem.Allocator) !matrix.general.Dense(N, layout) {
            const mat: matrix.general.Dense(N, layout) = try .init(allocator, self.rows, self.cols);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                            numeric.conjInto(&mat.data[j + i * mat.ld], self.data[i + j * self.ld]);
                        }

                        mat.data[j + j * mat.ld] = self.data[j + j * self.ld];
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        mat.data[j + j * mat.ld] = self.data[j + j * self.ld];

                        var i: usize = j + 1;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                            numeric.conjInto(&mat.data[j + i * mat.ld], self.data[i + j * self.ld]);
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        mat.data[i * mat.ld + i] = self.data[i * self.ld + i];

                        var j: usize = i + 1;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                            numeric.conjInto(&mat.data[j * mat.ld + i], self.data[i * self.ld + j]);
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j < i) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                            numeric.conjInto(&mat.data[j * mat.ld + i], self.data[i * self.ld + j]);
                        }

                        mat.data[i * mat.ld + i] = self.data[i * self.ld + i];
                    }
                }
            }

            return mat;
        }

        // pub fn copyToDenseArray(
        //     self: *const Dense(N, uplo, layout),
        //     allocator: std.mem.Allocator,
        //     ctx: anytype,
        // ) !array.Dense(N, layout) {
        //     var result: array.Dense(N, layout) = try .init(allocator, &.{ self.rows, self.cols });
        //     errdefer result.deinit(allocator);

        //     if (comptime !types.isArbitraryPrecision(N)) {
        //         comptime types.validateContext(@TypeOf(ctx), .{});

        //         if (comptime layout == .col_major) {
        //             if (comptime uplo == .upper) { // cu
        //                 var j: usize = 0;
        //                 while (j < self.cols) : (j += 1) {
        //                     var i: usize = 0;
        //                     while (i < j) : (i += 1) {
        //                         result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
        //                         result.data[j + i * result.strides[1]] = numeric.conj(self.data[i + j * self.ld], ctx) catch unreachable;
        //                     }

        //                     result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];
        //                 }
        //             } else {
        //                 var j: usize = 0;
        //                 while (j < self.cols) : (j += 1) {
        //                     result.data[j + j * result.strides[1]] = self.data[j + j * self.ld];

        //                     var i: usize = j + 1;
        //                     while (i < self.rows) : (i += 1) {
        //                         result.data[i + j * result.strides[1]] = self.data[i + j * self.ld];
        //                         result.data[j + i * result.strides[1]] = numeric.conj(self.data[i + j * self.ld], ctx) catch unreachable;
        //                     }
        //                 }
        //             }
        //         } else {
        //             if (comptime uplo == .upper) { // ru
        //                 var i: usize = 0;
        //                 while (i < self.rows) : (i += 1) {
        //                     result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];

        //                     var j: usize = i + 1;
        //                     while (j < self.cols) : (j += 1) {
        //                         result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
        //                         result.data[j * result.strides[0] + i] = numeric.conj(self.data[i * self.ld + j], ctx) catch unreachable;
        //                     }
        //                 }
        //             } else { // rl
        //                 var i: usize = 0;
        //                 while (i < self.rows) : (i += 1) {
        //                     var j: usize = 0;
        //                     while (j < i) : (j += 1) {
        //                         result.data[i * result.strides[0] + j] = self.data[i * self.ld + j];
        //                         result.data[j * result.strides[0] + i] = numeric.conj(self.data[i * self.ld + j], ctx) catch unreachable;
        //                     }

        //                     result.data[i * result.strides[0] + i] = self.data[i * self.ld + i];
        //                 }
        //             }
        //         }
        //     } else {
        //         @compileError("Arbitrary precision types not implemented yet");
        //     }

        //     return result;
        // }

        pub fn _index(self: *const Dense(N, uplo, layout), r: usize, c: usize) usize {
            return if (comptime layout == .col_major)
                r + c * self.ld
            else
                r * self.ld + c;
        }
    };
}
