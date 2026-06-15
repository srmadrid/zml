const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

const array = @import("../../array.zig");

/// Dense triangular matrix type, represented as a contiguous array of
/// `rows × cols` elements of type `N`, depending on `uplo`, stored in either
/// column-major or row-major order with a specified leading dimension. Only the
/// upper or lower triangular part of the matrix is accessed, depending on the
/// `uplo` parameter, and the diagonal can be either unit, meaning all diagonal
/// elements are assumed to be 1 and not accessed, or non-unit, meaning the
/// diagonal elements are accessed normally.
pub fn Dense(N: type, uplo: matrix.Uplo, diag: matrix.Diag, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.triangular.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        rows: usize,
        cols: usize,
        ld: usize,
        flags: matrix.Flags,

        // Type signatures
        pub const is_matrix = true;
        pub const is_dense = true;
        pub const is_triangular = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = uplo;
        pub const storage_diag: ?matrix.Diag = diag;

        // Numeric type
        pub const Numeric = N;

        pub const empty: Dense(N, uplo, diag, layout) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.triangular.Dense(N, uplo, diag, layout)`
        /// with the specified rows and columns.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the matrix.
        /// * `cols` (`usize`): The columns of the matrix.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The newly
        /// initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !matrix.triangular.Dense(N, uplo, diag, layout) {
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

        /// Initializes a new `matrix.triangular.Dense(N, uplo, diag, layout)`
        /// with the specified rows and columns, with all elements in the
        /// triangular part set to the specified value.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the matrix.
        /// * `cols` (`usize`): The columns of the matrix.
        /// * `value` (`N`): The value to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The newly
        /// initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn initValue(allocator: std.mem.Allocator, rows: usize, cols: usize, value: N) !matrix.triangular.Dense(N, uplo, diag, layout) {
            const mat: matrix.triangular.Dense(N, uplo, diag, layout) = try .init(allocator, rows, cols);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        var j: usize = 0;
                        while (j < cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j, rows)) : (i += 1) {
                                mat.data[i + j * mat.ld] = value;
                            }
                        }
                    } else { // cun
                        var j: usize = 0;
                        while (j < cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j + 1, rows)) : (i += 1) {
                                mat.data[i + j * mat.ld] = value;
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // clu
                        var j: usize = 0;
                        while (j < int.min(rows, cols)) : (j += 1) {
                            var i: usize = j + 1;
                            while (i < rows) : (i += 1) {
                                mat.data[i + j * mat.ld] = value;
                            }
                        }
                    } else { // cln
                        var j: usize = 0;
                        while (j < int.min(rows, cols)) : (j += 1) {
                            var i: usize = j;
                            while (i < rows) : (i += 1) {
                                mat.data[i + j * mat.ld] = value;
                            }
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // ruu
                        var i: usize = 0;
                        while (i < int.min(rows, cols)) : (i += 1) {
                            var j: usize = i + 1;
                            while (j < cols) : (j += 1) {
                                mat.data[i * mat.ld + j] = value;
                            }
                        }
                    } else { // run
                        var i: usize = 0;
                        while (i < int.min(rows, cols)) : (i += 1) {
                            var j: usize = i;
                            while (j < cols) : (j += 1) {
                                mat.data[i * mat.ld + j] = value;
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // rlu
                        var i: usize = 0;
                        while (i < rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i, cols)) : (j += 1) {
                                mat.data[i * mat.ld + j] = value;
                            }
                        }
                    } else { // rln
                        var i: usize = 0;
                        while (i < rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i + 1, cols)) : (j += 1) {
                                mat.data[i * mat.ld + j] = value;
                            }
                        }
                    }
                }
            }

            return mat;
        }

        /// Initializes a new `matrix.triangular.Dense(N, uplo, diag, layout)`
        /// with the specified rows and columns, with all elements in the
        /// triangular part set by calling the specified function with the given
        /// arguments.
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
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The newly
        /// initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn initFn(allocator: std.mem.Allocator, rows: usize, cols: usize, comptime @"fn": anytype, args: anytype) !matrix.triangular.Dense(N, uplo, diag, layout) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.triangular.Dense(N, uplo, diag, layout).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.triangular.Dense(N, uplo, diag, layout) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        var j: usize = 0;
                        while (j < cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j, rows)) : (i += 1) {
                                mat.data[i + j * mat.ld] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // cun
                        var j: usize = 0;
                        while (j < cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j + 1, rows)) : (i += 1) {
                                mat.data[i + j * mat.ld] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // clu
                        var j: usize = 0;
                        while (j < int.min(rows, cols)) : (j += 1) {
                            var i: usize = j + 1;
                            while (i < rows) : (i += 1) {
                                mat.data[i + j * mat.ld] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // cln
                        var j: usize = 0;
                        while (j < int.min(rows, cols)) : (j += 1) {
                            var i: usize = j;
                            while (i < rows) : (i += 1) {
                                mat.data[i + j * mat.ld] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // ruu
                        var i: usize = 0;
                        while (i < int.min(rows, cols)) : (i += 1) {
                            var j: usize = i + 1;
                            while (j < cols) : (j += 1) {
                                mat.data[i * mat.ld + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // run
                        var i: usize = 0;
                        while (i < int.min(rows, cols)) : (i += 1) {
                            var j: usize = i;
                            while (j < cols) : (j += 1) {
                                mat.data[i * mat.ld + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // rlu
                        var i: usize = 0;
                        while (i < rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i, cols)) : (j += 1) {
                                mat.data[i * mat.ld + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // rln
                        var i: usize = 0;
                        while (i < rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i + 1, cols)) : (j += 1) {
                                mat.data[i * mat.ld + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    }
                }
            }

            return mat;
        }

        /// Initializes a new identity
        /// `matrix.triangular.Dense(N, uplo, diag, layout)` of the specified
        /// size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The newly
        /// initialized identity matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initIdentity(allocator: std.mem.Allocator, size: usize) !matrix.triangular.Dense(N, uplo, diag, layout) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            const mat: matrix.triangular.Dense(N, uplo, diag, layout) = try .init(allocator, size, size);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            mat.data[i + j * mat.ld] = numeric.zero(N);
                        }

                        if (comptime diag == .non_unit) {
                            mat.data[j + j * mat.ld] = numeric.one(N);
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < size) : (j += 1) {
                        if (comptime diag == .non_unit) {
                            mat.data[j + j * mat.ld] = numeric.one(N);
                        }

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
                        if (comptime diag == .non_unit) {
                            mat.data[i * mat.ld + i] = numeric.one(N);
                        }

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

                        if (comptime diag == .non_unit) {
                            mat.data[i * mat.ld + i] = numeric.one(N);
                        }
                    }
                }
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.triangular.Dense(N, uplo, diag, layout)`): A
        ///   pointer to the matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.triangular.Dense(N, uplo, diag, layout), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Dense(N, uplo, diag, layout)`): The
        ///   matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.triangular.Dense(N, uplo, diag, layout), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (r > c)
                    return numeric.zero(N);
            } else {
                if (r < c)
                    return numeric.zero(N);
            }

            if (comptime diag == .unit) {
                if (r == c)
                    return numeric.one(N);
            }

            return self.data[self._index(r, c)];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Dense(N, uplo, diag, layout)`): The
        ///   matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds, on the correct triangular part, and outside the
        ///   diagonal if `diag` is unit.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds, on the correct triangular part, and outside the
        ///   diagonal if `diag` is unit.
        ///
        /// ## Returns
        /// `N`: The element at the specified position.
        pub fn getAssumeInBounds(self: matrix.triangular.Dense(N, uplo, diag, layout), r: usize, c: usize) N {
            return self.data[self._index(r, c)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.triangular.Dense(N, uplo, diag, layout)`): A
        ///   pointer to the matrix to set the element in.
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
        /// * `matrix.Error.BreaksStructure`: If `r == c` and `diag` is unit, or
        ///   if the index is outside the correct triangular part of the matrix.
        pub fn set(self: *matrix.triangular.Dense(N, uplo, diag, layout), r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (r > c)
                    return matrix.Error.BreaksStructure;
            } else {
                if (r < c)
                    return matrix.Error.BreaksStructure;
            }

            if (comptime diag == .unit) {
                if (r == c)
                    return matrix.Error.BreaksStructure;
            }

            self.data[self._index(r, c)] = value;
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.triangular.Dense(N, uplo, diag, layout)`): A
        ///   pointer to the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds, on the correct triangular part, and outside the
        ///   diagonal if `diag` is unit.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds, on the correct triangular part, and outside the
        ///   diagonal if `diag` is unit.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.triangular.Dense(N, uplo, diag, layout), r: usize, c: usize, value: N) void {
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
        pub fn setAll(self: *matrix.triangular.Dense(N, uplo, diag, layout), value: N) void {
            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        var j: usize = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j, self.rows)) : (i += 1) {
                                self.data[i + j * self.ld] = value;
                            }
                        }
                    } else { // cun
                        var j: usize = 0;
                        while (j < self.cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j + 1, self.rows)) : (i += 1) {
                                self.data[i + j * self.ld] = value;
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // clu
                        var j: usize = 0;
                        while (j < int.min(self.rows, self.cols)) : (j += 1) {
                            var i: usize = j + 1;
                            while (i < self.rows) : (i += 1) {
                                self.data[i + j * self.ld] = value;
                            }
                        }
                    } else { // cln
                        var j: usize = 0;
                        while (j < int.min(self.rows, self.cols)) : (j += 1) {
                            var i: usize = j;
                            while (i < self.rows) : (i += 1) {
                                self.data[i + j * self.ld] = value;
                            }
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // ruu
                        var i: usize = 0;
                        while (i < int.min(self.rows, self.cols)) : (i += 1) {
                            var j: usize = i + 1;
                            while (j < self.cols) : (j += 1) {
                                self.data[i * self.ld + j] = value;
                            }
                        }
                    } else { // run
                        var i: usize = 0;
                        while (i < int.min(self.rows, self.cols)) : (i += 1) {
                            var j: usize = i;
                            while (j < self.cols) : (j += 1) {
                                self.data[i * self.ld + j] = value;
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // rlu
                        var i: usize = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i, self.cols)) : (j += 1) {
                                self.data[i * self.ld + j] = value;
                            }
                        }
                    } else { // rln
                        var i: usize = 0;
                        while (i < self.rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i + 1, self.cols)) : (j += 1) {
                                self.data[i * self.ld + j] = value;
                            }
                        }
                    }
                }
            }
        }

        /// Creates a copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Dense(N, uplo, diag, layout)`): The
        ///   matrix to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The copied matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn clone(self: matrix.triangular.Dense(N, uplo, diag, layout), allocator: std.mem.Allocator) !matrix.triangular.Dense(N, uplo, diag, layout) {
            const mat: matrix.triangular.Dense(N, uplo, diag, layout) = try .init(allocator, self.rows, self.cols);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        var j: usize = 0;
                        while (j < mat.cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j, mat.rows)) : (i += 1) {
                                mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                            }
                        }
                    } else { // cun
                        var j: usize = 0;
                        while (j < mat.cols) : (j += 1) {
                            var i: usize = 0;
                            while (i < int.min(j + 1, mat.rows)) : (i += 1) {
                                mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // clu
                        var j: usize = 0;
                        while (j < int.min(mat.rows, mat.cols)) : (j += 1) {
                            var i: usize = j + 1;
                            while (i < mat.rows) : (i += 1) {
                                mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                            }
                        }
                    } else { // cln
                        var j: usize = 0;
                        while (j < int.min(mat.rows, mat.cols)) : (j += 1) {
                            var i: usize = j;
                            while (i < mat.rows) : (i += 1) {
                                mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                            }
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // ruu
                        var i: usize = 0;
                        while (i < int.min(mat.rows, mat.cols)) : (i += 1) {
                            var j: usize = i + 1;
                            while (j < mat.cols) : (j += 1) {
                                mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // run
                        var i: usize = 0;
                        while (i < int.min(mat.rows, mat.cols)) : (i += 1) {
                            var j: usize = i;
                            while (j < mat.cols) : (j += 1) {
                                mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // rlu
                        var i: usize = 0;
                        while (i < mat.rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i, mat.cols)) : (j += 1) {
                                mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                            }
                        }
                    } else { // rln
                        var i: usize = 0;
                        while (i < mat.rows) : (i += 1) {
                            var j: usize = 0;
                            while (j < int.min(i + 1, mat.cols)) : (j += 1) {
                                mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                            }
                        }
                    }
                }
            }

            return mat;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Dense(N, uplo, diag, layout)`): The
        ///   matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo.invert(), diag, layout.invert())`:
        /// The transposed matrix.
        pub fn transposeView(self: matrix.triangular.Dense(N, uplo, diag, layout)) matrix.triangular.Dense(N, uplo.invert(), diag, layout.invert()) {
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
        /// * `self` (`matrix.triangular.Dense(N, uplo, diag, layout)`): The
        ///   matrix to get the submatrix from.
        /// * `start` (`usize`): The starting diagonal index of the submatrix
        ///   (inclusive).
        /// * `row_end` (`usize`): The ending row index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        /// * `col_end` (`usize`): The ending column index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The submatrix.
        ///
        /// ## Errors
        /// * `matrix.Error.InvalidRange`: If the specified range is invalid.
        pub fn submatrixView(self: matrix.triangular.Dense(N, uplo, diag, layout), start: usize, row_end: usize, col_end: usize) !matrix.triangular.Dense(N, uplo, diag, layout) {
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            return .{
                .data = self.data + self._index(start, start),
                .rows = row_end - start,
                .cols = col_end - start,
                .ld = self.ld,
                .flags = .{ .owns_data = false },
            };
        }

        /// Copies the triangular matrix to a general dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Dense(N, uplo, diag, layout)`): The
        ///   matrix to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layour)`: The copied matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn cloneToGeneralDenseMatrix(self: matrix.triangular.Dense(N, uplo, diag, layout), allocator: std.mem.Allocator) !matrix.general.Dense(N, layout) {
            const mat: matrix.general.Dense(N, layout) = try .init(allocator, self.rows, self.cols);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                        }

                        if (comptime diag == .unit) {
                            mat.data[j + j * mat.ld] = numeric.one(N);
                        } else {
                            mat.data[j + j * mat.ld] = self.data[j + j * self.ld];
                        }

                        i = j + 1;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[i + j * mat.ld] = numeric.zero(N);
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i < j) : (i += 1) {
                            mat.data[i + j * mat.ld] = numeric.zero(N);
                        }

                        if (comptime diag == .unit) {
                            mat.data[j + j * mat.ld] = numeric.one(N);
                        } else {
                            mat.data[j + j * mat.ld] = self.data[j + j * self.ld];
                        }

                        i = j + 1;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.data[i + j * self.ld];
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j < i) : (j += 1) {
                            mat.data[i * mat.ld + j] = numeric.zero(N);
                        }

                        if (comptime diag == .unit) {
                            mat.data[i * mat.ld + i] = numeric.one(N);
                        } else {
                            mat.data[i * mat.ld + i] = self.data[i * self.ld + j];
                        }

                        j = i + 1;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j < i) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.data[i * self.ld + j];
                        }

                        if (comptime diag == .unit) {
                            mat.data[i * mat.ld + i] = numeric.one(N);
                        } else {
                            mat.data[i * mat.ld + i] = self.data[i * self.ld + i];
                        }

                        j = i + 1;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[i * mat.ld + j] = numeric.zero(N);
                        }
                    }
                }
            }

            return mat;
        }

        pub fn _index(self: *const Dense(N, uplo, diag, layout), r: usize, c: usize) usize {
            return if (comptime layout == .col_major)
                r + c * self.ld
            else
                r * self.ld + c;
        }
    };
}
