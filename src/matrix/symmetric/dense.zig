const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

const array = @import("../../array.zig");

const matutils = @import("../utils.zig");

/// Dense symmetric matrix type, represented as a contiguous array of
/// `size × size` elements of type `N`, stored in either column-major or
/// row-major order with a specified leading dimension. Only the upper or lower
/// triangular part of the matrix is accessed, depending on the `uplo`
/// parameter.
pub fn Dense(N: type, uplo: matrix.Uplo, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.symmetric.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        rows: usize,
        cols: usize,
        ld: usize,
        flags: matrix.Flags,

        // Type signatures
        pub const is_matrix = true;
        pub const is_dense = true;
        pub const is_symmetric = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = uplo;
        pub const storage_diag: ?matrix.Diag = null;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.symmetric.Dense(N, uplo, layout) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .ld = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.symmetric.Dense(N, uplo, layout)` with the
        /// specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn init(allocator: std.mem.Allocator, size: usize) !matrix.symmetric.Dense(N, uplo, layout) {
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

        /// Initializes a new `matrix.symmetric.Dense(N, uplo, layout)` with the
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
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initValue(allocator: std.mem.Allocator, size: usize, value: N) !matrix.symmetric.Dense(N, uplo, layout) {
            const mat: matrix.symmetric.Dense(N, uplo, layout) = try .init(allocator, size);

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

        /// Initializes a new `matrix.symmetric.Dense(N, uplo, layout)` with the
        /// specified size, with all elements in the stored triangular part set
        /// by calling the specified function with the given arguments.
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
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initFn(allocator: std.mem.Allocator, size: usize, comptime @"fn": anytype, args: anytype) !matrix.symmetric.Dense(N, uplo, layout) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.symmetric.Dense(N, uplo, layout).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.symmetric.Dense(N, uplo, layout) = try .init(allocator, size);
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

        /// Initializes a new identity `matrix.symmetric.Dense(N, uplo, layout)`
        /// of the specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The newly initialized
        /// identity matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initIdentity(allocator: std.mem.Allocator, size: usize) !matrix.symmetric.Dense(N, uplo, layout) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            const mat: matrix.symmetric.Dense(N, uplo, layout) = try .init(allocator, size);

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
        /// * `self` (`*matrix.symmetric.Dense(T, uplo, order)`): A pointer to
        ///   the matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///   `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.symmetric.Dense(N, uplo, layout), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.rows * self.cols]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Dense(N, uplo, layout)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.symmetric.Dense(N, uplo, layout), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: usize = r;
            var j: usize = c;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            }

            return if (self.flags.noconj)
                self.data[self._index(i, j)]
            else
                numeric.conj(self.data[self._index(i, j)]);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Dense(N, uplo, layout)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds and on the correct triangular part.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds and on the correct triangular part.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: matrix.symmetric.Dense(N, uplo, layout), r: usize, c: usize) N {
            return if (self.flags.noconj)
                self.data[self._index(r, c)]
            else
                numeric.conj(self.data[self._index(r, c)]);
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Dense(N, uplo, layout)`): A pointer to
        ///   the matrix to set the element in.
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
        pub fn set(self: *matrix.symmetric.Dense(N, uplo, layout), r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: usize = r;
            var j: usize = c;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            }

            self.data[self._index(i, j)] = if (self.flags.noconj) value else numeric.conj(value);
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Dense(N, uplo, layout)`): A pointer to
        ///   the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds and on the correct triangular part.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds and on the correct triangular part.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.symmetric.Dense(N, uplo, layout), r: usize, c: usize, value: N) void {
            self.data[self._index(r, c)] = if (self.flags.noconj) value else numeric.conj(value);
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
        pub fn setAll(self: *matrix.symmetric.Dense(N, uplo, layout), value: N) void {
            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            self.data[i + j * self.ld] = if (self.flags.noconj) value else numeric.conj(value);
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < self.cols) : (j += 1) {
                        var i: usize = j;
                        while (i < self.rows) : (i += 1) {
                            self.data[i + j * self.ld] = if (self.flags.noconj) value else numeric.conj(value);
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: usize = i;
                        while (j < self.cols) : (j += 1) {
                            self.data[i * self.ld + j] = if (self.flags.noconj) value else numeric.conj(value);
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < self.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            self.data[i * self.ld + j] = if (self.flags.noconj) value else numeric.conj(value);
                        }
                    }
                }
            }
        }

        /// Creates a copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Dense(N, uplo, layout)`): The matrix to
        ///   copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The copied matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn clone(self: matrix.symmetric.Dense(N, uplo, layout), allocator: std.mem.Allocator) !matrix.symmetric.Dense(N, uplo, layout) {
            const mat: matrix.symmetric.Dense(N, uplo, layout) = try .init(allocator, self.rows);

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.getAssumeInBounds(i, j);
                        }
                    }
                } else { // cl
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = j;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.getAssumeInBounds(i, j);
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = i;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.getAssumeInBounds(i, j);
                        }
                    }
                } else { // rl
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.getAssumeInBounds(i, j);
                        }
                    }
                }
            }

            return mat;
        }

        /// Creates a copy of the matrix with inverted `uplo`.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Dense(N, uplo, layout)`): The matrix to
        ///   copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo.invert(), layout)`: The copied
        /// matrix.
        ///
        /// ## Errors
        /// `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn cloneInverseUplo(self: matrix.symmetric.Dense(N, uplo, layout), allocator: std.mem.Allocator) !matrix.symmetric.Dense(N, uplo.invert(), layout) {
            const mat: matrix.symmetric.Dense(N, uplo.invert(), layout) = try .init(allocator, self.rows);

            if (comptime layout == .col_major) {
                if (comptime uplo.invert() == .upper) { // cl -> cu
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = 0;
                        while (i <= j) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.getAssumeInBounds(j, i);
                        }
                    }
                } else { // cu -> cl
                    var j: usize = 0;
                    while (j < mat.cols) : (j += 1) {
                        var i: usize = j;
                        while (i < mat.rows) : (i += 1) {
                            mat.data[i + j * mat.ld] = self.getAssumeInBounds(j, i);
                        }
                    }
                }
            } else {
                if (comptime uplo.invert() == .upper) { // rl -> ru
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = i;
                        while (j < mat.cols) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.getAssumeInBounds(j, i);
                        }
                    }
                } else { // ru -> rl
                    var i: usize = 0;
                    while (i < mat.rows) : (i += 1) {
                        var j: usize = 0;
                        while (j <= i) : (j += 1) {
                            mat.data[i * mat.ld + j] = self.getAssumeInBounds(j, i);
                        }
                    }
                }
            }

            return mat;
        }

        /// Returns a transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Dense(N, uplo, layout)`): The matrix to
        ///   transpose.
        ///
        /// ##  Returns
        /// `matrix.symmetric.Dense(N, uplo.invert(), layout.invert())`: The
        /// transposed matrix.
        pub fn transposeView(self: matrix.symmetric.Dense(N, uplo, layout)) matrix.symmetric.Dense(N, uplo.invert(), layout.invert()) {
            return .{
                .data = self.data,
                .rows = self.cols,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        /// Returns an adjoint view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Dense(N, uplo, layout)`): The matrix to
        ///   get the adjoint of.
        ///
        /// ##  Returns
        /// `matrix.symmetric.Dense(N, uplo.invert(), layout.invert())`: The
        /// adjoint matrix.
        pub fn adjointView(self: matrix.symmetric.Dense(N, uplo, layout)) matrix.symmetric.Dense(N, uplo.invert(), layout.invert()) {
            return .{
                .data = self.data,
                .rows = self.cols,
                .cols = self.cols,
                .ld = self.ld,
                .flags = .{ .owns_data = false, .noconj = !self.flags.noconj },
            };
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Dense(N, uplo, layout)`): The matrix to
        ///   get the submatrix from.
        /// * `start` (`usize`): The starting diagonal index of the submatrix
        ///   (inclusive).
        /// * `end` (`usize`): The ending diagonal index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The submatrix.
        pub fn submatrixView(self: matrix.symmetric.Dense(N, uplo, layout), start: usize, end: usize) !matrix.symmetric.Dense(N, uplo, layout) {
            if (start >= self.rows or end > self.rows or start >= end)
                return matrix.Error.InvalidRange;

            return .{
                .data = self.data + self._index(start, start),
                .rows = end - start,
                .cols = end - start,
                .ld = self.ld,
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        pub fn _index(self: *const Dense(N, uplo, layout), r: usize, c: usize) usize {
            return if (comptime layout == .col_major)
                r + c * self.ld
            else
                r * self.ld + c;
        }

        pub fn format(self: matrix.symmetric.Dense(N, uplo, layout), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const matrix.symmetric.Dense(N, uplo, layout), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .matrix = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                matrix: *const matrix.symmetric.Dense(N, uplo, layout),

                pub fn format(self: matrix.symmetric.Dense(N, uplo, layout).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    const rows = self.matrix.rows;
                    const cols = self.matrix.cols;

                    try writer.print("zsl.matrix.symmetric.Dense({s}, .{s}, .{s}) ({d} × {d}):\n\n", .{ @typeName(N), @tagName(uplo), @tagName(layout), rows, cols });

                    return matutils.format(self, num_fmt, rows, cols, writer);
                }
            };
        }
    };
}
