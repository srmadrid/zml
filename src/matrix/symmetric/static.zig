const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const matrix = @import("../../matrix.zig");

const matutils = @import("../utils.zig");

/// Static symmetric matrix type, represented as a contiguous array of
/// `size × size` elements of type `N`, stored in either column-major or
/// row-major order. Only the upper or lower triangular part of the matrix is
/// accessed, depending on the `uplo` parameter.
pub fn Static(size_: comptime_int, N: type, uplo: matrix.Uplo, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.symmetric.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [size * size]N,

        // Type signatures
        pub const is_matrix = true;
        pub const is_static = true;
        pub const is_symmetric = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = uplo;
        pub const storage_diag: ?matrix.Diag = null;

        pub const size = size_;
        pub const rows = size_;
        pub const cols = size_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.symmetric.Static(size, N, uplo, layout) = .{
            .data = undefined,
        };

        pub const init = empty;

        /// Initializes a new `matrix.symmetric.Static(size, N, uplo, layout)`,
        /// with all elements in the stored triangular part set to the specified
        /// value.
        ///
        /// ## Arguments
        /// * `value` (`N`): The value to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.symmetric.Static(size, N, uplo, layout)`: The newly
        /// initialized matrix.
        pub fn initValue(value: N) matrix.symmetric.Static(size, N, uplo, layout) {
            var mat: matrix.symmetric.Static(size, N, uplo, layout) = .init;

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    inline for (0..size) |j| {
                        inline for (0..j + 1) |i| {
                            mat.data[i + j * size] = value;
                        }
                    }
                } else { // cl
                    inline for (0..size) |j| {
                        inline for (j..size) |i| {
                            mat.data[i + j * size] = value;
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    inline for (0..size) |i| {
                        inline for (i..size) |j| {
                            mat.data[i * size + j] = value;
                        }
                    }
                } else { // rl
                    inline for (0..size) |i| {
                        inline for (0..i + 1) |j| {
                            mat.data[i * size + j] = value;
                        }
                    }
                }
            }

            return mat;
        }

        /// Initializes a new `matrix.symmetric.Static(size, N, uplo, layout)`,
        /// with all elements set from the specified values.
        ///
        /// ## Arguments
        /// * `values` (`[size * size]N`): The values to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.symmetric.Static(size, N, uplo, layout)`: The newly
        /// initialized matrix.
        pub fn initArray(values: [size * size]N) matrix.symmetric.Static(size, N, uplo, layout) {
            return .{ .data = values };
        }

        /// Initializes a new `matrix.symmetric.Static(size, N, uplo, layout)`,
        /// with all elements in the stored triangular part set by calling the
        /// specified function with the given arguments.
        ///
        /// ## Arguments
        /// * `@"fn"` (`anytype`): The function to call to fill the matrix.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `matrix.symmetric.Static(size, N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        pub fn initFn(comptime @"fn": anytype, args: anytype) !matrix.symmetric.Static(size, N, uplo, layout) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.symmetric.Static(size, N, uplo, layout).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.symmetric.Static(size, N, uplo, layout) = .init;

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) { // cu
                    inline for (0..size) |j| {
                        inline for (0..j + 1) |i| {
                            mat.data[i + j * size] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                } else { // cl
                    inline for (0..size) |j| {
                        inline for (j..size) |i| {
                            mat.data[i + j * size] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) { // ru
                    inline for (0..size) |i| {
                        inline for (i..size) |j| {
                            mat.data[i * size + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                } else { // rl
                    inline for (0..size) |i| {
                        inline for (0..i + 1) |j| {
                            mat.data[i * size + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                try @call(.auto, @"fn", args)
                            else
                                @call(.auto, @"fn", args);
                        }
                    }
                }
            }

            return mat;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Static(size, N, uplo, layout)`): The
        ///   matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of bounds.
        pub fn get(self: matrix.symmetric.Static(size, N, uplo, layout), r: usize, c: usize) !N {
            if (r >= size or c >= size)
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

            return self.data[_index(i, j)];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Static(size, N, uplo, layout)`): The
        ///   matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds and on the correct triangular part.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds and on the correct triangular part.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: matrix.symmetric.Static(size, N, uplo, layout), r: usize, c: usize) N {
            return self.data[_index(r, c)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Static(size, N, uplo, layout)`): A
        ///   pointer to the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set.
        /// * `c` (`usize`): The column index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of bounds.
        pub fn set(self: *matrix.symmetric.Static(size, N, uplo, layout), r: usize, c: usize, value: N) !void {
            if (r >= size or c >= size)
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

            self.data[_index(i, j)] = value;
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Static(size, N, uplo, layout)`): A
        ///   pointer to the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds and on the correct triangular part.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds and on the correct triangular part.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.symmetric.Static(size, N, uplo, layout), r: usize, c: usize, value: N) void {
            self.data[_index(r, c)] = value;
        }

        /// Creates a copy of the matrix with inverted `uplo`.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Static(size, N, uplo, layout)`): The matrix to
        ///   copy.
        ///
        /// ## Returns
        /// `matrix.symmetric.Static(size, N, uplo.invert(), layout)`: The copied
        /// matrix.
        pub fn copyInverseUplo(self: matrix.symmetric.Static(size, N, uplo, layout)) matrix.symmetric.Static(size, N, uplo.invert(), layout) {
            var mat: matrix.symmetric.Static(size, N, uplo.invert(), layout) = .init;

            if (comptime layout == .col_major) {
                if (comptime uplo.invert() == .upper) { // cl -> cu
                    inline for (0..size) |j| {
                        inline for (0..j + 1) |i| {
                            mat.data[i + j * size] = self.data[j + i * size];
                        }
                    }
                } else { // cu -> cl
                    inline for (0..size) |j| {
                        inline for (j..size) |i| {
                            mat.data[i + j * size] = self.data[j + i * size];
                        }
                    }
                }
            } else {
                if (comptime uplo.invert() == .upper) { // rl -> ru
                    inline for (0..size) |i| {
                        inline for (i..size) |j| {
                            mat.data[i * size + j] = self.data[j * size + i];
                        }
                    }
                } else { // ru -> rl
                    inline for (0..size) |i| {
                        inline for (0..i + 1) |j| {
                            mat.data[i * size + j] = self.data[j * size + i];
                        }
                    }
                }
            }

            return mat;
        }

        /// Returns the transposed copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Static(size, N, uplo, layout)`): The matrix to
        ///   transpose.
        ///
        /// ## Returns
        /// `matrix.symmetric.Static(size, N, uplo.invert(), layout.invert())`: The
        /// transposed matrix.
        pub fn transposeCopy(self: matrix.symmetric.Static(size, N, uplo, layout)) matrix.symmetric.Static(size, N, uplo.invert(), layout.invert()) {
            return .{ .data = self.data };
        }

        /// Returns a transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Static(size, N, uplo, layout)`): A
        ///   pointer to the matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo.invert(), layout.invert())`: The
        /// transposed dense view.
        pub fn transposeView(self: *matrix.symmetric.Static(size, N, uplo, layout)) matrix.symmetric.Dense(N, uplo.invert(), layout.invert()) {
            return .{
                .data = &self.data,
                .rows = size,
                .cols = size,
                .ld = size,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a submatrix copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Static(size, N, uplo, layout)`): The
        ///   matrix to get the submatrix from.
        /// * `start` (`comptime usize`): The starting diagonal index of the
        ///   submatrix (inclusive).
        /// * `end` (`comptime usize`): The ending diagonal index of the
        ///   submatrix (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.symmetric.Static(end - start, N, uplo, layout)`: The copied
        /// submatrix.
        pub fn submatrixCopy(
            self: matrix.symmetric.Static(size, N, uplo, layout),
            comptime start: usize,
            comptime end: usize,
        ) matrix.symmetric.Static(end - start, N, uplo, layout) {
            comptime {
                if (start >= size or end > size or start >= end)
                    @compileError("zsl.matrix.symmetric.Static(size, N, uplo, layout).submatrixCopy: Invalid range");
            }

            const new_size = end - start;
            var mat: matrix.symmetric.Static(new_size, N, uplo, layout) = .init;

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    inline for (0..new_size) |j| {
                        inline for (0..j + 1) |i| {
                            mat.data[i + j * new_size] = self.data[(start + i) + (start + j) * size];
                        }
                    }
                } else {
                    inline for (0..new_size) |j| {
                        inline for (j..new_size) |i| {
                            mat.data[i + j * new_size] = self.data[(start + i) + (start + j) * size];
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) {
                    inline for (0..new_size) |i| {
                        inline for (i..new_size) |j| {
                            mat.data[i * new_size + j] = self.data[(start + i) * size + (start + j)];
                        }
                    }
                } else {
                    inline for (0..new_size) |i| {
                        inline for (0..i + 1) |j| {
                            mat.data[i * new_size + j] = self.data[(start + i) * size + (start + j)];
                        }
                    }
                }
            }

            return mat;
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Static(size, N, uplo, layout)`): A
        ///   pointer to the matrix to get the submatrix from.
        /// * `start` (`usize`): The starting diagonal index of the submatrix
        ///   (inclusive).
        /// * `end` (`usize`): The ending diagonal index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The submatrix view.
        ///
        /// ## Errors
        /// * `matrix.Error.InvalidRange`: If the specified range is invalid.
        pub fn submatrixView(self: *matrix.symmetric.Static(size, N, uplo, layout), start: usize, end: usize) !matrix.symmetric.Dense(N, uplo, layout) {
            if (start >= size or end > size or start >= end)
                return matrix.Error.InvalidRange;

            return .{
                .data = @as([*]N, &self.data) + _index(start, start),
                .rows = end - start,
                .cols = end - start,
                .ld = size,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn _index(r: usize, c: usize) usize {
            return if (comptime layout == .col_major)
                r + c * size
            else
                r * size + c;
        }

        pub fn format(self: matrix.symmetric.Static(size, N, uplo, layout), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const matrix.symmetric.Static(size, N, uplo, layout), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .matrix = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                matrix: *const matrix.symmetric.Static(size, N, uplo, layout),

                pub fn format(self: matrix.symmetric.Static(size, N, uplo, layout).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    try writer.print("zsl.matrix.symmetric.Static({d}, {s}, .{s}, .{s}) ({d} × {d}):\n\n", .{ size, @typeName(N), @tagName(uplo), @tagName(layout), rows, cols });

                    return matutils.format(self, num_fmt, rows, cols, writer);
                }
            };
        }
    };
}
