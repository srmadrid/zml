const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

/// Static triangular matrix type, represented as a contiguous array of
/// `rows × cols` elements of type `N`, depending on `uplo`, stored in either
/// column-major or row-major order. Only the upper or lower triangular part of
/// the matrix is accessed, depending on the `uplo` parameter, and the diagonal
/// can be either unit, meaning all diagonal elements are assumed to be 1 and
/// not accessed, or non-unit, meaning the diagonal elements are accessed
/// normally.
pub fn Static(rows_: comptime_int, cols_: comptime_int, N: type, uplo: matrix.Uplo, diag: matrix.Diag, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.triangular.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [rows * cols]N,

        // Type signatures
        pub const is_matrix = true;
        pub const is_static = true;
        pub const is_triangular = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = uplo;
        pub const storage_diag: ?matrix.Diag = diag;

        pub const rows = rows_;
        pub const cols = cols_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.triangular.Static(rows, cols, N, uplo, diag, layout) = .{
            .data = undefined,
        };

        pub const init = empty;

        /// Initializes a new
        /// `matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`, with
        /// all elements set to the specified value.
        ///
        /// ## Arguments
        /// * `value` (`N`): The value to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`: The
        /// newly initialized matrix.
        pub fn initValue(value: N) matrix.triangular.Static(rows, cols, N, uplo, diag, layout) {
            var mat: matrix.triangular.Static(rows, cols, N, uplo, diag, layout) = .init;

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        inline for (0..cols) |j| {
                            inline for (0..(comptime int.min(j, rows))) |i| {
                                mat.data[i + j * rows] = value;
                            }
                        }
                    } else { // cun
                        inline for (0..cols) |j| {
                            inline for (0..(comptime int.min(j + 1, rows))) |i| {
                                mat.data[i + j * rows] = value;
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // clu
                        inline for (0..(comptime int.min(rows, cols))) |j| {
                            inline for ((comptime int.min(j + 1, rows))..rows) |i| {
                                mat.data[i + j * rows] = value;
                            }
                        }
                    } else { // cln
                        inline for (0..(comptime int.min(rows, cols))) |j| {
                            inline for ((comptime int.min(j, rows))..rows) |i| {
                                mat.data[i + j * rows] = value;
                            }
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // ruu
                        inline for (0..(comptime int.min(rows, cols))) |i| {
                            inline for ((comptime int.min(i + 1, cols))..cols) |j| {
                                mat.data[i * cols + j] = value;
                            }
                        }
                    } else { // run
                        inline for (0..(comptime int.min(rows, cols))) |i| {
                            inline for ((comptime int.min(i, cols))..cols) |j| {
                                mat.data[i * cols + j] = value;
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // rlu
                        inline for (0..rows) |i| {
                            inline for (0..(comptime int.min(i, cols))) |j| {
                                mat.data[i * cols + j] = value;
                            }
                        }
                    } else { // rln
                        inline for (0..rows) |i| {
                            inline for (0..(comptime int.min(i + 1, cols))) |j| {
                                mat.data[i * cols + j] = value;
                            }
                        }
                    }
                }
            }

            return mat;
        }

        /// Initializes a new
        /// `matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`, with
        /// all elements set from the specified values.
        ///
        /// ## Arguments
        /// * `values` (`[rows * cols]N`): The values to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`: The
        /// newly initialized matrix.
        pub fn initArray(values: [rows * cols]N) matrix.triangular.Static(rows, cols, N, uplo, diag, layout) {
            return .{ .data = values };
        }

        /// Initializes a new
        /// `matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`, with
        /// all elements in the stored triangular part set by calling the
        /// specified function with the given arguments.
        ///
        /// ## Arguments
        /// * `@"fn"` (`anytype`): The function to call to fill the matrix.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`: The
        /// newly initialized matrix.
        ///
        /// ## Errors
        pub fn initFn(comptime @"fn": anytype, args: anytype) !matrix.triangular.Static(rows, cols, N, uplo, diag, layout) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("matrix.triangular.Static(rows, cols, N, uplo, diag, layout).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.triangular.Static(rows, cols, N, uplo, diag, layout) = .init;

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        inline for (0..cols) |j| {
                            inline for (0..(comptime int.min(j, rows))) |i| {
                                mat.data[i + j * rows] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // cun
                        inline for (0..cols) |j| {
                            inline for (0..(comptime int.min(j + 1, rows))) |i| {
                                mat.data[i + j * rows] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // clu
                        inline for (0..(comptime int.min(rows, cols))) |j| {
                            inline for ((comptime int.min(j + 1, rows))..rows) |i| {
                                mat.data[i + j * rows] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // cln
                        inline for (0..(comptime int.min(rows, cols))) |j| {
                            inline for ((comptime int.min(j, rows))..rows) |i| {
                                mat.data[i + j * rows] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
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
                        inline for (0..(comptime int.min(rows, cols))) |i| {
                            inline for ((comptime int.min(i + 1, cols))..cols) |j| {
                                mat.data[i * cols + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // run
                        inline for (0..(comptime int.min(rows, cols))) |i| {
                            inline for ((comptime int.min(i, cols))..cols) |j| {
                                mat.data[i * cols + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // rlu
                        inline for (0..rows) |i| {
                            inline for (0..(comptime int.min(i, cols))) |j| {
                                mat.data[i * cols + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                                    try @call(.auto, @"fn", args)
                                else
                                    @call(.auto, @"fn", args);
                            }
                        }
                    } else { // rln
                        inline for (0..rows) |i| {
                            inline for (0..(comptime int.min(i + 1, cols))) |j| {
                                mat.data[i * cols + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
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

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Static(rows, cols, N, uplo, diag, order)`):
        ///   The matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.triangular.Static(rows, cols, N, uplo, diag, layout), r: usize, c: usize) !N {
            if (r >= rows or c >= cols)
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

            return self.data[_index(r, c)];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`):
        ///   The matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds, on the correct triangular part, and outside the
        ///   diagonal if `diag` is unit.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds, on the correct triangular part, and outside the
        ///   diagonal if `diag` is unit.
        ///
        /// ## Returns
        /// `N`: The element at the specified position.
        pub fn getAssumeInBounds(self: matrix.triangular.Static(rows, cols, N, uplo, diag, layout), r: usize, c: usize) N {
            return self.data[_index(r, c)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`):
        ///   A pointer to the matrix to set the element in.
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
        pub fn set(self: *matrix.triangular.Static(rows, cols, N, uplo, diag, layout), r: usize, c: usize, value: N) !void {
            if (r >= rows or c >= cols)
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

            self.data[_index(r, c)] = value;
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`):
        ///   A pointer to the matrix to set the element in.
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
        pub fn setAssumeInBounds(self: *matrix.triangular.Static(rows, cols, N, uplo, diag, layout), r: usize, c: usize, value: N) void {
            self.data[_index(r, c)] = value;
        }

        /// Returns the transposed copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`):
        ///   The matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.triangular.Static(cols, rows, N, uplo.invert(), diag, layout.invert())`:
        /// The transposed matrix.
        pub fn transposeCopy(self: matrix.triangular.Static(rows, cols, N, uplo, diag, layout)) matrix.triangular.Static(cols, rows, N, uplo.invert(), diag, layout.invert()) {
            return .{ .data = self.data };
        }

        /// Returns the transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`):
        ///   The matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo.invert(), diag, layout.invert())`:
        /// The transposed matrix.
        pub fn transposeView(self: matrix.triangular.Static(rows, cols, N, uplo, diag, layout)) matrix.triangular.Dense(N, uplo.invert(), diag, layout.invert()) {
            return .{
                .data = &self.data,
                .rows = cols,
                .cols = rows,
                .ld = if (comptime layout == .col_major) rows else cols,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a submatrix copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`):
        ///   The matrix to get the submatrix from.
        /// * `start` (`comptime usize`): The starting diagonal index of the submatrix
        ///   (inclusive).
        /// * `row_end` (`comptime usize`): The ending row index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        /// * `col_end` (`comptime usize`): The ending column index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.triangular.Static(row_end - start, col_end - start, N, uplo, diag, layout)`:
        /// The copied submatrix.
        pub fn submatrixCopy(
            self: matrix.triangular.Static(rows, cols, N, uplo, diag, layout),
            comptime start: usize,
            comptime row_end: usize,
            comptime col_end: usize,
        ) matrix.triangular.Static(row_end - start, col_end - start, N, uplo, diag, layout) {
            comptime if (start >= int.min(rows, cols) or
                row_end > rows or col_end > cols or
                row_end <= start or col_end <= start)
                @compileError("zsl.matrix.triangular.Static.submatrixCopy: Invalid range");

            const new_rows = row_end - start;
            const new_cols = col_end - start;
            var mat: matrix.triangular.Static(new_rows, new_cols, N, uplo, diag, layout) = .init;

            if (comptime layout == .col_major) {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // cuu
                        inline for (0..new_cols) |j| {
                            inline for (0..(comptime int.min(j, new_rows))) |i| {
                                mat.data[i + j * new_rows] = self.data[(start + i) + (start + j) * rows];
                            }
                        }
                    } else { // cun
                        inline for (0..new_cols) |j| {
                            inline for (0..(comptime int.min(j + 1, new_rows))) |i| {
                                mat.data[i + j * new_rows] = self.data[(start + i) + (start + j) * rows];
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // clu
                        inline for (0..(comptime int.min(new_rows, new_cols))) |j| {
                            inline for ((comptime int.min(j + 1, new_rows))..new_rows) |i| {
                                mat.data[i + j * new_rows] = self.data[(start + i) + (start + j) * rows];
                            }
                        }
                    } else { // cln
                        inline for (0..(comptime int.min(new_rows, new_cols))) |j| {
                            inline for ((comptime int.min(j, new_rows))..new_rows) |i| {
                                mat.data[i + j * new_rows] = self.data[(start + i) + (start + j) * rows];
                            }
                        }
                    }
                }
            } else {
                if (comptime uplo == .upper) {
                    if (comptime diag == .unit) { // ruu
                        inline for (0..(comptime int.min(new_rows, new_cols))) |i| {
                            inline for ((comptime int.min(i + 1, new_cols))..new_cols) |j| {
                                mat.data[i * new_cols + j] = self.data[(start + i) * cols + (start + j)];
                            }
                        }
                    } else { // run
                        inline for (0..(comptime int.min(new_rows, new_cols))) |i| {
                            inline for ((comptime int.min(i, new_cols))..new_cols) |j| {
                                mat.data[i * new_cols + j] = self.data[(start + i) * cols + (start + j)];
                            }
                        }
                    }
                } else {
                    if (comptime diag == .unit) { // rlu
                        inline for (0..new_rows) |i| {
                            inline for (0..(comptime int.min(i, new_cols))) |j| {
                                mat.data[i * new_cols + j] = self.data[(start + i) * cols + (start + j)];
                            }
                        }
                    } else { // rln
                        inline for (0..new_rows) |i| {
                            inline for (0..(comptime int.min(i + 1, new_cols))) |j| {
                                mat.data[i * new_cols + j] = self.data[(start + i) * cols + (start + j)];
                            }
                        }
                    }
                }
            }

            return mat;
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.triangular.Static(rows, cols, N, uplo, diag, layout)`):
        ///   The matrix to get the submatrix from.
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
        pub fn submatrixView(self: matrix.triangular.Static(rows, cols, N, uplo, diag, layout), start: usize, row_end: usize, col_end: usize) !matrix.triangular.Dense(N, uplo, diag, layout) {
            if (start >= int.min(rows, cols) or
                row_end > rows or col_end > cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            return .{
                .data = @as([*]N, &self.data) + _index(start, start),
                .rows = row_end - start,
                .cols = col_end - start,
                .ld = if (comptime layout == .col_major) rows else cols,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn _index(r: usize, c: usize) usize {
            return if (comptime layout == .col_major)
                r + c * rows
            else
                r * cols + c;
        }
    };
}
