const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");

/// Static general matrix type, represented as a contiguous array of
/// `rows × cols` elements of type `N`, stored in either column-major or
/// row-major order.
pub fn Static(rows_: comptime_int, cols_: comptime_int, N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.general.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [rows * cols]N,

        // Type signatures
        pub const is_matrix = true;
        pub const is_static = true;
        pub const is_general = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = null;
        pub const storage_diag: ?matrix.Diag = null;
        pub const rows = rows_;
        pub const cols = cols_;

        // Numeric type
        pub const Numeric = N;

        pub const empty = matrix.general.Static(rows, cols, N, layout){
            .data = undefined,
        };

        pub const init = empty;

        /// Initializes a new `matrix.general.Static(rows, cols, N, layout)`,
        /// with all elements set to the specified value.
        ///
        /// ## Arguments
        /// * `value` (`N`): The value to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.general.Static(rows, cols, N, layout)`: The newly
        /// initialized matrix.
        pub fn initValue(value: N) matrix.general.Static(rows, cols, N, layout) {
            return .{ .data = @splat(value) };
        }

        /// Initializes a new `matrix.general.Static(rows, cols, N, layout)`,
        /// with all elements set from the specified values.
        ///
        /// ## Arguments
        /// * `values` (`[rows * cols]N`): The values to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.general.Static(rows, cols, N, layout)`: The newly
        /// initialized matrix.
        pub fn initArray(values: [rows * cols]N) matrix.general.Static(rows, cols, N, layout) {
            return .{ .data = values };
        }

        /// Initializes a new `matrix.general.Static(rows, cols,N, layout)`,
        /// with all elements set by calling the specified function with the
        /// given arguments.
        ///
        /// ## Arguments
        /// * `@"fn"` (`anytype`): The function to call to fill the matrix.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `matrix.general.Static(rows, cols,N, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        pub fn initFn(comptime @"fn": anytype, args: anytype) !matrix.general.Static(rows, cols, N, layout) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.general.matrix.general.Static(rows, cols,N, layout).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.general.Static(rows, cols, N, layout) = .init;

            if (comptime layout == .col_major) {
                inline for (0..cols) |j| {
                    inline for (0..rows) |i| {
                        mat.data[i + j * rows] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                            try @call(.auto, @"fn", args)
                        else
                            @call(.auto, @"fn", args);
                    }
                }
            } else {
                inline for (0..rows) |i| {
                    inline for (0..cols) |j| {
                        mat.data[i * cols + j] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                            try @call(.auto, @"fn", args)
                        else
                            @call(.auto, @"fn", args);
                    }
                }
            }

            return mat;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
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
        pub fn get(self: matrix.general.Static(rows, cols, N, layout), r: usize, c: usize) !N {
            if (r >= rows or c >= cols)
                return matrix.Error.PositionOutOfBounds;

            return self.data[_index(r, c)];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: matrix.general.Static(rows, cols, N, layout), r: usize, c: usize) N {
            return self.data[_index(r, c)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)`): A
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
        pub fn set(self: *matrix.general.Static(rows, cols, N, layout), r: usize, c: usize, value: N) !void {
            if (r >= rows or c >= cols)
                return matrix.Error.PositionOutOfBounds;

            self.data[_index(r, c)] = value;
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)`): A
        ///   pointer to the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.general.Static(rows, cols, N, layout), r: usize, c: usize, value: N) void {
            self.data[_index(r, c)] = value;
        }

        /// Returns the transposed copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.general.Static(cols, rows, N, layout.invert())`: The
        /// transposed matrix.
        pub fn transposeCopy(self: matrix.general.Static(rows, cols, N, layout)) matrix.general.Static(cols, rows, N, layout.invert()) {
            return .{ .data = self.data };
        }

        /// Returns the transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout.invert())`: The transposed matrix.
        pub fn transposeView(self: matrix.general.Static(rows, cols, N, layout)) matrix.general.Dense(N, layout.invert()) {
            return .{
                .data = &self.data,
                .rows = cols,
                .cols = rows,
                .ld = if (comptime layout == .col_major) rows else cols,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a submatrix copy of the matrix
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The matrix
        ///   to get the submatrix from.
        /// * `row_start` (`comptime usize`): The starting row index (inclusive).
        /// * `row_end` (`comptime usize`): The ending row index (exclusive).
        /// * `col_start` (`comptime usize`): The starting column index (inclusive).
        /// * `col_end` (`comptime usize`): The ending column index (exclusive).
        ///
        /// ## Returns
        /// `matrix.general.Static(row_end - row_start, col_end - col_start, N, layout)`:
        /// The copied submatrix.
        pub fn submatrixCopy(
            self: matrix.general.Static(rows, cols, N, layout),
            comptime row_start: usize,
            comptime row_end: usize,
            comptime col_start: usize,
            comptime col_end: usize,
        ) matrix.general.Static(row_end - row_start, col_end - col_start, N, layout) {
            comptime {
                if (row_start >= rows or col_start >= cols or
                    row_end > rows or col_end > cols or
                    row_start >= row_end or col_start >= col_end)
                    @compileError("zsl.matrix.general.Static(rows, cols, N, layout).submatrixCopy: Invalid range");
            }

            const new_rows = row_end - row_start;
            const new_cols = col_end - col_start;
            var mat: matrix.general.Static(new_rows, new_cols, N, layout) = .init;

            if (comptime layout == .col_major) {
                inline for (0..new_cols) |j| {
                    inline for (0..new_rows) |i| {
                        mat.data[i + j * new_rows] = self.data[row_start + i + (col_start + j) * rows];
                    }
                }
            } else {
                inline for (0..new_rows) |i| {
                    inline for (0..new_cols) |j| {
                        mat.data[i * new_cols + j] = self.data[(row_start + i) * cols + col_start + j];
                    }
                }
            }

            return mat;
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)``): A
        ///   pointer to the matrix to get the submatrix from.
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
        pub fn submatrixView(self: *matrix.general.Static(rows, cols, N, layout), row_start: usize, row_end: usize, col_start: usize, col_end: usize) !matrix.general.Dense(N, layout) {
            if (row_start >= rows or col_start >= cols or
                row_end > rows or col_end > cols or
                row_start >= row_end or col_start >= col_end)
                return matrix.Error.InvalidRange;

            return .{
                .data = @as([*]N, &self.data) + _index(row_start, col_start),
                .rows = row_end - row_start,
                .cols = col_end - col_start,
                .ld = if (comptime layout == .col_major) rows else cols,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a copy of the specified row as a static vector.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to get the row from.
        /// * `r` (`usize`): The row index to get.
        ///
        /// ## Returns
        /// `vector.Static(cols, N)`: The specified row as a static vector.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` is out of bounds.
        pub fn rowCopy(self: matrix.general.Static(rows, cols, N, layout), r: usize) !vector.Static(cols, N) {
            if (r >= rows)
                return matrix.Error.PositionOutOfBounds;

            var result: vector.Static(cols, N) = .init;

            inline for (0..cols) |j| {
                result.data[j] = self.data[_index(r, j)];
            }

            return result;
        }

        /// Returns a view of the specified row as a dense vector.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)`): A
        ///   pointer to the matrix to get the row from.
        /// * `r` (`usize`): The row index to get.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The specified row as a dense vector.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` is out of bounds.
        pub fn rowView(self: *matrix.general.Static(rows, cols, N, layout), r: usize) !vector.Dense(N) {
            if (r >= rows)
                return matrix.Error.PositionOutOfBounds;

            return .{
                .data = @as([*]N, &self.data) + _index(r, 0),
                .len = cols,
                .inc = if (comptime layout == .col_major)
                    numeric.cast(isize, rows)
                else
                    1,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a copy of the specified column as a static vector.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to get the column from.
        /// * `c` (`usize`): The column index to get.
        ///
        /// ## Returns
        /// `vector.Static(rows, N)`: The specified column as a static vector.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `c` is out of bounds.
        pub fn colCopy(self: matrix.general.Static(rows, cols, N, layout), c: usize) !vector.Static(rows, N) {
            if (c >= cols)
                return matrix.Error.PositionOutOfBounds;

            var result: vector.Static(rows, N) = .init;

            inline for (0..rows) |i| {
                result.data[i] = self.data[_index(i, c)];
            }

            return result;
        }

        /// Returns a view of the specified column as a dense vector.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)`): A
        ///   pointer to the matrix to get the column from.
        /// * `c` (`usize`): The column index to get.
        ///
        /// ## Returns
        /// `vector.Static(rows, N)`: The specified column as a dense vector.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `c` is out of bounds.
        pub fn colView(self: *matrix.general.Static(rows, cols, N, layout), c: usize) !vector.Dense(N) {
            if (c >= cols)
                return matrix.Error.PositionOutOfBounds;

            return .{
                .data = @as([*]N, &self.data) + _index(0, c),
                .len = rows,
                .inc = if (comptime layout == .col_major)
                    1
                else
                    numeric.cast(isize, cols),
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a copy of the matrix as a symmetric static matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to copy.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        ///
        /// ## Returns
        /// `matrix.symmetric.Static(rows, N, uplo, layout)`: The symmetric
        /// static matrix copy.
        pub fn symmetricCopy(self: matrix.general.Static(rows, cols, N, layout), comptime uplo: matrix.Uplo) !matrix.symmetric.Static(rows, cols, N, uplo, layout) {
            if (comptime rows != cols)
                @compileError("zsl.matrix.general.Static(rows, cols, N, layout).symmetricCopy: Matrix must be square");

            return .{ .data = self.data };
        }

        /// Returns a view of the matrix as a symmetric dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)`): A
        ///   pointer to the matrix to get the view of.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        ///
        /// ## Returns
        /// `matrix.symmetric.Dense(N, uplo, layout)`: The symmetric dense
        /// matrix view.
        pub fn symmetricView(self: *matrix.general.Static(rows, cols, N, layout), comptime uplo: matrix.Uplo) !matrix.symmetric.Dense(N, uplo, layout) {
            if (comptime rows != cols)
                @compileError("zsl.matrix.general.Static(rows, cols, N, layout).symmetricView: Matrix must be square");

            return .{
                .data = &self.data,
                .rows = rows,
                .cols = cols,
                .ld = if (comptime layout == .col_major) rows else cols,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a copy of the matrix as a Hermitian static matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to copy.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        ///
        /// ## Returns
        /// `matrix.hermitian.Static(rows, N, uplo, layout)`: The Hermitian
        /// static matrix copy.
        pub fn hermitianCopy(self: matrix.general.Static(rows, cols, N, layout), comptime uplo: matrix.Uplo) !matrix.hermitian.Static(rows, N, uplo, layout) {
            if (comptime rows != cols)
                @compileError("zsl.matrix.general.Static(rows, cols, N, layout).hermitianCopy: Matrix must be square");

            return .{ .data = self.data };
        }

        /// Returns a view of the matrix as a Hermitian dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)`): A
        ///   pointer to the matrix to get the view of.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        ///
        /// ## Returns
        /// `matrix.hermitian.Dense(N, uplo, layout)`: The Hermitian dense
        /// matrix view.
        pub fn hermitianView(self: *matrix.general.Static(rows, cols, N, layout), comptime uplo: matrix.Uplo) !matrix.hermitian.Dense(N, uplo, layout) {
            if (comptime rows != cols)
                @compileError("zsl.matrix.general.Static(rows, cols, N, layout).hermitianView: Matrix must be square");

            return .{
                .data = &self.data,
                .rows = rows,
                .cols = cols,
                .ld = if (comptime layout == .col_major) rows else cols,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a copy of the matrix as a triangular static matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Static(rows, cols, N, layout)`): The
        ///   matrix to copy.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        /// * `diag` (`comptime matrix.Diag`): Specifies whether the matrix is
        ///   unit triangular (diagonal elements are assumed to be 1 and are
        ///   ignored) or non-unit triangular.
        ///
        /// ## Returns
        /// `matrix.triangular.Static(rows, N, uplo, diag, layout)`: The
        /// triangular static matrix copy.
        pub fn triangularCopy(self: matrix.general.Static(rows, cols, N, layout), comptime uplo: matrix.Uplo, comptime diag: matrix.Diag) matrix.triangular.Static(rows, cols, N, uplo, diag, layout) {
            return .{ .data = self.data };
        }

        /// Returns a view of the matrix as a triangular dense matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Static(rows, cols, N, layout)`): A
        ///   pointer to the matrix to get the view of.
        /// * `uplo` (`comptime matrix.Uplo`): Specifies whether the upper or
        ///   lower triangle of the matrix is used, the other triangle is
        ///   ignored.
        /// * `diag` (`comptime matrix.Diag`): Specifies whether the matrix is
        ///   unit triangular (diagonal elements are assumed to be 1 and are
        ///   ignored) or non-unit triangular.
        ///
        /// ## Returns
        /// `matrix.triangular.Dense(N, uplo, diag, layout)`: The truangular
        /// dense matrix view.
        pub fn triangularView(self: *matrix.general.Static(rows, cols, N, layout), comptime uplo: matrix.Uplo, comptime diag: matrix.Diag) matrix.triangular.Dense(N, uplo, diag, layout) {
            return .{
                .data = &self.data,
                .rows = rows,
                .cols = cols,
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
