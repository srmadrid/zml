const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

const matutils = @import("../utils.zig");

/// Static diagonal matrix type, represented as a contiguous array of
/// `min(rows, cols)` elements of type `N`.
pub fn Static(rows_: comptime_int, cols_: comptime_int, N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.diagonal.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [int.min(rows, cols)]N,

        // Type signatures
        pub const is_matrix = true;
        pub const is_static = true;
        pub const is_diagonal = true;
        pub const storage_layout: ?matrix.Layout = null;
        pub const storage_uplo: ?matrix.Uplo = null;
        pub const storage_diag: ?matrix.Diag = null;
        pub const rows = rows_;
        pub const cols = cols_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.diagonal.Static(rows, cols, N) = .{
            .data = undefined,
        };

        pub const init = empty;

        /// Initializes a new `matrix.diagonal.Static(rows, cols, N)` with the specified
        /// rows and columns, with all diagonal elements set to the specified
        /// value.
        ///
        /// ## Arguments
        /// * `value` (`N`): The value to fill the matrix with.
        ///
        /// ## Returns
        /// `matrix.diagonal.Static(rows, cols, N)`: The newly initialized matrix.
        pub fn initValue(value: N) matrix.diagonal.Static(rows, cols, N) {
            const mat: matrix.diagonal.Static(N) = .init;

            inline for (0..(comptime int.min(rows, cols))) |i| {
                mat.data[i] = value;
            }

            return mat;
        }

        /// Initializes a new `matrix.diagonal.Static(rows, cols, N)` with the
        /// specified rows and columns, with all diagonal elements set by
        /// calling the specified function with the given arguments.
        ///
        /// ## Arguments
        /// * `@"fn"` (`anytype`): The function to call to fill the matrix.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `matrix.diagonal.Static(rows, cols, N)`: The newly initialized
        /// matrix.
        pub fn initFn(comptime @"fn": anytype, args: anytype) !matrix.diagonal.Static(rows, cols, N) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.diagonal.Static(rows, cols, N).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.diagonal.Static(rows, cols, N) = .init;

            inline for (0..(comptime int.min(rows, cols))) |i| {
                mat.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                    try @call(.auto, @"fn", args)
                else
                    @call(.auto, @"fn", args);
            }

            return mat;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Static(rows, cols, N)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified position.
        ///
        /// ## Errors
        /// `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of bounds.
        pub fn get(self: matrix.diagonal.Static(rows, cols, N), r: usize, c: usize) !N {
            if (r >= rows or c >= cols)
                return matrix.Error.PositionOutOfBounds;

            if (r != c)
                return numeric.cast(N, 0);

            return self.data[r];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Static(rows, cols, N)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds and equal to `c`.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds and equal to `r`.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: matrix.diagonal.Static(rows, cols, N), r: usize, c: usize) N {
            _ = c;
            return self.data[r];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.diagonal.Static(rows, cols, N)`): A pointer to
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
        /// * `matrix.Error.BreaksStructure`: If `r` is not equal to `c`.
        pub fn set(self: *matrix.diagonal.Static(rows, cols, N), r: usize, c: usize, value: N) !void {
            if (r >= rows or c >= cols)
                return matrix.Error.PositionOutOfBounds;

            if (r != c)
                return matrix.Error.BreaksStructure;

            self.data[r] = value;
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.diagonal.Static(rows, cols, N)`): A pointer to
        ///   the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds and equal to `c`.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds and equal to `r`.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.diagonal.Static(rows, cols, N), r: usize, c: usize, value: N) void {
            _ = c;
            self.data[r] = value;
        }

        /// Returns a submatrix copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Static(rows, cols, N)`): The matrix to
        ///   get the submatrix from.
        /// * `start` (`comptime usize`): The starting diagonal index of the
        ///   submatrix (inclusive).
        /// * `row_end` (`comptime usize`): The ending row index of the
        ///   submatrix (exclusive). Must be greater than `start`.
        /// * `col_end` (`comptime usize`): The ending column index of the
        ///   submatrix (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.diagonal.Static(row_end - start, col_end - start, N)`: The submatrix.
        pub fn submatrixCopy(self: matrix.diagonal.Static(rows, cols, N), comptime start: usize, comptime row_end: usize, comptime col_end: usize) matrix.diagonal.Static(row_end - start, col_end - start, N) {
            comptime if (start >= int.min(rows, cols) or
                row_end > rows or col_end > cols or
                row_end < start or col_end < start)
                @compileError("zsl.matrix.diagonal.Static(rows, cols, N).submatrixCopy: Invalid range");

            const new_rows = row_end - start;
            const new_cols = col_end - start;
            var mat: matrix.diagonal.Static(new_rows, new_cols, N) = .init;

            inline for (0..(comptime int.min(new_rows, new_cols))) |i| {
                mat.data[i] = self.data[start + i];
            }

            return mat;
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Static(rows, cols, N)`): The matrix to
        ///   get the submatrix from.
        /// * `start` (`usize`): The starting diagonal index of the submatrix
        ///   (inclusive).
        /// * `row_end` (`usize`): The ending row index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        /// * `col_end` (`usize`): The ending column index of the submatrix
        ///   (exclusive). Must be greater than `start`.
        ///
        /// ## Returns
        /// `matrix.diagonal.Sparse(N)`: The submatrix.
        ///
        /// ## Errors
        /// * `matrix.Error.InvalidRange`: If the specified range is invalid.
        pub fn submatrixView(self: matrix.diagonal.Static(rows, cols, N), start: usize, row_end: usize, col_end: usize) !matrix.diagonal.Sparse(N) {
            if (start >= int.min(rows, cols) or
                row_end > rows or col_end > cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            return .{
                .data = self.data + start,
                .rows = row_end - start,
                .cols = col_end - start,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn format(self: matrix.diagonal.Static(rows, cols, N), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const matrix.diagonal.Static(rows, cols, N), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .matrix = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                matrix: *const matrix.diagonal.Static(rows, cols, N),

                pub fn format(self: matrix.diagonal.Static(rows, cols, N).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    try writer.print("zsl.matrix.diagonal.Static({d}, {d}, {s}) ({d} × {d}):\n\n", .{ rows, cols, @typeName(N), rows, cols });

                    return matutils.format(self, num_fmt, rows, cols, writer);
                }
            };
        }
    };
}
