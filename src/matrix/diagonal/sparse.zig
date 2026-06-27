const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

/// Sparse diagonal matrix type, represented as a contiguous array of
/// `min(rows, cols)` elements of type `N`.
pub fn Sparse(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.diagonal.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        rows: usize,
        cols: usize,
        flags: matrix.Flags = .{},

        // Type signatures
        pub const is_matrix = true;
        pub const is_sparse = true;
        pub const is_diagonal = true;
        pub const storage_layout: ?matrix.Layout = null;
        pub const storage_uplo: ?matrix.Uplo = null;
        pub const storage_diag: ?matrix.Diag = null;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.diagonal.Sparse(N) = .{
            .data = &.{},
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.diagonal.Sparse(N)` with the specified
        /// rows and columns.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the matrix.
        /// * `cols` (`usize`): The columns of the matrix.
        ///
        /// ## Returns
        /// `matrix.diagonal.Sparse(N)`: The newly initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !matrix.diagonal.Sparse(N) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .data = (try allocator.alloc(N, int.min(rows, cols))).ptr,
                .rows = rows,
                .cols = cols,
                .flags = .{ .owns_data = true },
            };
        }

        // pub fn initBuffer

        /// Initializes a new `matrix.diagonal.Sparse(N)` with the specified
        /// rows and columns, with all diagonal elements set to the specified
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
        /// `matrix.diagonal.Sparse(N)`: The newly initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn initValue(allocator: std.mem.Allocator, rows: usize, cols: usize, value: N) !matrix.diagonal.Sparse(N) {
            const mat: matrix.diagonal.Sparse(N) = try .init(allocator, rows, cols);

            var i: usize = 0;
            while (i < int.min(rows, cols)) : (i += 1) {
                mat.data[i] = value;
            }

            return mat;
        }

        /// Initializes a new `matrix.diagonal.Sparse(N)` with the specified
        /// rows and columns, with all diagonal elements set by calling the
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
        /// `matrix.diagonal.Sparse(N)`: The newly initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        pub fn initFn(allocator: std.mem.Allocator, rows: usize, cols: usize, comptime @"fn": anytype, args: anytype) !matrix.diagonal.Sparse(N) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.matrix.diagonal.Sparse(N).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var mat: matrix.diagonal.Sparse(N) = try .init(allocator, rows, cols);
            errdefer mat.deinit(allocator);

            var i: usize = 0;
            while (i < int.min(rows, cols)) : (i += 1) {
                mat.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                    try @call(.auto, @"fn", args)
                else
                    @call(.auto, @"fn", args);
            }

            return mat;
        }

        /// Initializes a new identity `matrix.diagonal.Sparse(N)` of the
        /// specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.diagonal.Sparse(N)`: The newly initialized identity matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initIdentity(allocator: std.mem.Allocator, size: usize) !matrix.diagonal.Sparse(N) {
            const mat: matrix.diagonal.Sparse(N) = try .init(allocator, size, size);

            var i: usize = 0;
            while (i < size) : (i += 1) {
                mat.data[i] = numeric.one(N);
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.diagonal.Sparse(N)`): A pointer to the matrix to
        ///   deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.diagonal.Sparse(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..(int.min(self.rows, self.cols))]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Sparse(N)`): The matrix to get the
        ///   element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified position.
        ///
        /// ## Errors
        /// `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of bounds.
        pub fn get(self: matrix.diagonal.Sparse(N), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (r != c)
                return numeric.zero(N);

            return if (self.flags.noconj)
                self.data[r]
            else
                numeric.conj(self.data[r]);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Sparse(N)`): The matrix to get the
        ///   element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds and equal to `c`.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds and equal to `r`.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: matrix.diagonal.Sparse(N), r: usize, c: usize) N {
            _ = c;
            return if (self.flags.noconj)
                self.data[r]
            else
                numeric.conj(self.data[r]);
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.diagonal.Sparse(N)`): A pointer to the matrix to
        ///   set the
        ///   element in.
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
        pub fn set(self: *matrix.diagonal.Sparse(N), r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (r != c)
                return matrix.Error.BreaksStructure;

            self.data[r] = if (self.flags.noconj) value else numeric.conj(value);
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.diagonal.Sparse(N)`): A pointer to the matrix to
        ///   set the element in.
        /// * `r` (`usize`): The row index of the element to set. Assumed to be
        ///   within bounds and equal to `c`.
        /// * `c` (`usize`): The column index of the element to set. Assumed to
        ///   be within bounds and equal to `r`.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *matrix.diagonal.Sparse(N), r: usize, c: usize, value: N) void {
            _ = c;
            self.data[r] = if (self.flags.noconj) value else numeric.conj(value);
        }

        /// Creates a copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Sparse(N)`): The matrix to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.diagonal.Sparse(N)`: The copied matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn clone(self: matrix.diagonal.Sparse(N), allocator: std.mem.Allocator) !matrix.diagonal.Sparse(N) {
            const mat: matrix.diagonal.Sparse(N) = try .init(allocator, self.rows, self.cols);

            var i: usize = 0;
            while (i < int.min(mat.rows, mat.cols)) : (i += 1) {
                mat.data[i] = if (self.flags.noconj) self.data[i] else numeric.conj(self.data[i]);
            }

            return mat;
        }

        /// Returns a submatrix view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.diagonal.Sparse(N)`): The matrix to get the
        ///   submatrix from.
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
        pub fn submatrixView(self: matrix.diagonal.Sparse(N), start: usize, row_end: usize, col_end: usize) !matrix.diagonal.Sparse(N) {
            if (start >= int.min(self.rows, self.cols) or
                row_end > self.rows or col_end > self.cols or
                row_end < start or col_end < start)
                return matrix.Error.InvalidRange;

            return .{
                .data = self.data + start,
                .rows = row_end - start,
                .cols = col_end - start,
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }
    };
}
